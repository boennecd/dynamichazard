#include "dynamichazard.h"

extern "C"
{
  /*
  *
  * The Fortran code is:
  *       SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D,
                            *+      RBAR, THETAB, SSERR, IER)
                            * ...
                            *
                            *C**** WARNING: The elements of XROW are overwritten  ****
                            *C
                            *      INTEGER NP, NRBAR, IER
                            *      DOUBLE PRECISION WEIGHT, XROW(NP), YELEM, D(NP), RBAR(*),
                            *      +    THETAB(NP), SSERR
                            * ...
                            */
  void includ_(int *, int *, double *, double *, double *,
               double *, double *, double *, double *, int *);

  /*
  * The Fortran code is:
  *      SUBROUTINE REGCF(NP, NRBAR, D, RBAR, THETAB, TOL, BETA,
                          *      +     NREQ, IER)
                          * ...
                          *      INTEGER NP, NRBAR, NREQ, IER
                          *      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), TOL(NP),
                          *      +     BETA(NP)
                          * ...
                          */
  void regcf_(int *, int *, double *, double *, double *,
              double *, double *, int *, int *);
}


int binomialCoeff(int n, int k)
{
  // Base Cases
  if (k==0 || k==n)
    return 1;

  // Recur
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}



template<class T>
arma::vec bigglm_updateQR<T>::linkinv(const arma::vec &eta){
  arma::vec out(eta.n_elem);

  for(arma::uword i = 0; i < eta.n_elem; i++)
    out[i] = t.link_func_inv(eta[i]);

  return(out);
}

template<class T>
arma::vec bigglm_updateQR<T>::d_mu_d_eta(const arma::vec &eta){
  arma::vec out(eta.n_elem);

  for(arma::uword i = 0; i < eta.n_elem; i++)
    out[i] = t.d_mu_d_eta(eta[i]);

  return(out);
}

template<class T>
arma::vec bigglm_updateQR<T>::variance(const arma::vec &mu){
  arma::vec out(mu.n_elem);

  for(arma::uword i = 0; i < mu.n_elem; i++)
    out[i] = t.variance(mu[i]);

  return(out);
}


template<class T>
void bigglm_updateQR<T>::update(qr_obj &qr, // Previous/starting value. Will be overwritten
            const arma::mat &X, const arma::vec &eta,
            const arma::vec &offset, arma::vec &y) // y will not be altered
{
  arma::vec eta_plus_off = eta + offset;
  arma::vec mu = linkinv(eta_plus_off);
  arma::vec dmu = d_mu_d_eta(eta_plus_off);
  arma::vec z = eta + (y - mu) / dmu; // note that offset is not added as in bigglm.function
  arma::vec ww = dmu % dmu / variance(mu);

  int n_parems = X.n_rows;
  int nrbar = qr.rbar->n_elem;
  int ier = 0;


  double *y_ptr = z.memptr(); // a bit confussion with the chance of notion! This is "the same y as in the C code"
  double *w_ptr = ww.memptr();

  for(unsigned int i = 0; i < ww.n_elem; ++i, ++y_ptr, ++w_ptr){
    // copy taken as subroutine will overwrite values so we take a copy
    arma::vec x_row = X.col(i);

    //Here is the "stack trace" from bigglm:
    //qr<-update(qr,mm,z-off,ww)
    //function (bigQR, X, y, w = NULL, singcheck = FALSE, add.intercept = FALSE)
    //bigQR = qr  X = mm  y = z-off w = ww
    //  .Call("updateQR", X, y, w, bigQR, add.intercept)
    //  updateQR(SEXP X, SEXP y, SEXP w, SEXP bigQR, SEXP intercept)
    //    F77_CALL(includ)(&p, &nrbar, REAL(w)+i, row,
    //    REAL(y)+i, REAL(D), REAL(Rbar), REAL(thetab),
    //    REAL(sse), &ier);

    includ_(&n_parems, &nrbar, w_ptr, x_row.memptr(),
            y_ptr, qr.D->memptr(), qr.rbar->memptr(), qr.thetab->memptr(),
            &qr.ss, &ier);
  }
}

arma::vec bigglm_regcf(qr_obj &qr){
  int p = qr.D->n_elem;
  arma::vec beta(p, arma::fill::ones);

  int dum_arg1 = (p * p / 2);
  int dum_arg2 = 1;

  regcf_(&p, &dum_arg1, qr.D->memptr(), qr.rbar->memptr(),
         qr.thetab->memptr(), qr.tol->memptr(), beta.memptr(),
         &p, &dum_arg2);

  return beta;
}


// Only exported for tests
// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset, arma::vec &y){
  qr_obj qr;
  qr.D = &D; // the latter make sure that the same memory is used
  qr.rbar = &rbar;
  qr.thetab = &thetab;
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = &tol;

  if(model == "logit"){
    return(bigglm_updateQR_logit().update(qr, X, eta, offset, y));
  } else if (model == "exponential_combined" ||
             model == "exponential_binary_only"  ||
             model == "exponential_trunc_time_only"){
    return(bigglm_updateQR_poisson().update(qr, X, eta, offset, y));
  }
}

// Only exported for tests
// [[Rcpp::export]]
arma::vec bigglm_regcf_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                            double &ss, bool &checked, arma::vec &tol){
  qr_obj qr;
  qr.D = &D;
  qr.rbar = &rbar;
  qr.thetab = &thetab;
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = &tol;

  return(bigglm_regcf(qr));
}
