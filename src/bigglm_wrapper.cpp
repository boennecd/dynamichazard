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
class bigglm_updateQR{
  // match logic update.bigqr
  arma::vec linkinv(const arma::vec &eta){
    arma::vec out(eta.n_elem);

    for(arma::uword i = 0; i < eta.n_elem; i++)
      out[i] = t.link_func_inv(eta[i]);

    return(out);
  }

  arma::vec d_mu_d_eta(const arma::vec &eta){
    arma::vec out(eta.n_elem);

    for(arma::uword i = 0; i < eta.n_elem; i++)
      out[i] = t.d_mu_d_eta(eta[i]);

    return(out);
  }

  arma::vec variance(const arma::vec &mu){
    arma::vec out(mu.n_elem);

    for(arma::uword i = 0; i < mu.n_elem; i++)
      out[i] = t.variance(mu[i]);

    return(out);
  }

protected:
  T t;

public:
  bigglm_updateQR<T>(): t() {}

  void update(qr_obj &qr, // Previous/starting value. Will be overwritten
              const arma::mat &X, const arma::vec &eta,
              const arma::vec &offset, arma::vec &y) // y will not be altered
  {
    //TODO: look into fortran memory storage versus C++ / armadillo

    arma::vec eta_plus_off = eta + offset;
    arma::vec mu = linkinv(eta_plus_off);
    arma::vec dmu = d_mu_d_eta(eta_plus_off);
    arma::vec z = eta + (y - mu) / dmu; // note that offset is not added as in bigglm.function
    arma::vec ww = dmu % dmu % variance(mu);

    int n_parems = X.n_rows;
    int nrbar = qr.rbar->n_elem;
    int ier = 0;


    auto y_ptr = y.memptr();
    auto z_ptr = z.memptr();
    for(int i = 0; i < ww.n_cols; ++i, ++y_ptr, ++z_ptr){
      // created as subroutine will overwrite values so we take a copy
      arma::vec x_row(X.colptr(i), X.n_rows);

      //Here is the stack trace from bigglm:
      //qr<-update(qr,mm,z-off,ww)
      //function (bigQR, X, y, w = NULL, singcheck = FALSE, add.intercept = FALSE)
      //  .Call("updateQR", X, y, w, bigQR, add.intercept)
      //  updateQR(SEXP X, SEXP y, SEXP w, SEXP bigQR, SEXP intercept)
      //    F77_CALL(includ)(&p, &nrbar, REAL(w)+i, row,
      //    REAL(y)+i, REAL(D), REAL(Rbar), REAL(thetab),
      //    REAL(sse), &ier);

      includ_(&n_parems, &nrbar, z_ptr, x_row.memptr(),
              y_ptr, qr.D->memptr(), qr.rbar->memptr(), qr.thetab->memptr(),
              &qr.ss, &ier);
    }
  }
};

using bigglm_updateQR_logit = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;


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
  qr.D = &D;
  qr.rbar = &rbar;
  qr.thetab = &thetab;
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = &tol;

  if(model == "logit"){
    return(bigglm_updateQR_logit().update(qr, X, eta, offset, y));
  } else if (model == "exponential"){
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