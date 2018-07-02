#include "bigglm_wrapper.h"
#include "R_ext/RS.h"		/* for F77_... */

extern "C"
{
  /*
  *
  * The Fortran code is:
  *      SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D,
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
  void F77_NAME(includ)(int *, int *, double *, double *, double *,
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
  void F77_NAME(regcf)(int *, int *, double *, double *, double *,
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



arma::vec bigglm_updateQR::linkinv(
    const arma::vec &eta, const arma::vec &exp_eta,
    const arma::vec &at_risk_length, family_base &fam){
  arma::vec out(eta.n_elem);

  double *i_out = out.memptr();
  const double *i_eta = eta.memptr();
  const double *i_exp_eta = exp_eta.memptr();

  if(fam.uses_at_risk_length()){
    const double *i_at = at_risk_length.memptr();
    for(arma::uword i = 0; i < out.n_elem;
        ++i, ++i_out, ++i_eta, ++i_exp_eta, ++i_at)
      *i_out = fam.linkinv(*i_eta, *i_exp_eta, *i_at);

  } else {
    for(arma::uword i = 0; i < out.n_elem; ++i, ++i_out, ++i_eta, ++i_exp_eta)
      *i_out = fam.linkinv(*i_eta, *i_exp_eta, 0);

  }

  return(out);
}

arma::vec bigglm_updateQR::d_mu_d_eta(
    const arma::vec &eta, const arma::vec &exp_eta,
    const arma::vec &at_risk_length, family_base &fam){
  arma::vec out(eta.n_elem);

  double *i_out = out.memptr();
  const double *i_eta = eta.memptr();
  const double *i_exp_eta = exp_eta.memptr();

  if(fam.uses_at_risk_length()){
    const double *i_at = at_risk_length.memptr();
    for(arma::uword i = 0; i < out.n_elem;
    ++i, ++i_out, ++i_eta, ++i_exp_eta, ++i_at)
      *i_out = fam.mu_eta(*i_eta, *i_exp_eta, *i_at);

  } else {
    for(arma::uword i = 0; i < out.n_elem; ++i, ++i_out, ++i_eta, ++i_exp_eta)
      *i_out = fam.mu_eta(*i_eta, *i_exp_eta, 0);

  }

  return(out);
}

arma::vec bigglm_updateQR::variance(
    const arma::vec &eta, const arma::vec &exp_eta,
    const arma::vec &at_risk_length, family_base &fam){
  arma::vec out(eta.n_elem);

  double *i_out = out.memptr();
  const double *i_eta = eta.memptr();
  const double *i_exp_eta = exp_eta.memptr();

  if(fam.uses_at_risk_length()){
    const double *i_at = at_risk_length.memptr();
    for(arma::uword i = 0; i < out.n_elem;
    ++i, ++i_out, ++i_eta, ++i_exp_eta, ++i_at)
      *i_out = fam.var(*i_eta, *i_exp_eta, *i_at);

  } else {
    for(arma::uword i = 0; i < out.n_elem; ++i, ++i_out, ++i_eta, ++i_exp_eta)
      *i_out = fam.var(*i_eta, *i_exp_eta, 0);

  }

  return(out);
}


void bigglm_updateQR::update(qr_obj &qr, // Previous/starting value. Will be overwritten
            const arma::mat &X, const arma::vec &eta,
            const arma::vec &offset, const arma::vec &at_risk_length,
            arma::vec &y, // y will not be altered
            const arma::vec &w, family_base &fam)
{
  arma::vec eta_plus_off = eta + offset;
  arma::vec exp_eta = arma::exp(eta_plus_off);
  arma::vec mu = linkinv(eta_plus_off, exp_eta, at_risk_length, fam);
  arma::vec dmu = d_mu_d_eta(eta_plus_off, exp_eta, at_risk_length, fam);
  arma::vec z = eta + (y - mu) / dmu; // note that offset is not added as in bigglm.function
  arma::vec ww = w % dmu % dmu /
    variance(eta_plus_off, exp_eta, at_risk_length, fam);

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

    F77_CALL(includ)(
        &n_parems, &nrbar, w_ptr, x_row.memptr(),
        y_ptr, qr.D->memptr(), qr.rbar->memptr(), qr.thetab->memptr(),
        &qr.ss, &ier);
  }
}

arma::vec bigglm_regcf(qr_obj &qr){
  int p = qr.D->n_elem;
  arma::vec beta(p, arma::fill::ones);

  int dum_arg1 = (p * p / 2);
  int dum_arg2 = 1;

  F77_CALL(regcf)(&p, &dum_arg1, qr.D->memptr(), qr.rbar->memptr(),
                  qr.thetab->memptr(), qr.tol->memptr(), beta.memptr(),
                  &p, &dum_arg2);

  return beta;
}

// Only exported for tests
// [[Rcpp::export]]
arma::vec bigglm_regcf_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                            double &ss, bool &checked, arma::vec &tol){
  qr_obj qr;
  qr.D = std::shared_ptr<arma::vec>(&D, [](arma::vec*x) -> void { });
  qr.rbar = std::shared_ptr<arma::vec>(&rbar, [](arma::vec*x) -> void { });
  qr.thetab = std::shared_ptr<arma::vec>(&thetab, [](arma::vec*x) -> void { });
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = std::shared_ptr<arma::vec>(&tol, [](arma::vec*x) -> void { });

  return(bigglm_regcf(qr));
}
