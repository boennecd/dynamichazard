#include "dynamichazard.h"


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




template<class T>
class bigglm_updateQR{
  // match logic update.bigqr
  static arma::vec linkinv(const arma::vec &eta){
    return(eta.for_each([](double &val ){ return(T::link_func_inv(val)); }));
  }

  static arma::vec d_mu_d_eta(const arma::vec &eta){
    return(eta.for_each([](double &val ){ return(T::d_mu_d_eta(val)); }));
  }

  static arma::vec variance(const arma::vec &mu){
    return(mu.for_each([](double &val ){ return(T::variance(val)); }));
  }

public:
  void update(qr_obj &qr, // Previous/starting value. Will be overwritten
              const arma::mat &X, const arma::vec &eta,
              const arma::vec &offset, const arma::vec &y){
    //TODO: look into fortran memory storage versus C++ / armadillo

    arma::vec eta_plus_off = eta + offset;
    arma::vec mu = linkinv(eta_plus_off);
    arma::vec dmu = d_mu_d_eta(eta_plus_off);
    arma::vec z = eta + (y - mu) / dmu; // note that offset is not added as in bigglm.function
    arma::vec ww = dmu % dmu % variance(mu);

    int n_parems = X.n_rows;
    int nrbar = qr.rbar->n_elem;
    int ier = 0;



    for(int i = 0; i < ww.n_cols; ++i){
      // created as subroutine will overwrite values
      arma::vec x_row(X.colptr(i), X.n_rows);

      F77_CALL(includ)(&n_parems, &nrbar, z(i), x_row.memptr(),
               y.memptr(), qr.D->memptr(), qr.rbar->memptr(), qr.thetab->memptr(),
               &qr.ss, &ier)
    }
  }
};

using bigglm_updateQR_logit = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;


arma::vec bigglm_regcf(qr_obj &qr){
  int p = qr.D->n_elem;
  arma::vec beta(p, arma::fill::ones);
  F77_CALL(regcf)(p, p * p / 2, qr.D->memptr(), qr.rbar->memptr(),
           qr.thetab->memptr(), qr.tol->memptr(), beta.memptr(),
           p, 1L);

  return beta;
}


// Only exported for tests
// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset, const arma::vec &y){
  qr_obj qr;
  qr.D = &D;
  qr.rbar = &rbar;
  qr.thetab = &thetab;
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = &tol;

  if(model == "logit"){
    return(bigglm_updateQR_logit::update(qr, X, eta, offset, y));
  } else if (model == "exponential"){
    return(bigglm_updateQR_poisson::update(qr, X, eta, offset, y));
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
  qr.ss = &ss;
  qr.checked = &checked;
  qr.tol = &tol;

  return(bigglm_regcf(qr));
}
