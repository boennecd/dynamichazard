// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <thread>
#include <future>
#include <algorithm>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#if defined(USE_OPEN_BLAS)
// Used to set the number of threads later
#include "cblas.h"
extern void openblas_set_num_threads(int num_threads);
extern int openblas_get_num_threads();
#define ARMA_USE_OPENBLAS
#define ARMA_DONT_USE_WRAPPER
#else
#define ARMA_USE_BLAS
#endif

// we know these are avialble with all R installations
#define ARMA_USE_LAPACK

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF

// Rcpp has its own stream object which cooperates more nicely
// with R's i/o -- and as of Armadillo 2.4.3, we can use this
// stream object as well
#if !defined(ARMA_DEFAULT_OSTREAM)
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif

//#define ARMA_NO_DEBUG
// Note: This also disables the check in inv(A, B) for whether inversion is succesfull it seems http://arma.sourceforge.net/docs.html#inv
// from armadillo config.hpp
//// Uncomment the above line if you want to disable all run-time checks.
//// This will result in faster code, but you first need to make sure that your code runs correctly!
//// We strongly recommend to have the run-time checks enabled during development,
//// as this greatly aids in finding mistakes in your code, and hence speeds up development.
//// We recommend that run-time checks be disabled _only_ for the shipped version of your program.


#define ARMA_DONT_PRINT_ERRORS
// from armadillo config.hpp:
//// Comment out the above line if you don't want errors and warnings printed (eg. failed decompositions)

#include <RcppArmadillo.h> // has to come after defines: http://artax.karlin.mff.cuni.cz/r-help/library/RcppArmadillo/html/RcppArmadillo-package.html
// [[Rcpp::depends("RcppArmadillo")]]

// Classes for R like distributions families
class dist_family {
public:
  virtual double link_func(const double&) const = 0;
  virtual double link_func_inv(const double&) const = 0;
  virtual double variance(const double&) const = 0;
  virtual double d_mu_d_eta(const double&) const = 0; // d mu / d eta
  virtual double dev_resids(const double&, const double&, const double&) const = 0;

  virtual double time_offset(const double&) const = 0;
};

class logit_fam : public dist_family {
private:
  static constexpr double THRESH = 30.;
  static constexpr double MTHRESH = -30.;
  static constexpr double INVEPS = 1 / DOUBLE_EPS;

public:
  double link_func(const double&) const;
  double link_func_inv(const double&) const;
  double variance(const double&) const;
  double d_mu_d_eta(const double&) const; // d mu / d eta
  double dev_resids(const double&, const double&, const double&) const;

  double time_offset(const double&) const;
};

class poisson_fam : public dist_family
{
public:
  double link_func(const double&) const;
  double link_func_inv(const double&) const;
  double variance(const double&) const;
  double d_mu_d_eta(const double&) const; // d mu / d eta
  double dev_resids(const double&, const double&, const double&) const;

  double time_offset(const double&) const;
};


// Functions and classes for fixed effects. Similar to object in bigglm
int binomialCoeff(int n, int k);

class qr_obj{
public:
  qr_obj(unsigned int p):
  D(new arma::vec(p, arma::fill::zeros)), rbar(new arma::vec((p == 1)? 0 : binomialCoeff(p, 2), arma::fill::zeros)),
  thetab(new arma::vec(p, arma::fill::zeros)), ss(0.), checked(false),
  tol(new arma::vec(p, arma::fill::zeros))
  {}
  qr_obj() = default;

  std::shared_ptr<arma::vec> D;
  std::shared_ptr<arma::vec> rbar;
  std::shared_ptr<arma::vec> thetab;
  double ss;
  bool checked;
  std::shared_ptr<arma::vec> tol;
};





template<class T>
class bigglm_updateQR{
  // match logic in update.bigqr
  arma::vec linkinv(const arma::vec &eta);

  arma::vec d_mu_d_eta(const arma::vec &eta);

  arma::vec variance(const arma::vec &mu);

protected:
  T t;

public:
  bigglm_updateQR<T>(): t() {}

  void update(qr_obj &qr, // Previous/starting value. Will be overwritten
              const arma::mat &X, const arma::vec &eta,
              const arma::vec &offset, arma::vec &y, // y will not be altered
              const arma::vec &w);
};

using bigglm_updateQR_logit   = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;

arma::vec bigglm_regcf(qr_obj &qr);

bool is_exponential_model(std::string model);
