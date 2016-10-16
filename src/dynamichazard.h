// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <thread>
#include <future>
#include <algorithm>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#if defined(USE_OPEN_BLAS) // Used to set the number of threads later
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

// DEBUG flags
// Only define this if the observational equation and state space equation
// have reasonable dimensions
//#define MYDEBUG_UKF

// #define MYDEBUG_EKF
// #define MYDEBUG_M_STEP
