#ifndef DD_ARMA_N_RCPP
#define DD_ARMA_N_RCPP

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

// Print function to print out vectors as rows
template<typename T>
inline
  void
  my_print(const T& X, std::string msg = "")
  {
    if(msg != "")
      Rcpp::Rcout << msg << std::endl;

    if(X.n_cols > 1){
      X.print();

    } else{
      for(arma::uword col=0; col < X.n_cols; ++col)
      {
        for(arma::uword row=0; row < X.n_rows; ++row)
          Rcpp::Rcout << std::setw(10) << std::setprecision(5) <<  X(row,col) << ' ';

        Rcpp::Rcout << std::endl;
      }
    }
  }



// Inversion methods
template<typename eT, typename T2>
inline
  void inv(arma::Mat<eT>& out, T2&& X,
           const bool use_pinv = false,
           const std::string err_msg = ""){
    if(use_pinv){
      // Compute the pseudo inverse using SVD
      // Tolerance controls the value for which eigenvalues and vectors are
      // dropped
      // the default tolerance is max(m,n)*max_sv*datum::eps, where:
      //   m = number of rows and n = number of columns
      //   max_sv = maximal singular value
      //   datum::eps = difference between 1 and the least value greater than 1 that is representable

      // Two method are avialable for the SVD: "dc" and "std" (former is default)
      // See armadillo-#.###.#\include\armadillo_bits\op_pinv_meat.hpp
      // "dc" uses auxlib::svd_dc_econ which calls lapack::cx_gesdd
      // "std" uses auxlib::svd_econ which calls lapack::cx_gesvd
      // see armadillo-#.###.#\include\armadillo_bits\auxlib_meat.hpp
      // and armadillo-#.###.#\include\armadillo_bits\wrapper_lapack.hpp

      // For double, former calls zgesdd and latter calls arma_zgesvd

      // Acording to this post: https://groups.google.com/forum/#!topic/julia-dev/mmgO65i6-fA
      // "The routine dgesdd uses a divide and conquer algorithm for the
      //   bidiagonal SVD, whereas dgesvd uses a QR algorithm. The former is
      //   faster in cases where there is significant deflation and should be
      //   generally preferred for large matrices."

      // A next question is regarding whether the tolerance. Wikipedia suggest
      // this tolerance is the default in other langagues as well https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
      if(!arma::pinv(out, std::forward<T2>(X))){
        Rcpp::stop(err_msg);
      }
    } else{
      if(!arma::inv(out, std::forward<T2>(X))){
        Rcpp::stop(err_msg);
      }
    }
  }

template<typename eT, typename T2>
inline
  void inv_sympd(arma::Mat<eT>& out, T2&& X,
                 const bool use_pinv = false,
                 const std::string err_msg = ""){

    if(use_pinv){
      inv(out, std::forward<T2>(X), true, err_msg);
    } else if(!arma::inv_sympd(out, std::forward<T2>(X))){
      Rcpp::stop(err_msg);
    }
  }

#endif
