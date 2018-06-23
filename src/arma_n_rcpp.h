#ifndef DD_ARMA_N_RCPP
#define DD_ARMA_N_RCPP

#define _USE_MATH_DEFINES
#include <cmath>

// Don't use openMP for elementwise operations. It is automatically if -fopenmp
// is present. Seems to cause issues on some platforms
#ifndef ARMA_DONT_USE_OPENMP
  #define ARMA_DONT_USE_OPENMP 1
#endif

#define ARMA_NO_DEBUG
// Note: This also disables the check in inv(A, B) for whether inversion is succesfull it seems http://arma.sourceforge.net/docs.html#inv
// from armadillo config.hpp
//// Uncomment the above line if you want to disable all run-time checks.
//// This will result in faster code, but you first need to make sure that your code runs correctly!
//// We strongly recommend to have the run-time checks enabled during development,
//// as this greatly aids in finding mistakes in your code, and hence speeds up development.
//// We recommend that run-time checks be disabled _only_ for the shipped version of your program.

#include <RcppArmadillo.h>

// Print function to print out vectors as rows
template<typename T>
inline
void
my_print(const T &X, std::string msg = "", std::string prefix = "")
{
  // See armadillo-code/include/armadillo_bits/arma_ostream_meat.hpp
  std::stringstream os;

  const std::streamsize cell_width = 14;

  if(msg != "")
    os << prefix << msg << std::endl;

  arma::mat to_print = X;
  if(X.n_cols == 1)
    to_print = to_print.t();

  for(arma::uword row=0; row < to_print.n_rows; ++row)
  {
    os << prefix;

    for(arma::uword col=0; col < to_print.n_cols; ++col){
      os.width(cell_width);
      arma::arma_ostream::print_elem(os, to_print.at(row,col), true);
    }

    os << std::endl;
  }

  Rcpp::Rcout << os.str();
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
      if(!arma::inv_sympd(out, std::forward<T2>(X)))
        inv(out, std::forward<T2>(X), true, err_msg);
    } else if(!arma::inv_sympd(out, std::forward<T2>(X))){
      Rcpp::stop(err_msg);
    }
  }

#endif
