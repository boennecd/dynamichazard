#include "arma_utils.h"
#include "LAPACK_BLAS_wrapper.h"

// Exported for tests
// [[Rcpp::export]]
void chol_rank_one_update(arma::mat &R, arma::vec x){
  // WARNING:
  // R will be overwritten!
  arma::vec copy_x = x; // will be destructed after Fortran call

  int n = R.n_rows;

  ddhazard_dchur(R.memptr(), copy_x.memptr(), n, n);
}

// Exported for tests
// [[Rcpp::export]]
void square_tri_inv(const arma::mat &R, arma::mat &out){
  // Take copy
  out = R;
  int n = out.n_cols;

  square_tri_inv(out.memptr(), n, n);
};

// Exported for tests
// [[Rcpp::export]]
void symmetric_mat_chol(const arma::mat& A, arma::mat & out){
  // Take copy
  out = A;

  int n = out.n_cols;

  symmetric_mat_chol(out.memptr(), n, n);
  out = arma::trimatl(out); // keep only lower diagonal
}

// Exported for tests
// [[Rcpp::export]]
void tri_mat_times_vec(arma::mat &A, const arma::vec &x, arma::vec &out, bool is_transpose){
  // Computes x <-- A * x where A is a triangular matrix
  // Computes x <-- A^T * x if TRANS = 'T' or 't' and not it 'N' or 'n'
  // A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  // X is DOUBLE PRECISION array of dimension at least
  // ( 1 + ( n - 1 )*abs( INCX ) ).

  int n = A.n_cols;

  // Take copy
  const double *x_ptr = x.memptr();
  double *out_ptr = out.memptr();
  unsigned int i = 0;
  for(; i < x.n_elem; i++, out_ptr++, x_ptr++){
    *out_ptr = *x_ptr;
  }
  for(; i < out.n_elem; i++, out_ptr++)
    *out_ptr = 0.;

  tri_mat_times_vec(A.memptr(), out.memptr(), n, n, is_transpose);
}

void tri_mat_times_vec(arma::mat &A, arma::vec &out, bool is_transpose){
  int n = A.n_cols;
  tri_mat_times_vec(A.memptr(), out.memptr(), n, n, is_transpose);
}

// Exported for tests
// [[Rcpp::export]]
void sym_mat_rank_one_update(const double alpha, const arma::vec &x, arma::mat &A){
  int n = A.n_cols;

  sym_mat_rank_one_update(&n, &alpha, x.memptr(), A.memptr());
};

arma::vec sym_mat_times_vec(const arma::vec &x, const arma::mat &A){
  int n = A.n_cols;
  arma::vec out(A.n_cols, arma::fill::zeros);

  sym_mat_vec_mult(&n, x.memptr(), A.memptr(), out.memptr());

  return out;
}
