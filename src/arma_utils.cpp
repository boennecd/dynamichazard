#include "BLAS_and_LAPACK/arma_utils.h"
#include "BLAS_and_LAPACK/LAPACK_BLAS_wrapper.h"

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

arma::mat out_mat_prod(const arma::mat &A){
  int n = A.n_rows;
  int k = A.n_cols;

  arma::mat out(n, n, arma::fill::zeros);

  symmetric_rank_k_update(&n, &k, A.memptr(), out.memptr());

  return out;
}

/*
  solve_w_precomputed_chol computes
    A X = B
  where D = chol_decomp is a Cholesky decomposition such that D^T D = A

 Maybe use BLAS DTRSM as backsolve and forward solve in R (see
 r-source/src/main/array.c)
*/
template<>
arma::vec solve_w_precomputed_chol(const arma::mat &chol_decomp, const arma::vec& B){
  arma::vec out = B; /* take copy */
  int n = out.n_elem;

  triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, true  /* transpose */, n, 1);
  triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, false /* don't transpose */, n, 1);

  return out;
}
