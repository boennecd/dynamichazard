#include "../arma_BLAS_LAPACK.h"
#include "../R_BLAS_LAPACK.h"
#include <sstream>
#define MIN(a,b) (((a)<(b))?(a):(b))

void chol_rank_one_update(arma::mat &R, arma::vec x){
  // WARNING:
  // R will be overwritten!
  arma::vec copy_x = x; // will be destructed after Fortran call

  int n = R.n_rows;

  R_BLAS_LAPACK::ddhazard_dchur(R.memptr(), copy_x.memptr(), n, n);
}

void square_tri_inv(const arma::mat &R, arma::mat &out){
  // Take copy
  out = R;
  int n = out.n_cols;

  R_BLAS_LAPACK::square_tri_inv(out.memptr(), n, n);
}

void symmetric_mat_chol(const arma::mat& A, arma::mat & out){
  // Take copy
  out = A;

  int n = out.n_cols;

  R_BLAS_LAPACK::symmetric_mat_chol(out.memptr(), n, n);
  out = arma::trimatl(out); // keep only lower diagonal
}

void tri_mat_times_vec(const arma::mat &A, const arma::vec &x, arma::vec &out, bool is_transpose){
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

  R_BLAS_LAPACK::tri_mat_times_vec(A.memptr(), out.memptr(), n, n, is_transpose);
}

void tri_mat_times_vec(const arma::mat &A, arma::vec &out, bool is_transpose){
  int n = A.n_cols;
  R_BLAS_LAPACK::tri_mat_times_vec(A.memptr(), out.memptr(), n, n, is_transpose);
}

void sym_mat_rank_one_update(const double alpha, const arma::vec &x, arma::mat &A){
  int n = A.n_cols;

  R_BLAS_LAPACK::sym_mat_rank_one_update(&n, &alpha, x.memptr(), A.memptr());
}

arma::vec sym_mat_times_vec(const arma::vec &x, const arma::mat &A){
  int n = A.n_cols;
  arma::vec out(A.n_cols, arma::fill::zeros);

  R_BLAS_LAPACK::sym_mat_vec_mult(&n, x.memptr(), A.memptr(), out.memptr());

  return out;
}

arma::mat out_mat_prod(const arma::mat &A){
  int n = A.n_rows;
  int k = A.n_cols;

  arma::mat out(n, n, arma::fill::zeros);

  R_BLAS_LAPACK::symmetric_rank_k_update(&n, &k, A.memptr(), out.memptr());

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

  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, true  /* transpose */, n, 1);
  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, false /* don't transpose */, n, 1);

  return out;
}


LU_factorization::LU_factorization(const arma::mat& A):
  M(A.n_rows), N(A.n_cols), has_elem(M > 0 && N > 0),
  A_(new double[M * N]), IPIV_(new int[MIN(M, N)]){
  if(!has_elem)
    return;

  // copy A
  const double *a = A.memptr();
  for(int i = 0; i < M * N; ++i, ++a)
    A_[i] = *a;

  int LDA = M, INFO;
  R_BLAS_LAPACK::dgetrf(&M, &N, &A_[0], &LDA, &IPIV_[0], &INFO);

  if(INFO < 0){
    std::stringstream ss;
    ss << "The " << -INFO << "-th argument to dgetrf had an illegal value";
    Rcpp::stop(ss.str());
  }
  if(INFO > 0){
    std::stringstream ss;
    ss << "U(" << INFO << ", " << INFO << ") is exactly zero in dgetrf";
    Rcpp::stop(ss.str());
  }
}

arma::mat LU_factorization::solve(const arma::mat &B, bool tranpose) const {
  if(!has_elem && B.n_elem == 0)
    return arma::mat();

  // take copy
  arma::mat out = B;

  int NRHS = B.n_cols, LDA = N, LDB = B.n_rows, INFO;

  R_BLAS_LAPACK::dgetrs(
    tranpose ? "T" : "N", &N, &NRHS, &A_[0], &LDA, &IPIV_[0], out.memptr(),
    &LDB, &INFO);

  if(INFO < 0){
    std::stringstream ss;
    ss << "The " << -INFO << "-th argument to dgetrs had an illegal value";
    Rcpp::stop(ss.str());
  }

  return out;
}

arma::vec LU_factorization::solve(const arma::vec& B, bool tranpose) const {
  if(!has_elem && B.n_elem == 0)
    return arma::vec();

  // take copy
  arma::vec out = B;

  int NRHS = 1, LDA = N, LDB = B.n_elem, INFO;

  R_BLAS_LAPACK::dgetrs(
    tranpose ? "T" : "N", &N, &NRHS, &A_[0], &LDA, &IPIV_[0], out.memptr(),
    &LDB, &INFO);

  if(INFO < 0){
    std::stringstream ss;
    ss << "The " << -INFO << "-th argument to dgetrs had an illegal value";
    Rcpp::stop(ss.str());
  }

  return out;
}
