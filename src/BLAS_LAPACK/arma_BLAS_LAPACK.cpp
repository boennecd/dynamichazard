#include "../arma_BLAS_LAPACK.h"
#include "../R_BLAS_LAPACK.h"
#include <sstream>

#define MIN(a,b) (((a)<(b))?(a):(b))

# define LAPACK_CHECK_ILLEGAL(info_sym, meth_name)                                           \
if(info_sym < 0){                                                                            \
  std::stringstream ss;                                                                      \
  ss << "The " << -info_sym << "-th argument to " << #meth_name  << " had an illegal value"; \
  Rcpp::stop(ss.str());                                                                      \
}

void chol_rank_one_update(arma::mat &R, arma::vec x){
  // WARNING:
  // R will be overwritten!
  arma::vec copy_x = x; // will be destructed after Fortran call

  int n = R.n_rows;

  R_BLAS_LAPACK::ddhazard_dchur(R.memptr(), copy_x.memptr(), n, n);
}

void square_tri_inv(arma::mat &out){
  int n = out.n_cols;
  R_BLAS_LAPACK::square_tri_inv(out.memptr(), n, n);
}

void symmetric_mat_chol(const arma::mat& A, arma::mat &out){
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
arma::vec solve_w_precomputed_chol
  (const arma::mat &chol_decomp, const arma::vec& B)
{
  arma::vec out = B; /* take copy */
  int n = out.n_elem;

  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, true  /* transpose */, n, 1);
  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, false /* don't transpose */, n, 1);

  return out;
}


arma::mat solve_w_precomputed_chol
  (const arma::mat &chol_decomp, const arma::mat& B)
{
  arma::mat out = B; /* take copy */
  int n = out.n_rows, p = B.n_cols;

  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, true  /* transpose */, n, p);
  R_BLAS_LAPACK::triangular_sys_solve(
    chol_decomp.memptr(), out.memptr(), true, false /* don't transpose */, n, p);

  return out;
}


LU_factorization::LU_factorization(const arma::mat& A):
  M(A.n_rows), N(A.n_cols), has_elem(M > 0 && N > 0),
  A(new double[M * N]), IPIV(new int[MIN(M, N)]){
  if(!has_elem)
    return;

  // copy A
  const double *a = A.memptr();
  for(int i = 0; i < M * N; ++i, ++a)
    this->A[i] = *a;

  int LDA = M, INFO;
  R_BLAS_LAPACK::dgetrf(&M, &N, &this->A[0], &LDA, &IPIV[0], &INFO);

  LAPACK_CHECK_ILLEGAL(INFO, dgetrf)
  if(INFO > 0){
    std::stringstream ss;
    ss << "U(" << INFO << ", " << INFO << ") is exactly zero in dgetrf";
    Rcpp::stop(ss.str());
  }
}

arma::mat LU_factorization::solve() const {
  if(!has_elem)
    return arma::mat();

  if(M != N)
    Rcpp::stop("Non-square matrix in `LU_factorization::solve()`");

  // take copy
  arma::mat out(&A[0], N, N);

  int LWORK = N * N, INFO, LDA = N;
  std::unique_ptr<double []>  dwo(new double[LWORK]);

  // see https://stackoverflow.com/a/3520106/5861244 for example
  R_BLAS_LAPACK::dgetri(
    &N, &out[0], &LDA, &IPIV[0], &dwo[0], &LWORK, &INFO);
  LAPACK_CHECK_ILLEGAL(INFO, dgetri)

  return out;
}

arma::mat LU_factorization::solve(
    const arma::mat &B, const bool transpose) const {
  if(!has_elem && B.n_elem == 0)
    return arma::mat();

  // take copy
  arma::mat out = B;

  int NRHS = B.n_cols, LDA = N, LDB = B.n_rows, INFO;

  R_BLAS_LAPACK::dgetrs(
    transpose ? "T" : "N", &N, &NRHS, &A[0], &LDA, &IPIV[0], out.memptr(),
    &LDB, &INFO);
  LAPACK_CHECK_ILLEGAL(INFO, dgetrs)

  return out;
}

arma::vec LU_factorization::solve(
    const arma::vec& B, const bool transpose) const {
  if(!has_elem && B.n_elem == 0)
    return arma::vec();

  // take copy
  arma::vec out = B;

  int NRHS = 1, LDA = N, LDB = B.n_elem, INFO;

  R_BLAS_LAPACK::dgetrs(
    transpose ? "T" : "N", &N, &NRHS, &A[0], &LDA, &IPIV[0], out.memptr(),
    &LDB, &INFO);
  LAPACK_CHECK_ILLEGAL(INFO, dgetrs)

  return out;
}



QR_factorization::QR_factorization(const arma::mat &A):
  M(A.n_rows), N(A.n_cols), qr(new double[M * N]),
  qraux(new double[MIN(M, N)]), pivot_(new int[N]){
  // copy A
  const double *a = A.memptr();
  for(int i = 0; i < M * N; ++i, ++a)
    qr[i] = *a;

  /* initalize */
  for(int i = 0; i < N; ++i)
    pivot_[i] = 0;

  /* compute QR */
  int info, lwork = -1;
  double tmp;
  R_BLAS_LAPACK::dgeqp3(
    &M, &N, &qr[0], &M, &pivot_[0], &qraux[0], &tmp, &lwork, &info);
  LAPACK_CHECK_ILLEGAL(info, dgeqp3)

  lwork = (int) tmp;
  std::unique_ptr<double []> dwo(new double[lwork]);
  R_BLAS_LAPACK::dgeqp3(
    &M, &N, &qr[0], &M, &pivot_[0], &qraux[0], &dwo[0], &lwork, &info);
  LAPACK_CHECK_ILLEGAL(info, dgeqp3)

  rank = MIN(M, N);
}

arma::mat QR_factorization::qy(
    const arma::mat &B, const bool transpose) const {
  // take copy
  arma::mat out = B;
  int NRHS = B.n_cols, K = MIN(M, N);
  if(B.n_rows != (unsigned int)M)
    Rcpp::stop("Invalid `B` matrix in `QR_factorization::qy`");

  /* compute QR */
  int info, lwork = -1;
  double tmp;
  R_BLAS_LAPACK::dormqr(
    "L", transpose ? "T" : "N", &M, &NRHS, &K, &qr[0], &M, &qraux[0],
    out.memptr(), &M, &tmp, &lwork, &info);
  LAPACK_CHECK_ILLEGAL(info, dormqr)

  lwork = (int) tmp;
  std::unique_ptr<double []> work(new double[lwork]);
  R_BLAS_LAPACK::dormqr(
    "L", transpose ? "T" : "N", &M, &NRHS, &K, &qr[0], &M, &qraux[0],
    out.memptr(), &M, &work[0], &lwork, &info);
  LAPACK_CHECK_ILLEGAL(info, dormqr)

  return out;
}

arma::vec QR_factorization::qy(
    const arma::vec &B, const bool transpose) const {
  arma::mat out = B;
  return qy(out, transpose);
}

arma::mat QR_factorization::R() const {
  arma::mat out(&qr[0], M, N);
  out = out.rows(0, MIN(M, N) - 1);

  return arma::trimatu(out);
}

arma::uvec QR_factorization::pivot() const {
  arma::uvec out(N);
  for(int i = 0; i < N; ++i)
    out[i] = pivot_[i] - 1; /* want zero index */

  return out;
}
