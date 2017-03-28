#include "LAPACK_wrapper.h"

extern "C"
{
  // Non LAPACK function to make rank one update of of chol decomp
  void dchur_(
      char*,   // UPLO
      char*,   // TRANS
      int*,    // N
      int*,    // M
      double*, // R
      int*,    // LDR
      double*, // X
      double*, // Z
      int*,    // LDZ
      double*, // Y
      double*, // RHO
      double*, // C
      double*, // S
      int*     // INFO
  );

  // computes the inverse of a real upper or lower triangular matrix A
  void dtrtri_(
      char*,   // UPLO
      char*,   // DIAG
      int*,    // N
      double*, // A
      int*,    // LDA
      int*     // INFO
  );

  // computes the Cholesky factorization of a real symmetric positive definite
  // matrix
  void dpotrf_(
    char*,    // UPLO
    int*,     // N
    double*,  // A
    int*,     // LDA
    int*      // INFO
  );

  // // Computes traingular matrix vector product
  // void strmv_(
  //     char*,    // UPLO
  //     char*,    // TRANS
  //     char*,    // DIAG
  //     int*,     // N
  //     double*,  // A
  //     int*,     // LDA
  //     double*,  // X
  //     int*      // INCX
  // );
}

// Exported for tests
// [[Rcpp::export]]
void chol_rank_one_update(arma::mat &R, arma::vec x){
  // WARNING:
  // R will be overwritten!
  arma::vec copy_x = x; // will be destructed after Fortran call

  char uplo[] = "L";  // lower triangular
  char trans[] = "N"; // does not matter. Relates to Z

  int n = R.n_rows;
  int ldr = n;
  int info;

  double *c = new double[n];
  double *s = new double[n];

  // unsused dummies
  int m = 0;
  int ldz = 1;
  double z, y, rho;

  dchur_(uplo, trans, &n, &m, R.memptr(), &ldr,
         copy_x.memptr(), &z, &ldz, &y, &rho, c, s, &info);

  delete[] c;
  delete[] s;

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
    Rcpp::stop(str.str());
  }
}

// Exported for tests
// [[Rcpp::export]]
void square_tri_inv(const arma::mat &R, arma::mat &out){
  // Take copy
  out = R;

  char uplo[] = "L"; // lower triangular
  char diag[] = "N"; //  non-unit triangular

  int n = out.n_cols;
  int lda = n;
  int info;

  dtrtri_(uplo, diag, &n, out.memptr(), &lda, &info);

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
    Rcpp::stop(str.str());
  }
};

// Exported for tests
// [[Rcpp::export]]
void symmetric_mat_chol(const arma::mat& A, arma::mat & out){
  // Take copy
  out = A;

  char uplo[] = "L"; // lower triangular
  int n = out.n_cols;
  int lda = n;
  int info;

  dpotrf_(uplo, &n, out.memptr(), &lda, &info);

  out = arma::trimatl(out); // keep only lower diagonal

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making cholesky decomposition of symmetric matrix";
    Rcpp::stop(str.str());
  }
}

// // Product betwen square
// inline void square_tri_vec_prod(arma::mat &A, const arma::vec &x, arma::vec &out){
//   char uplo[] = "L";  // lower triangular
//   char trans[] = "N"; // not transposed
//   char diag[] = "N"; //  non-unit triangular
//
//   int n = A.n_cols;
//   int lda = n;
//   int incx = 1;
//
//   // Take copy
//   out = x;
//
//   F77_CALL(strmv)(uplo, trans, diag, &n, A.memptr(), &lda, out.memptr(), &incx);
// }
//
// // Exported for tests
// // [[Rcpp::export]]
// void square_tri_vec_prod_test(arma::mat &A, const arma::vec &x, arma::vec &out){
//   square_tri_vec_prod(A, x, out);
//}

