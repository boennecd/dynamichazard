#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

extern "C"
{
  // Non LAPACK function to make rank one update of of chol decomp
  // We could use the LINPACK function dchex. See
  //  http://www.netlib.org/linpack/dchex.f
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
}

void ddhazard_dchur(double *R, double *x, int n, int ldr){
  char uplo[] = "L";  // lower triangular
  char trans[] = "N"; // does not matter. Relates to Z

  int info;

  double *c = new double[n];
  double *s = new double[n];

  // unsused dummies
  int m = 0;
  int ldz = 1;
  double z, y, rho;

  dchur_(uplo, trans, &n, &m, R, &ldr,
         x, &z, &ldz, &y, &rho, c, s, &info);

  delete[] c;
  delete[] s;

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
    Rcpp::stop(str.str());
  }
};

void square_tri_inv(double *out, int n, int ldr){
  char uplo[] = "L"; // lower triangular
  char diag[] = "N"; //  non-unit triangular

  int info;

  F77_CALL(dtrtri)(uplo, diag, &n, out, &ldr, &info);

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
    Rcpp::stop(str.str());
  }
};

void symmetric_mat_chol(double *out, int n, int lda){
  char uplo[] = "L"; // lower triangular
  int info;

  F77_CALL(dpotrf)(uplo, &n, out, &lda, &info);

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making cholesky decomposition of symmetric matrix";
    Rcpp::stop(str.str());
  }
};

void tri_mat_times_vec(double *A, double *x, int n, int lda, bool is_transpose){
  // Computes x <-- A * x where A is a triangular matrix
  // Computes x <-- A^T * x if TRANS = 'T' or 't' and not it 'N' or 'n'
  // A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  // X is DOUBLE PRECISION array of dimension at least
  // ( 1 + ( n - 1 )*abs( INCX ) ).

  char uplo[] = "L"; // Not upper triangular
  const char *trans;
  if(is_transpose){
    trans = "T";
  } else{
    trans = "N";
  }
  char diag[] = "N";  //  non-unit triangular

  int incx = 1;
  F77_CALL(dtrmv)(uplo, trans, diag, &n, A, &lda, x, &incx);
}

void sym_mat_rank_one_update(
    const int *n, const double *alpha, const double *x, double *A){
 // computes A := alpha * x * x^T + A
 // where A is a n x n is a square matrix
 //       x is a n matrix

  int inx = 1;

  //F77_NAME(dger)(const int *m, const int *n, const double *alpha,
  // const double *x, const int *incx,
  // const double *y, const int *incy,
  // double *a, const int *lda);
  F77_NAME(dger)(n, n, alpha,
                 x, &inx, x, &inx,
                 A, n);
}
