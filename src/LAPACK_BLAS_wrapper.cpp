#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

extern "C"
{
  // Non LAPACK function to make rank one update of of chol decomp
  // We could use the LINPACK function dchex. See
  //  http://www.netlib.org/linpack/dchex.f
  // See http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=2646
  // I use the macro from r-source/src/include/R_ext/RS.h
  void F77_NAME(dchur)(
      int *,   // UPPER
      int *,   // DOTRAN
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

  int info;

  double *c = new double[n];
  double *s = new double[n];

  // unused dummies
  int m = 0;
  int ldz = 1;
  double z, y, rho;

  int UPPER = false;
  int DOTRAN = false;

  F77_CALL(dchur)(
      &UPPER,   // lower triangular
      &DOTRAN,   // does not matter. Relates to Z
      &n, &m, R, &ldr,
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
  int info;

  F77_CALL(dtrtri)("L", "N", &n, out, &ldr, &info);

  if(info != 0){
    std::stringstream str;
    str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
    Rcpp::stop(str.str());
  }
};

void symmetric_mat_chol(double *out, int n, int lda){
  int info;

  F77_CALL(dpotrf)("L", &n, out, &lda, &info);

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

  int incx = 1;
  F77_CALL(dtrmv)("L", is_transpose ? "T" : "N", "N", &n, A, &lda, x, &incx);
}

void sym_mat_rank_one_update(
    const int *n, const double *alpha, const double *x, double *A){
  // computes A := alpha * x * x^T + A
  // where A is a n x n is a square matrix
  //       x is a n matrix
  // Note that only the upper part will be updated!
  // See http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga35ca25bb135cd7bfdd5d6190b1aa4d07.html#ga35ca25bb135cd7bfdd5d6190b1aa4d07

  int inx = 1;

  //F77_NAME(dger)(const int *m, const int *n, const double *alpha,
  // const double *x, const int *incx,
  // const double *y, const int *incy,
  // double *a, const int *lda);
  F77_NAME(dsyr)(
      "U", n, alpha,
      x, &inx,
      A, n);
}

void sym_mat_vec_mult(
    const int *n, const double *x,
    const double *A, double *y){
  /*
   DSYMV  performs the matrix-vector  operation

   y := alpha*A*x + beta*y,

   where alpha and beta are scalars, x and y are n element vectors and
   A is an n by n symmetric matrix.
   */

  // This is only for alpha = 1, beta = 1. We use the upper part of A
  // See http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga6ab49c8fa5e2608d9545483045bf3d03.html#ga6ab49c8fa5e2608d9545483045bf3d03

  const double dum_d = 1.0;
  const int dum_i = 1L;

  F77_NAME(dsymv)(
      "U", n,
      &dum_d, A, n,
      x, &dum_i, &dum_d,
      y, &dum_i);
}

void symmetric_rank_k_update(
    const int *n, const int *k, const double *A, double *C){
  /*
   DSYRK  performs one of the symmetric rank k operations

   C := alpha*A*A**T + beta*C,

   or

   C := alpha*A**T*A + beta*C,

   where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   and  A  is an  n by k  matrix in the first case and a  k by n  matrix
   in the second case.
   */

  // We use C := alpha*A*A**T + beta*C with alpha = beta = 1.
  // Notice that we set UPLO to U so only the upper part is updated

  const double dum_d = 1.0;

  F77_NAME(dsyrk)(
      "U", "N",
      n, k, &dum_d,
      A, n,
      &dum_d, C, n);
}
