#include <Rcpp.h>
#include "../Rconfig-wrap.h"
#include <R_ext/Lapack.h>
#include "../R_BLAS_LAPACK.h"

static const char C_N = 'N', C_L = 'L';

namespace R_BLAS_LAPACK {
  extern "C"
  {
    void F77_NAME(dchur)(
        const char*,   // UPLO
        const char*,   // TRANS
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
        int*     // INFO,
        FCLEN FCLEN
    );

    void F77_NAME(dsyr)(
      const char *uplo, const int *n, const double *alpha,
      const double *x, const int *incx,
      double *a, const int *lda FCLEN);
  }

  void ddhazard_dchur(double *R, double *x, int n, int ldr){

    int info;

    double *c = new double[n];
    double *s = new double[n];

    // unused dummies
    int m = 0;
    int ldz = 1;
    double z, y, rho;

    F77_CALL(dchur)(
        &C_L,   // lower triangular
        &C_N,   // does not matter. Relates to Z
        &n, &m, R, &ldr,
        x, &z, &ldz, &y, &rho, c, s, &info FCONE FCONE);

    delete[] c;
    delete[] s;

    if(info != 0){
      std::stringstream str;
      str << "Got error code '" << info << "' when making rank one update of cholesky decomposition";
      Rcpp::stop(str.str());
    }
  }

  void square_tri_inv(double *out, int n, int ldr){
    int info;

    char c1 = 'L', c2 = 'N';
    F77_CALL(dtrtri)(&c1, &c2, &n, out, &ldr, &info FCONE FCONE);

    if(info != 0){
      std::stringstream str;
      str << "Got error code '" << info << "' from 'dtrtri'";
      Rcpp::stop(str.str());
    }
  }

  void symmetric_mat_chol(double *out, int n, int lda){
    int info;

    char c1 = 'L';
    F77_CALL(dpotrf)(&c1, &n, out, &lda, &info FCONE);

    if(info != 0){
      std::stringstream str;
      str << "Got error code '" << info << "' when making cholesky decomposition of symmetric matrix";
      Rcpp::stop(str.str());
    }
  }

  void tri_mat_times_vec(const double *A, double *x, int n, int lda, bool is_transpose){
    // Computes x <-- A * x where A is a triangular matrix
    // Computes x <-- A^T * x if TRANS = 'T' or 't' and not it 'N' or 'n'
    // A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
    // X is DOUBLE PRECISION array of dimension at least
    // ( 1 + ( n - 1 )*abs( INCX ) ).

    int incx = 1;
    char c1 = 'L', c2 = is_transpose ? 'T' : 'N', c3 = 'N';
    F77_CALL(dtrmv)(
        &c1, &c2, &c3, &n, A, &lda, x, &incx FCONE FCONE FCONE);
  }

  void sym_mat_rank_one_update(
      const int *n, const double *alpha, const double *x, double *A){
    // computes A := alpha * x * x^T + A
    // where A is a n x n is a square matrix
    //       x is a n matrix
    // Note that only the upper part will be updated!
    // See http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga35ca25bb135cd7bfdd5d6190b1aa4d07.html#ga35ca25bb135cd7bfdd5d6190b1aa4d07

    int inx = 1;
    char c1 = 'U';

    //F77_CALL(dger)(const int *m, const int *n, const double *alpha,
    // const double *x, const int *incx,
    // const double *y, const int *incy,
    // double *a, const int *lda);
    F77_CALL(dsyr)(
        &c1, n, alpha,
        x, &inx,
        A, n FCONE);
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
    char c1 = 'U';

    F77_CALL(dsymv)(
        &c1, n,
        &dum_d, A, n,
        x, &dum_i, &dum_d,
        y, &dum_i FCONE);
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
    char c1 = 'U', c2 = 'N';

    F77_CALL(dsyrk)(
        &c1, &c2,
        n, k, &dum_d,
        A, n,
        &dum_d, C, n FCONE FCONE);
  }

  void triangular_sys_solve(
      const double *A, double *B, const bool is_upper, const bool trans, const int n,
      const int nrhs){
    /*
       DTRTRS solves a triangular system of the form

       A * X = B  or  A**T * X = B,

       where A is a triangular matrix of order N, and B is an N-by-NRHS
       matrix.  A check is made to verify that A is nonsingular.
     */

    int info;
    char c1 = is_upper ? 'U' : 'L', c2 = trans ? 'T' : 'N', c3 = 'N';

    F77_CALL(dtrtrs)(
        &c1, &c2, &c3,
        &n, &nrhs,
        A, &n,
        B, &n,
        &info FCONE FCONE FCONE);

    if(info != 0){
      std::stringstream str;
      str << "Got error code '" << info << "' when using LAPACK dtrtrs";
      Rcpp::stop(str.str());
    }
  }

  void dtrmm(
      const char *side, const char *uplo, const char *transa,
      const char *diag, const int *m, const int *n,
      const double *alpha, const double *a, const int *lda,
      double *b, const int *ldb){
    F77_CALL(dtrmm)(
        side, uplo, transa, diag, m, n, alpha, a, lda,
        b, ldb FCONE FCONE FCONE FCONE);
  }

  void dtrmv(
      const char *uplo, const char *trans, const char *diag,
      const int *n, const double *a, const int *lda,
      double *x, const int *incx){
    F77_CALL(dtrmv)(
        uplo, trans, diag, n, a, lda, x, incx FCONE FCONE FCONE);
  }

  void dger(
      const int *m, const int *n, const double *alpha,
      const double *x, const int *incx,
      const double *y, const int *incy,
      double *a, const int *lda){
    F77_CALL(dger)(
        m, n, alpha,
        x, incx,
        y, incy,
        a, lda);
  }

  double ddot(
      const int *n, const double *dx, const int *incx,
      const double *dy, const int *incy){
    return F77_CALL(ddot)(n, dx, incx, dy, incy);
  }

  void dgetrf(
      const int* m, const int* n, double* a, const int* lda,int* ipiv,
      int* info){
    F77_CALL(dgetrf)(m, n, a, lda, ipiv, info);
  }

  void dgetrs(
      const char* trans, const int* n, const int* nrhs, const double* a,
      const int* lda, const int* ipiv, double* b, const int* ldb, int* info){
    F77_CALL(dgetrs)(trans, n, nrhs, a, lda, ipiv, b, ldb, info FCONE);
  }

  void dormqr(const char* side, const char* trans,
              const int* m, const int* n, const int* k,
              const double* a, const int* lda,
              const double* tau, double* c, const int* ldc,
              double* work, const int* lwork, int* info){
    F77_CALL(dormqr)(
        side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info
        FCONE FCONE);
  }

  void dgeqp3(const int* m, const int* n, double* a, const int* lda,
              int* jpvt, double* tau, double* work, const int* lwork,
              int* info){
    F77_CALL(dgeqp3)(m, n, a, lda, jpvt, tau, work, lwork, info);
  }

  void dgetri(const int* n, double* a, const int* lda,
              int* ipiv, double* work, const int* lwork, int* info){
    F77_CALL(dgetri)(n, a, lda, ipiv, work, lwork, info);
  }

  void dsyr(const char *uplo, const int *n, const double *alpha,
            const double *x, const int *incx,
            double *a, const int *lda){
    F77_CALL(dsyr)(uplo, n, alpha, x, incx, a, lda FCONE);
  }
}
