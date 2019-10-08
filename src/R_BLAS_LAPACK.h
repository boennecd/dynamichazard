namespace R_BLAS_LAPACK {
  void ddhazard_dchur(double*, double*, int, int);

  // wrapper for dtrtri
  void square_tri_inv(double*, int, int);

  // wrapper for dpotrf
  void symmetric_mat_chol(double*, int, int);

  // wrapper for dtrmv
  void tri_mat_times_vec(const double*, double*, int , int, bool);

  // wrapper for dger
  void sym_mat_rank_one_update(const int*, const double*, const double*, double *);

  // wrapper for DSYMV
  void sym_mat_vec_mult(const int *n, const double *x, const double *A, double *y);

  // wrapper for DSYRK
  void symmetric_rank_k_update(
      const int *n, const int *k, const double *A, double *C);

  // wrapper for dtrtrs
  void triangular_sys_solve(
      const double*, double*, const bool, const bool, const int , const int);

  void dtrmm(
      const char *side, const char *uplo, const char *transa,
      const char *diag, const int *m, const int *n,
      const double *alpha, const double *a, const int *lda,
      double *b, const int *ldb);

  void dtrmv(
      const char *uplo, const char *trans, const char *diag,
      const int *n, const double *a, const int *lda,
      double *x, const int *incx);

  void dger(
      const int *m, const int *n, const double *alpha,
      const double *x, const int *incx,
      const double *y, const int *incy,
      double *a, const int *lda);

  double ddot(
      const int *n, const double *dx, const int *incx,
      const double *dy, const int *incy);

  void dgetrf(
      const int* m, const int* n, double* a, const int* lda,int* ipiv,
      int* info);

  void dgetrs(
      const char* trans, const int* n, const int* nrhs, const double* a,
      const int* lda, const int* ipiv, double* b, const int* ldb, int* info);

  void dormqr(const char* side, const char* trans,
              const int* m, const int* n, const int* k,
              const double* a, const int* lda,
              const double* tau, double* c, const int* ldc,
              double* work, const int* lwork, int* info);

  void dgeqp3(const int* m, const int* n, double* a, const int* lda,
              int* jpvt, double* tau, double* work, const int* lwork,
              int* info);

  void dgetri(const int*, double*, const int*,
              int*, double*, const int*, int*);

  void dsyr(const char*, const int*, const double*, const double*,
            const int*,double*, const int*);
}
