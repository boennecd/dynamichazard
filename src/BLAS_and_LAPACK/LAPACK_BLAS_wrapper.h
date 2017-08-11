void ddhazard_dchur(double*, double*, int, int);

// wrapper for dtrtri
void square_tri_inv(double*, int, int);

// wrapper for dpotrf
void symmetric_mat_chol(double*, int, int);

// wrapper for dtrmv
void tri_mat_times_vec(double*, double*, int , int, bool);

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
