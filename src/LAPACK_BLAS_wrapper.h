void ddhazard_dchur(double*, double*, int, int);

// wrapper for dtrtri
void square_tri_inv(double*, int, int);

// wrapper for dpotrf
void symmetric_mat_chol(double*, int, int);

// wrapper for dtrmv
void tri_mat_times_vec(double*, double*, int , int, bool);

// wrapper for dger
void sym_mat_rank_one_update(const int*, const double*, const double*, double *);

// wrapper for
void sym_mat_vec_mult(const int *n, const double *x, const double *A, double *y);
