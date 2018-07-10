#ifndef LIN_MAPS
#define LIN_MAPS

#include "arma_n_rcpp.h"
#include "arma_BLAS_LAPACK.h"

/* Object used for map function. The unique pointer need to be to a new object
 * created in the functions to make sure the orginal object is not destructed.
 */

template <typename T_view, typename T_type>
class map_res {
  using ptr_T = std::unique_ptr<T_type>;

public:
  T_view sv; /* sv for (s)ub(v)iew */
  ptr_T org_ptr;

  map_res(T_view sv): sv(sv) {}
  map_res(T_view sv, ptr_T &ptr): sv(sv) {
    org_ptr = std::move(ptr);
  }
};

using map_res_col = map_res<arma::subview_col<double>, arma::vec>;
using map_res_mat = map_res<arma::subview<double>, arma::mat>;

enum side { left, both, right };
enum do_trans { dont_trans = 0, trans = 1 };

/* Compute products with matrix B */
class linear_mapper {
protected:
  using ptr_vec = std::unique_ptr<arma::vec>;
  using ptr_mat = std::unique_ptr<arma::mat>;

  virtual
    map_res_col map_(
        const arma::vec&, do_trans, ptr_vec&) const = 0;

public:
  virtual const arma::mat& map() const = 0;

  // create a virtual, default destructor
  virtual ~linear_mapper() = default;

  map_res_col map(arma::vec &x, do_trans transpose = dont_trans) const {
    ptr_vec ptr;
    return map_(x, transpose, ptr);
  }
  map_res_col map(arma::subview_col<double> x,
                  do_trans transpose = dont_trans) const {
    ptr_vec ptr;
    return map_(x, transpose, ptr);
  }
  map_res_col map(const arma::vec &x, do_trans transpose = dont_trans) const {
    /* TODO: re-implement to avoid copy */
    ptr_vec ptr(new arma::vec(x));
    return map_(x, transpose, ptr);
  }

  virtual map_res_mat map(const arma::mat&, side s = both,
                          do_trans transpose = dont_trans) const = 0;

  /* returns indicies for which rows and columns have non-zero entries */
  virtual const arma::uvec& non_zero_row_idx() const = 0;
  virtual const arma::uvec& non_zero_col_idx() const = 0;
};




#define PROTECTED_OVERIDES                                           \
map_res_col map_(const arma::vec&, do_trans, ptr_vec&) const override;

#define PUBLIC_OVERIDES                                           \
const arma::mat& map() const override;                            \
using linear_mapper::map;                                         \
map_res_mat map(const arma::mat&, side, do_trans) const override; \
const arma::uvec& non_zero_row_idx() const override;              \
const arma::uvec& non_zero_col_idx() const override;


/* Dens matrix B = A */
class dens_mapper: public linear_mapper {
  const arma::mat A;

  PROTECTED_OVERIDES

public:
  dens_mapper(const arma::mat &A) : A(A) {}

  PUBLIC_OVERIDES
};




/* Matrix B = R has a subset of identity matrix columns */
class select_mapper : public linear_mapper {
  const selection_matrix A;

  PROTECTED_OVERIDES

public:
  select_mapper(const selection_matrix &A) : A(A) {}

  PUBLIC_OVERIDES
};




/* Dense matrix B = A^{-1} */
class inv_mapper : public linear_mapper {
  const LU_factorization A_LU;
  const arma::mat A_inv;

  PROTECTED_OVERIDES

public:
  inv_mapper(const arma::mat &A) :
  A_LU(LU_factorization(A)), A_inv(A_LU.solve()) {}

  PUBLIC_OVERIDES
};




/* Dense matrix B = R^\top A^{-1} R where A is dense and R has a subset
 * of identity matrix columns */
class inv_sub_mapper : public linear_mapper {
  const LU_factorization A_LU;
  const selection_matrix R;
  const arma::mat inv_mat;

  PROTECTED_OVERIDES

public:
  inv_sub_mapper(const arma::mat &A, const selection_matrix &R):
  A_LU(LU_factorization(A)), R(R), inv_mat(R.map_inv(A_LU.solve())) {}

  PUBLIC_OVERIDES
};

#undef PROTECTED_OVERIDES
#undef PUBLIC_OVERIDES

#endif // LIN_MAPSs
