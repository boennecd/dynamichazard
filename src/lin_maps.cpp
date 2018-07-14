#include "lin_maps.h"

using map_res_col = map_res<arma::subview_col<double>, arma::vec>;
using map_res_mat = map_res<arma::subview<double>, arma::mat>;

/* dens_mapper */
map_res_col
  dens_mapper::map_(const arma::vec &x, do_trans transpose, ptr_vec &ptr) const
  {
    if(transpose == trans){
      ptr.reset(new arma::vec(A.t() * x));
    } else
      ptr.reset(new arma::vec(A * x));

    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

const arma::mat& dens_mapper::map() const
{
  return A;
}

map_res_mat dens_mapper::map
  (const arma::mat &X, side s, do_trans transpose) const
  {
    ptr_mat ptr;

    if(transpose == trans){
      switch(s){
      case  left:
        ptr.reset(new arma::mat(A.t() * X));
        break;
      case both:
        ptr.reset(new arma::mat(A.t() * X * A));
        break;
      case right:
        ptr.reset(new arma::mat(        X * A));
        break;
      default:
        Rcpp::stop("'Side' not implemented");
      }
    } else
      switch(s){
      case  left:
        ptr.reset(new arma::mat(A * X));
        break;
      case both:
        ptr.reset(new arma::mat(A * X * A.t()));
        break;
      case right:
        ptr.reset(new arma::mat(    X * A.t()));
        break;
      default:
        Rcpp::stop("'Side' not implemented");
      }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

const arma::uvec& dens_mapper::non_zero_row_idx() const{
  Rcpp::stop("'dens_mapper::non_zero_row_idx' is not implemented");
}

const arma::uvec& dens_mapper::non_zero_col_idx() const{
  Rcpp::stop("'dens_mapper::non_zero_col_idx' is not implemented");
}




/* select_mapper */
map_res_col
  select_mapper::map_(const arma::vec &x, do_trans transpose, ptr_vec &ptr) const
  {
    if(transpose == trans){
      ptr.reset(new arma::vec(A.map_inv(x)));
    } else
      ptr.reset(new arma::vec(A.map    (x)));

    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

const arma::mat& select_mapper::map() const
{
  return A.A;
}

map_res_mat select_mapper::map
  (const arma::mat &X, side s, do_trans transpose) const
  {
    ptr_mat ptr;

    if(transpose == trans){
      switch(s){
      case  left:
        ptr.reset(new arma::mat(          A.map_inv(X)));
        break;
      case both:
        ptr.reset(new arma::mat(A.map_inv(A.map_inv(X), true)));
        break;
      case right:
        ptr.reset(new arma::mat(A.map_inv(          X , true)));
        break;
      default:
        Rcpp::stop("'Side' not implemented");
      }
    } else
      switch(s){
      case  left:
        ptr.reset(new arma::mat(      A.map(X)));
        break;
      case both:
        ptr.reset(new arma::mat(A.map(A.map(X), true)));
        break;
      case right:
        ptr.reset(new arma::mat(A.map(      X , true)));
        break;
      default:
        Rcpp::stop("'Side' not implemented");
      }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

const arma::uvec& select_mapper::non_zero_row_idx() const {
  return A.non_zero_row_idx();
}

const arma::uvec& select_mapper::non_zero_col_idx() const {
  return A.non_zero_col_idx();
}




/* inv_mapper */
map_res_col
  inv_mapper::map_(const arma::vec &x, do_trans transpose, ptr_vec &ptr) const
  {
    ptr.reset(new arma::vec(A_LU.solve(x, transpose)));

    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

const arma::mat& inv_mapper::map() const
{
  return A_inv;
}

map_res_mat inv_mapper::map
  (const arma::mat &X, side s, do_trans transpose) const
  {
    ptr_mat ptr;

    if(transpose == trans){
      switch(s){
      case  left:
        ptr.reset(new arma::mat(A_LU.solve(X    , true)));
        break;
      case both: {
        // A^{-\top} X A^{^-1} = (A^{-\top}X^\top A^{-1})^\top = (A^{-\top} (A^{-\top} X)^\top)^\top
        arma::mat tmp = A_LU.solve(X, true).t();
        ptr.reset(new arma::mat(A_LU.solve(tmp  , true).t()));
      } break;
      case right: {
        // XA^{-1} = (A^{-\top}X^\top)^\top
        arma::mat X_T = X.t();
        ptr.reset(new arma::mat(A_LU.solve(X_T, true).t()));
      } break;
      default:
        Rcpp::stop("'Side' not implemented");
      }
    } else
      switch(s){
      case  left:
        ptr.reset(new arma::mat(A_LU.solve(X)));
        break;
      case both: {
        //  A^{-1}XA^{-\top} = (A^{-1}X^\top A^{-\top})^\top = (A^{-1}(A^{-1}X)^\top)^\top
        arma::mat tmp = A_LU.solve(X).t();
        ptr.reset(new arma::mat(A_LU.solve(tmp).t()));
      } break;
      case right: {
        // XA^{-\top} = (A^{-1}X^\top)^\top
        arma::mat X_T = X.t();
        ptr.reset(new arma::mat(A_LU.solve(X_T).t()));
      } break;
      default:
        Rcpp::stop("'Side' not implemented");
      }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

const arma::uvec& inv_mapper::non_zero_row_idx() const {
  Rcpp::stop("'inv_mapper::non_zero_row_idx' is not implemented");
}

const arma::uvec& inv_mapper::non_zero_col_idx() const {
  Rcpp::stop("'inv_mapper::non_zero_col_idx' is not implemented");
}

/* inv_sub_mapper
 * TODO: This could like be done smarter... */
map_res_col
  inv_sub_mapper::map_(const arma::vec &x, do_trans transpose, ptr_vec &ptr) const
  {
    ptr.reset(new arma::vec(R.map_inv(A_LU.solve(x, transpose))));

    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

const arma::mat& inv_sub_mapper::map() const
{
  return inv_mat;
}

map_res_mat inv_sub_mapper::map
  (const arma::mat &X, side s, do_trans transpose) const
  {
    ptr_mat ptr;

    if(transpose == trans){
      Rcpp::stop("transpose not implemented in `inv_sub_mapper::map`");

    } else {
      arma::mat tmp = X; // copy
      if(s == left or s == both){
        /* \begin{align*}
         *    B  &= R^\top A^{-1}  \\
         *    BX &= R^\top A^{-1}X
         * \end{align*} */
        tmp = R.map_inv(A_LU.solve(tmp));

      }
      if(s == right or s == both){
        /* \begin{align*}
         *    B       &= R^\top A^{-1}              \\
         *    XB^\top &= (R^\top A^{-1}X^\top)^\top
         * \end{align*} */

        // XR^\top A^{-\top}R = (R^\top A^{-1} RX^\top)^\top
        arma::inplace_trans(tmp);
        tmp = R.map_inv(A_LU.solve(tmp)).t();

      }

      ptr.reset(new arma::mat(std::move(tmp)));
    }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

const arma::uvec& inv_sub_mapper::non_zero_row_idx() const {
  Rcpp::stop("'inv_sub_mapper::non_zero_row_idx' is not implemented");
}

const arma::uvec& inv_sub_mapper::non_zero_col_idx() const {
  Rcpp::stop("'inv_sub_mapper::non_zero_col_idx' is not implemented");
}

