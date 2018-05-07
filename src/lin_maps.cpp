#include "lin_maps.h"

using map_res_col = map_res<arma::subview_col<double>, arma::vec>;
using map_res_mat = map_res<arma::subview<double>, arma::mat>;

/* dens_mapper */
map_res_col
  dens_mapper::map_(const arma::vec &x, const bool transpose, ptr_vec &ptr) const
  {
    if(transpose){
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
  (const arma::mat &X, side s, const bool transpose) const
  {
    ptr_mat ptr;

    if(transpose){
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




/* select_mapper */
map_res_col
  select_mapper::map_(const arma::vec &x, const bool transpose, ptr_vec &ptr) const
  {
    if(transpose){
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
  (const arma::mat &X, side s, const bool transpose) const
  {
    ptr_mat ptr;

    if(transpose){
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





/* inv_mapper */
map_res_col
  inv_mapper::map_(const arma::vec &x, const bool transpose, ptr_vec &ptr) const
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
  (const arma::mat &X, side s, const bool transpose) const
  {
    ptr_mat ptr;

    if(transpose){
      switch(s){
      case  left:
        ptr.reset(new arma::mat(A_LU.solve(X    , true)));
        break;
      case both: {
        // A^{-\top} X A^-1 = (A^{-\top} (A^{-\top} X)^\top)^\top
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
        // A^{-1}XA^{-\top} = (A^{-1}(A^{-1}X)^\top)^\top
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

/* inv_sub_mapper
 * TODO: This could like be done smarter... */
map_res_col
  inv_sub_mapper::map_(const arma::vec &x, const bool transpose, ptr_vec &ptr) const
  {
    ptr.reset(new arma::vec(R.map_inv(A_LU.solve(R.map(x), transpose))));

    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

const arma::mat& inv_sub_mapper::map() const
{
  return inv_mat;
}

map_res_mat inv_sub_mapper::map
  (const arma::mat &X, side s, const bool transpose) const
  {
    ptr_mat ptr;

    if(transpose){
      Rcpp::stop("transpose not implemented in `inv_sub_mapper::map`");

    } else {
      arma::mat tmp = X; // copy
      if(s == left or s == both){
        tmp = R.map_inv(A_LU.solve(R.map(tmp)));

      }
      if(s == right or s == both){
        // XR^\top A^{-\top}R = (R^\top A^{-1} RX^\top)^\top
        tmp = tmp.t();
        tmp = R.map_inv(A_LU.solve(R.map(tmp))).t();

      }

      ptr.reset(new arma::mat(std::move(tmp)));
    }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }



