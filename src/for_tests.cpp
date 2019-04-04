/*
  Rcpp does not search for attributes in sub directories which is why some
  functions are in here. See:
    http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2015-March/008473.html
*/

#include "sample_funcs.h"
#include "PF/densities.h"

// [[Rcpp::export]]
arma::uvec sample_indices_test(int size, arma::vec probs){
  return(sample_indices(size, probs));
}

// [[Rcpp::export]]
arma::uvec systematic_resampling_test(const int size, arma::vec probs){
  return systematic_resampling(size, probs);
}

// [[Rcpp::export]]
arma::umat
  sample_n_count_replicas_indices_test(const int size, arma::vec probs)
  {
    std::map<arma::uword, arma::uword> tmp =
      sample_n_count_replicas<sample_indices>(size, probs);
    arma::umat out(tmp.size(), 2L);

    arma::uword j = 0L;
    for (auto it = tmp.begin(); it != tmp.end(); it++, ++j)
    {
      out(j, 0L) = it->first;
      out(j, 1L) = it->second;
    }

    return out;
  }

// [[Rcpp::export]]
arma::umat
  sample_n_count_replicas_systematic_test(const int size, arma::vec probs)
  {
    std::map<arma::uword, arma::uword> tmp =
      sample_n_count_replicas<systematic_resampling>(size, probs);
    arma::umat out(tmp.size(), 2L);

    arma::uword j = 0L;
    for (auto it = tmp.begin(); it != tmp.end(); it++, ++j)
    {
      out(j, 0L) = it->first;
      out(j, 1L) = it->second;
    }

    return out;
  }

// [[Rcpp::export]]
arma::vec mvrnorm_test(const arma::vec mu, const arma::mat sigma_chol){
  return mvrnorm(mu, sigma_chol);
}

// [[Rcpp::export]]
double dmvnrm_log_test(
    const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  return(dmvnrm_log(x, mean, sigma_chol_inv));
}

// [[Rcpp::export]]
arma::vec mvtrnorm_test(
    const arma::vec mu, const arma::mat sigma_chol, const int nu)
  {
    return mvtrnorm(mu, sigma_chol, nu);
  }

// [[Rcpp::export]]
double dmvtrm_log_test(
    const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv,
    const int nu)
  {
    return(dmvtrm_log(x, mean, sigma_chol_inv, nu));
  }

// -------------------------------------------------- //

#include "ddhazard.h"
#include "family.h"
#include "estimate_fixed_effects_M_step.h"
#include "bigglm_wrapper.h"

// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset,
                          const arma::vec &at_risk_length,
                          arma::vec &y,
                          const arma::vec &w){
  qr_obj qr;
  qr.D = std::shared_ptr<arma::vec>(&D, [](arma::vec*x) -> void { });
  qr.rbar = std::shared_ptr<arma::vec>(&rbar, [](arma::vec*x) -> void { });
  qr.thetab = std::shared_ptr<arma::vec>(&thetab, [](arma::vec*x) -> void { });
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = std::shared_ptr<arma::vec>(&tol, [](arma::vec*x) -> void { });

  if(model == "logit"){
    logistic fam;
    return(bigglm_updateQR::update(
        qr, X, eta, offset, at_risk_length, y, w, fam));
  } else if (is_exponential_model(model)){
    exponential fam;
    return(bigglm_updateQR::update(
        qr, X, eta, offset, at_risk_length, y, w, fam));
  }
}

// [[Rcpp::export]]
double SMA_hepler_logit_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y){
  logistic f;
  return SMA::compute_length(
    offset, coef1, coef2, w, y, 0., f);
}

// [[Rcpp::export]]
double SMA_hepler_exp_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y, const double length){
  exponential f;
  return SMA::compute_length(
    offset, coef1, coef2, w, y, length, f);
}

// -------------------------------------------------- //

#include "PF/PF_utils.h"

// [[Rcpp::export]]
Rcpp::List PF_cloud_to_rcpp_and_back(const Rcpp::List &rcpp_list){
  auto cpp_result = get_clouds_from_rcpp_list(rcpp_list);

  return get_rcpp_list_from_cloud(cpp_result);
}

// [[Rcpp::export]]
Rcpp::List test_get_ancestors(const Rcpp::List &rcpp_list){
  auto cpp_result = get_clouds_from_rcpp_list(rcpp_list);


  std::vector<std::set<arma::uword> > output =
    get_ancestors(cpp_result.forward_clouds);

  Rcpp::List out(output.size());
  unsigned int i = 0;
  for(auto tp : output)
    out[i++] = Rcpp::wrap(tp);

  return out;
}

double get_weight_from_dpair(const std::pair<double, double> &x){
  return std::get<0>(x);
}

double get_resample_weight_from_dpair(const std::pair<double, double> &x){
  return std::get<1>(x);
}

// [[Rcpp::export]]
Rcpp::List test_get_resample_idx_n_log_weight
  (const arma::vec &log_weights, const arma::vec &log_resample_weights,
   const arma::uvec &resample_idx)
{
  std::vector<std::pair<double, double> > data_holder;
  auto lrw = log_resample_weights.begin();
  for(auto x : log_weights)
    data_holder.emplace_back(x, *(lrw++));

  std::map<arma::uword, double> out =
    get_resample_idx_n_log_weight
    <std::vector<std::pair<double, double> >, get_weight_from_dpair,
     get_resample_weight_from_dpair>
    (data_holder, resample_idx);

  unsigned int n_elem = out.size();
  arma::uvec idx(n_elem);
  arma::vec log_weights_out(n_elem);

  auto i = idx.begin();
  auto lw = log_weights_out.begin();
  for(auto x : out){
    *(i++)  = x.first;
    *(lw++) = x.second;

  }

  return Rcpp::List::create(
    Rcpp::Named("idx") = std::move(idx),
    Rcpp::Named("log_weights") = std::move(log_weights_out));
}

// -------------------------------------------------- //

#include "arma_BLAS_LAPACK.h"

// [[Rcpp::export]]
void chol_rank_one_update_test(arma::mat &R, arma::vec x){
  return chol_rank_one_update(R, x);
}

// [[Rcpp::export]]
arma::mat square_tri_inv_test(const arma::mat &R){
  arma::mat cp = R;
  square_tri_inv(cp);
  return cp;
}

// [[Rcpp::export]]
void symmetric_mat_chol_test(const arma::mat& A, arma::mat &out){
  return symmetric_mat_chol(A, out);
}

// [[Rcpp::export]]
void tri_mat_times_vec_test(arma::mat &A, const arma::vec &x, arma::vec &out, bool is_transpose){
  return tri_mat_times_vec(A, x, out, is_transpose);
}

// [[Rcpp::export]]
void sym_mat_rank_one_update_test(const double alpha, const arma::vec &x, arma::mat &A){
  return sym_mat_rank_one_update(alpha, x, A);
}

// [[Rcpp::export]]
Rcpp::List solve_w_precomputed_chol_test
  (const arma::mat &chol_decomp, const arma::vec& B, const arma::mat &D){
  return Rcpp::List::create(
    Rcpp::Named("B") = solve_w_precomputed_chol(chol_decomp, B),
    Rcpp::Named("D") = solve_w_precomputed_chol(chol_decomp, D));
}

// [[Rcpp::export]]
arma::mat solve_LU_inv(const arma::mat &A){
  LU_factorization fac(A);
  return fac.solve();
}

// [[Rcpp::export]]
arma::mat solve_LU_mat(const arma::mat &A, const arma::mat &B){
  LU_factorization fac(A);
  return fac.solve(B);
}

// [[Rcpp::export]]
arma::vec solve_LU_vec(const arma::mat &A, const arma::vec &B){
  LU_factorization fac(A);
  return fac.solve(B);
}

// [[Rcpp::export]]
arma::mat qr_qty_mat_test(const arma::mat &A, const arma::mat &B){
  QR_factorization fac(A);
  return fac.qy(B, true);
}

// [[Rcpp::export]]
arma::vec qr_qty_vec_test(const arma::mat &A, const arma::vec &B){
  QR_factorization fac(A);
  return fac.qy(B, true);
}

// [[Rcpp::export]]
arma::mat qr_R_test(const arma::mat &A){
  return QR_factorization(A).R();
}

// [[Rcpp::export]]
arma::mat selection_matrix_map_mat_test(arma::mat L, arma::mat X, bool is_right, bool is_inv){
  selection_matrix S_L(L);
  if(is_inv)
    return S_L.map_inv(X, is_right);

  return S_L.map(X, is_right);
}

// [[Rcpp::export]]
arma::vec selection_matrix_map_vec_test(arma::mat L, arma::vec X, bool is_inv){
  selection_matrix S_L(L);
  if(is_inv)
    return S_L.map_inv(X);

  return S_L.map(X);
}

// -------------------------------------------------- //

#include "utils.h"

// [[Rcpp::export]]
double lambert_W0_test(const double x){
  return lambert_W0(x);
}

// [[Rcpp::export]]
Rcpp::List trunc_eta_exponential_test(
  const double eta, const double at_risk_length, const bool is_event)
{
  auto ans = trunc_eta_exponential(is_event, eta, exp(eta), at_risk_length);

  return Rcpp::List::create(
    Rcpp::Named("eta_trunc") = ans.eta_trunc,
    Rcpp::Named("exp_eta_trunc") = ans.exp_eta_trunc);
}

// [[Rcpp::export]]
double trunc_eta_exponential_test_log_eps(){
  return trunc_eta_exponential_log_eps;
}

// -------------------------------------------------- //

#include "lin_maps.h"

// [[Rcpp::export]]
Rcpp::List linear_mapper_test(
    const arma::mat &A, const arma::vec x, const arma::mat X,
    const arma::vec z, const arma::mat Z, std::string type,
    const arma::mat R){
  bool has_R = R.n_elem > 0;
  std::unique_ptr<linear_mapper> ptr;

         if (type == "dens_mapper"){
    ptr.reset(new dens_mapper(A));
  } else if (type == "select_mapper"){
    ptr.reset(new select_mapper(A));
  } else if (type == "inv_mapper"){
    ptr.reset(new inv_mapper(A));
  } else if (type == "inv_sub_mapper"){
    ptr.reset(new inv_sub_mapper(A, R));
  } else
    Rcpp::stop("type not implemented");

  if(!has_R)
    return Rcpp::List::create(
      Rcpp::Named("A") = ptr->map(),

      Rcpp::Named("A_x")     = arma::vec(ptr->map(x       , dont_trans).sv),
      Rcpp::Named("A_T_z")   = arma::vec(ptr->map(z       , trans     ).sv),

      Rcpp::Named("A_X")     = arma::mat(ptr->map(X, left , dont_trans).sv),
      Rcpp::Named("X_A_T")   = arma::mat(ptr->map(X, right, dont_trans).sv),
      Rcpp::Named("A_X_A_T") = arma::mat(ptr->map(X, both , dont_trans).sv),

      Rcpp::Named("A_T_Z")   = arma::mat(ptr->map(Z, left , trans     ).sv),
      Rcpp::Named("Z_A")     = arma::mat(ptr->map(Z, right, trans     ).sv),
      Rcpp::Named("A_T_Z_A") = arma::mat(ptr->map(Z, both , trans     ).sv)
    );

  return Rcpp::List::create(
    Rcpp::Named("A") = ptr->map(),

    Rcpp::Named("A_x")     = arma::vec(ptr->map(x       , dont_trans).sv),

    Rcpp::Named("A_X")     = arma::mat(ptr->map(X, left , dont_trans).sv),
    Rcpp::Named("X_A_T")   = arma::mat(ptr->map(X, right, dont_trans).sv),
    Rcpp::Named("A_X_A_T") = arma::mat(ptr->map(X, both , dont_trans).sv)
  );
}

// -------------------------------------------------- //

#include "PF/dists.h"

// [[Rcpp::export]]
Rcpp::List check_state_fw(
  arma::vec parent, arma::vec parent1,
  arma::vec child, arma::vec child1,
  arma::mat F, arma::mat Q){
  covarmat cQ(Q);
  state_fw obj(parent, F, cQ);

  return Rcpp::List::create(
    Rcpp::Named("log_dens_func") =
      state_fw::log_dens_func(child, parent, F, cQ),

    Rcpp::Named("is_mvn") = obj.is_mvn(),
    Rcpp::Named("is_grad_z_hes_const") = obj.is_grad_z_hes_const(),
    Rcpp::Named("dim") = obj.dim(),

    Rcpp::Named("log_dens")  = obj.log_dens(child),
    Rcpp::Named("log_dens1") = obj.log_dens(child1),

    Rcpp::Named("gradient")  = obj.gradient(child),
    Rcpp::Named("gradient1") = obj.gradient(child1),

    Rcpp::Named("gradient_zero")  = obj.gradient_zero(&parent),
    Rcpp::Named("gradient_zero1") = obj.gradient_zero(&parent1),

    Rcpp::Named("neg_Hessian")  = obj.neg_Hessian(child),
    Rcpp::Named("neg_Hessian1") = obj.neg_Hessian(child1));
}

// [[Rcpp::export]]
Rcpp::List check_state_bw(
    arma::vec parent, arma::vec parent1,
    arma::vec child, arma::vec child1,
    arma::mat F, arma::mat Q){
  covarmat cQ(Q);
  state_bw obj(child, F, cQ);

  return Rcpp::List::create(
    Rcpp::Named("log_dens_func") =
      state_bw::log_dens_func(parent, child, F, cQ),

      Rcpp::Named("is_mvn") = obj.is_mvn(),
      Rcpp::Named("is_grad_z_hes_const") = obj.is_grad_z_hes_const(),
      Rcpp::Named("dim") = obj.dim(),

      Rcpp::Named("log_dens")  = obj.log_dens(parent),
      Rcpp::Named("log_dens1") = obj.log_dens(parent1),

      Rcpp::Named("gradient")  = obj.gradient(parent),
      Rcpp::Named("gradient1") = obj.gradient(parent1),

      Rcpp::Named("gradient_zero")  = obj.gradient_zero(&child),
      Rcpp::Named("gradient_zero1") = obj.gradient_zero(&child1),

      Rcpp::Named("neg_Hessian")  = obj.neg_Hessian(parent),
      Rcpp::Named("neg_Hessian1") = obj.neg_Hessian(parent1));
}

// [[Rcpp::export]]
Rcpp::List check_artificial_prior(
    arma::vec state,
    arma::mat F, arma::mat Q, arma::vec m_0, arma::mat Q_0,
    unsigned int t1, unsigned int t2, unsigned int t3){
  covarmat cQ(Q), cQ_0(Q_0);
  artificial_prior_generator gen(F, cQ, m_0, cQ_0);

  auto func = [&](unsigned int i){
    auto prior = gen.get_artificial_prior(i);

    return Rcpp::List::create(
        Rcpp::Named("is_mvn") = prior.is_mvn(),
        Rcpp::Named("is_grad_z_hes_const") = prior.is_grad_z_hes_const(),
        Rcpp::Named("dim") = prior.dim(),

        Rcpp::Named("log_dens")  = prior.log_dens(state),
        Rcpp::Named("gradient")  = prior.gradient(state),
        Rcpp::Named("gradient_zero")  = prior.gradient_zero(nullptr),
        Rcpp::Named("neg_Hessian")  = prior.neg_Hessian(state));
  };

  return Rcpp::List::create(
    Rcpp::Named(std::to_string(t1)) = func(t1),
    Rcpp::Named(std::to_string(t2)) = func(t2),
    Rcpp::Named(std::to_string(t3)) = func(t3));
}

// [[Rcpp::export]]
Rcpp::List check_observational_cdist(
    const arma::mat &X, const arma::uvec &is_event,
    const arma::vec &offsets, const arma::vec &tstart,
    const arma::vec &tstop, const double bin_start, const double bin_stop,
    const bool multithreaded, std::string fam,
    const arma::vec state, const arma::vec state1){
  std::shared_ptr<PF_cdist> dist =
    get_observational_cdist(
      fam, X, is_event, offsets, tstart, tstop, bin_start, bin_stop,
      multithreaded);

  return Rcpp::List::create(
    Rcpp::Named("is_mvn") = dist->is_mvn(),
    Rcpp::Named("dim") = dist->dim(),

    Rcpp::Named("log_dens")  = dist->log_dens(state),
    Rcpp::Named("log_dens1") = dist->log_dens(state1),

    Rcpp::Named("gradient")  = dist->gradient(state),
    Rcpp::Named("gradient1") = dist->gradient(state1),

    Rcpp::Named("neg_Hessian")  = dist->neg_Hessian(state),
    Rcpp::Named("neg_Hessian1") = dist->neg_Hessian(state1));
}

// -------------------------------------------------- //

#include "PF/cond_approx.h"

// [[Rcpp::export]]
Rcpp::List check_fw_bw_comb(
    arma::mat F, arma::mat Q,
    arma::vec parent, arma::vec parent1,
    arma::vec grand_child, arma::vec grand_child1,
    arma::vec x, int nu){
  covarmat cQ(Q);

  state_fw fw(parent     , F, cQ);
  state_bw bw(grand_child, F, cQ);

  std::vector<PF_cdist*> objs = { &fw, &bw };
  cdist_comb_generator combi_gen(objs, parent, nu);

  auto func = [&](arma::vec &p, arma::vec &gc){
    std::unique_ptr<dist_comb> comb = combi_gen.
      get_dist_comb({ &p, &gc });

    return Rcpp::List::create(
      Rcpp::Named("log_density") = comb->log_density(x),
      Rcpp::Named("mean") = comb->get_mean(),
      Rcpp::Named("covar") = comb->get_covar());
  };

  return Rcpp::List::create(
    Rcpp::Named("comb_1_1") = func(parent, grand_child),
    Rcpp::Named("comb_2_2") = func(parent1, grand_child1));
}

// [[Rcpp::export]]
Rcpp::List check_prior_bw_comb(
    arma::mat F, arma::mat Q, arma::vec m_0, arma::mat Q_0,
    arma::vec child, arma::vec child1, arma::vec parent,
    unsigned int t1, unsigned int t2){
  covarmat cQ(Q), cQ_0(Q_0);
  state_bw bw(child, F, cQ);
  artificial_prior_generator gen(F, cQ, m_0, cQ_0);

  auto func = [&](unsigned int i){
    auto prior = gen.get_artificial_prior(i);
    std::vector<PF_cdist*> objs = { &prior, &bw };
    cdist_comb_generator comb(objs);
    std::unique_ptr<dist_comb>
      d1 = comb.get_dist_comb({ &child  }),
      d2 = comb.get_dist_comb({ &child1 });

    return Rcpp::List::create(
      Rcpp::Named("mean1") = d1->get_mean(),
      Rcpp::Named("mean2") = d2->get_mean(),

      Rcpp::Named("covar1") = d1->get_covar(),
      Rcpp::Named("covar2") = d2->get_covar(),

      Rcpp::Named("log_dens1")  = d1->log_density(parent),
      Rcpp::Named("log_dens2")  = d2->log_density(parent));
  };

  return Rcpp::List::create(
    Rcpp::Named(std::to_string(t1)) = func(t1),
    Rcpp::Named(std::to_string(t2)) = func(t2));
}

// [[Rcpp::export]]
Rcpp::List check_prior_bw_state_comb(
    const arma::mat &X, const arma::uvec &is_event,
    const arma::vec &offsets, const arma::vec &tstart,
    const arma::vec &tstop, const double bin_start, const double bin_stop,
    std::string fam,
    arma::mat F, arma::mat Q, arma::vec m_0, arma::mat Q_0,
    arma::vec child, arma::vec child1, arma::vec parent,
    unsigned int t1, arma::mat Q_xtra, unsigned int nu = -1,
    const double covar_fac = -1, const double ftol_rel = 1e-8){
  const bool multithreaded = 1;
  covarmat cQ(Q), cQ_0(Q_0);

  state_bw bw(child, F, cQ);
  artificial_prior_generator gen(F, cQ, m_0, cQ_0);
  std::shared_ptr<PF_cdist> dist =
    get_observational_cdist(
      fam, X, is_event, offsets, tstart, tstop, bin_start, bin_stop,
      multithreaded);

  auto prior = gen.get_artificial_prior(t1);
  std::vector<PF_cdist*> objs_sta = { &prior, &bw };
  cdist_comb_generator comb_start(objs_sta);

  std::unique_ptr<dist_comb> d1, d2;
  if(Q_xtra.n_cols == Q.n_cols){
    std::vector<PF_cdist*> objs = { &prior, &bw, dist.get() };
    cdist_comb_generator comb(
        objs, comb_start.get_dist_comb({ &child })->get_mean(), nu, &Q_xtra,
        covar_fac, ftol_rel);
    d1 = comb.get_dist_comb({ &child   });
    d2 = comb.get_dist_comb({ &child1  });

  } else {
    std::vector<PF_cdist*> objs = { &prior, &bw, dist.get() };
    cdist_comb_generator comb(
        objs, comb_start.get_dist_comb({ &child })->get_mean(), nu);
    d1 = comb.get_dist_comb({ &child   });
    d2 = comb.get_dist_comb({ &child1  });

  }

  return Rcpp::List::create(
    Rcpp::Named("mean1") = d1->get_mean(),
    Rcpp::Named("mean2") = d2->get_mean(),

    Rcpp::Named("covar1") = d1->get_covar(),
    Rcpp::Named("covar2") = d2->get_covar(),

    Rcpp::Named("log_dens1")  = d1->log_density(parent),
    Rcpp::Named("log_dens2")  = d2->log_density(parent));
}

