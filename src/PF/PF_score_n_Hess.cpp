#include "PFs.h"
#include "../thread_pool.h"
#include "importance_samplers.h"
#include "resamplers.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

void print_means_score_n_hess
  (const std::vector<const score_n_hess_base*> &dat, const unsigned int t){
  arma::vec obs;
  arma::vec state;

  /* normalize */
  double norm = 0;
  for(auto i : dat)
    norm += i->get_weight();
  arma::vec ws(dat.size());
  auto w = ws.begin();
  for(auto i : dat)
    *(w++) = i->get_weight() / norm;

  /* compute weighted means */
  bool first = true;
  w = ws.begin();
  for(auto i : dat){
    if(first){
      obs   = *w * i->get_a_obs();
      state = *(w++) * i->get_a_state();
      first = false;
      continue;

    }

    obs   += *w * i->get_a_obs();
    state += *(w++) * i->get_a_state();

  }

  Rcpp::Rcout << "Obs   score at " << t << ": " << obs.t()
              << "State score at " << t << ": " << state.t();

}

std::vector<std::set<arma::uword> > get_ancestors
  (const std::vector<cloud> &clouds)
{
  std::vector<std::set<arma::uword> > out(clouds.size());

  std::vector<std::set<arma::uword> >::reverse_iterator oi = out.rbegin();
  bool is_first = true;
  for(auto ci = clouds.rbegin(); ci != clouds.rend(); ++ci, ++oi){
    if(is_first){
      /* insert all */
      for(auto p : *ci)
        oi->insert(p.get_cloud_idx());

      is_first = false;
      continue;

    }

    /* insert parents */
    const std::set<arma::uword> &prev = *(oi - 1L);
    const cloud &children = *(ci - 1L);
    for(auto o : prev){
      const particle &p = children[o];
      if(p.parent)
        oi->insert(p.parent->get_cloud_idx());

    }
  }

  return out;
}


/* functions and classes to make approximations as in
 * Cappé, Olivier, and Eric Moulines. "Recursive computation of the score and observed information matrix in hidden Markov models." In IEEE/SP 13th Workshop on Statistical Signal Processing, 2005, pp. 703-708. IEEE, 2005. */

struct score_n_hess_dat;
class score_n_hess : public score_n_hess_base {
  arma::vec a_state;
  arma::vec a_obs;
  arma::mat B_state;
  arma::mat B_obs;
  bool is_set;
  double weight;
public:
  const arma::mat &get_a_state() const override;
  const arma::mat &get_a_obs() const override;
  const arma::mat &get_B_state() const override;
  const arma::mat &get_B_obs() const override;
  const double get_weight() const override;

  score_n_hess();
  score_n_hess(const score_n_hess_dat&, const particle&, const particle&,
               const bool);

  score_n_hess& operator+=(const score_n_hess&);
};

struct score_n_hess_dat {
  const arma::mat X;
  const arma::vec eta;
  const arma::vec y;
  const arma::vec dts;
  const arma::mat ran_vars;
  std::unique_ptr<family_base> family;
  const arma::mat &F;
  const arma::mat &Q;
  const arma::mat Q_chol;
  const arma::mat K;
  const arma::mat K_3_2;

  score_n_hess_dat(arma::mat &&X, arma::vec &&y, arma::vec &&dts,
                   arma::mat &&ran_vars, const arma::vec &fixed_params,
                   const std::string family, const arma::mat &F,
                   const arma::mat &Q):
    X(X), eta(X.t() * fixed_params), y(y), dts(dts), ran_vars(ran_vars),
    family(get_fam<family_base>(family)), F(F), Q(Q), Q_chol(arma::chol(Q)),
    K(Q.i()), K_3_2(K * 1.5)
  { }
};

score_n_hess_dat get_score_n_hess_dat
  (const arma::mat &ran_vars, const arma::mat &fixed_terms,
   const arma::ivec &is_event_in, const arma::vec &event_times,
   const arma::vec &tstart, const arma::vec &tstop,
   const arma::uvec &r_obj, const int ti, const arma::vec &fixed_params,
   const std::string family, const arma::mat &F, const arma::mat &Q)
{
  arma::mat ran_vars_i = ran_vars   .cols(r_obj);
  arma::mat X_i        = fixed_terms.cols(r_obj);
  arma::vec y_i;
  {
    arma::ivec tmp = is_event_in.elem(r_obj);
    tmp.for_each([ti](arma::ivec::elem_type &val) { val = val == (int)ti; } );
    y_i = arma::conv_to<arma::vec>::from(tmp);
  }

  arma::vec dts;
  {
    double int_stop  = event_times[ti + 1L],
           int_start = event_times[ti     ];
    arma::vec sta = tstart.elem(r_obj), sto = tstop.elem(r_obj);
    sta.for_each([int_start](arma::vec::elem_type &val) {
      val = std::max(val, int_start); } );
    sto.for_each([int_stop ](arma::vec::elem_type &val) {
      val = std::min(val, int_stop ); } );
    dts = sto - sta;
  }

  return {
    std::move(X_i), std::move(y_i), std::move(dts),
    std::move(ran_vars_i), fixed_params, family, F, Q
  };
}

const arma::mat &score_n_hess::get_a_state() const {
  return a_state;
}
const arma::mat &score_n_hess::get_a_obs() const {
  return a_obs;
}
const arma::mat &score_n_hess::get_B_state() const {
  return B_state;
}
const arma::mat &score_n_hess::get_B_obs() const {
  return B_obs;
}
const double score_n_hess::get_weight() const {
  return weight;
}


struct state_derivs_output {
  const arma::vec a_state;
  const arma::mat B_state;
};

state_derivs_output get_state_derivs_output
  (const score_n_hess_dat &dat, const particle &p, const particle &parent,
   const bool only_score)
{
  arma::uword k = p.get_state().n_elem,
    B_sta_end_1 = k * k - 1L, B_sta_end_2 = 2L * k * k - 1L;
  arma::vec a_state(k * k * 2L, arma::fill::zeros);
  arma::mat B_state;
  if(!only_score)
    B_state.zeros(k * k * 2L, k * k * 2L);

  /* TODO: move some of this computation to somwhere else (e.g., the kronecker
   *       product)... */

  arma::vec innovation = p.get_state() - dat.F * parent.get_state();
  arma::vec innovation_std = solve_w_precomputed_chol(dat.Q_chol, innovation);
  arma::mat Fd = parent.get_state() * innovation_std.t(),
    Qd = (innovation / 2) * innovation_std.t();
  Qd.diag() -= .5;
  Qd = solve_w_precomputed_chol(dat.Q_chol, Qd);

  double *a_state_it = a_state.begin();
  for(auto f = Fd.begin(); f != Fd.end(); ++f, ++a_state_it)
    *a_state_it = *f;
  for(auto qd = Qd.begin(); qd != Qd.end(); ++qd, ++a_state_it)
    *a_state_it = *qd;

  if(!only_score){
    B_state.submat(0L, 0L, B_sta_end_1, B_sta_end_1) =
    arma::kron(dat.K, (-parent.get_state()) * parent.get_state().t());

    B_state.submat(0L, B_sta_end_1 + 1L, B_sta_end_1, B_sta_end_2).zeros();
    B_state.submat
      (B_sta_end_1 + 1L, B_sta_end_1 + 1L, B_sta_end_2, B_sta_end_2) =
        arma::kron(innovation_std * innovation_std.t() - dat.K_3_2, dat.K);
  }

  return { std::move(a_state), std::move(B_state) };
}

struct obs_derivs_output {
  const arma::vec a_obs;
  const arma::mat B_obs;
};

obs_derivs_output get_obs_derivs_output
  (const score_n_hess_dat &dat, const particle &p, const bool only_score)
{
  arma::uword q = dat.X.n_rows, n = dat.X.n_cols;
  int qi = q;
  arma::vec a_obs(q, arma::fill::zeros);
  arma::mat B_obs;
  if(!only_score)
    B_obs.zeros(q, q);

  for(arma::uword i = 0; i < n; ++i){
    double eta = dat.eta[i] + arma::dot(dat.ran_vars.col(i), p.get_state());

    /* terms from observation equation */
    auto trunc_eta = dat.family->truncate_eta(
      dat.y[i], eta, exp(eta), dat.dts[i]);
    double d  = dat.family->d_log_like (dat.y[i], trunc_eta, dat.dts[i]);
    a_obs += d * dat.X.col(i);

    if(!only_score){
      double dd = dat.family->dd_log_like(dat.y[i], trunc_eta, dat.dts[i]);
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&qi, &dd, dat.X.colptr(i), B_obs.memptr());
    }
  }

  return { std::move(a_obs), std::move(B_obs) };
}

score_n_hess::score_n_hess():
  is_set(false), weight(std::numeric_limits<double>::quiet_NaN()) { }

score_n_hess::score_n_hess
  (const score_n_hess_dat &dat, const particle &p, const particle &parent,
   const bool only_score): is_set(true), weight(exp(p.log_weight))
{
  {
    state_derivs_output o =
      get_state_derivs_output(dat, p, parent, only_score);
    a_state = std::move(o.a_state);
    B_state = std::move(o.B_state);
  }

  {
    obs_derivs_output o =
      get_obs_derivs_output(dat, p, only_score);
    a_obs = std::move(o.a_obs);
    B_obs = std::move(o.B_obs);
  }
}

score_n_hess& score_n_hess::operator+=(const score_n_hess& rhs){
  if(is_set){
    a_state += rhs.a_state;
    a_obs += rhs.a_obs;
    B_state += rhs.B_state;
    B_obs += rhs.B_obs;

  } else {
    a_state = rhs.a_state;
    a_obs = rhs.a_obs;
    B_state = rhs.B_state;
    B_obs = rhs.B_obs;

  }

  return *this;
}

class PF_get_score_n_hess_worker {
public:
  const score_n_hess_dat &dat;
  const particle &p;
  const particle &parent;
  const bool only_score;
  PF_get_score_n_hess_worker(
    const score_n_hess_dat &dat, const particle &p,
    const particle &parent, const bool only_score):
    dat(dat), p(p), parent(parent), only_score(only_score) { }

  score_n_hess operator()(){
    return score_n_hess(dat, p, parent, only_score);
  }
};

std::vector<std::unique_ptr<score_n_hess_base> > PF_get_score_n_hess
  (const std::vector<cloud> &forward_clouds, const arma::mat &Q,
   const arma::mat &F, const std::vector<arma::uvec> &risk_obj,
   const arma::ivec &is_event_in, const arma::vec &event_times,
   const arma::mat &ran_vars, const arma::mat &fixed_terms,
   const arma::vec &tstart, const arma::vec &tstop,
   const arma::vec &fixed_params, const std::string family,
   const int max_threads, const bool debug, const bool only_score)
{
  thread_pool pool(std::max((int)1L, max_threads));
  std::vector<std::set<arma::uword> > needed_particles =
    get_ancestors(forward_clouds);

  std::map<arma::uword, score_n_hess> old_res;
  auto particle_idxs = needed_particles.begin() + 1L;
  auto cl_i = forward_clouds.begin() + 1L;
  auto r_obj = risk_obj.begin();
  bool is_first = true;

  unsigned int i = 0;
  for(; particle_idxs != needed_particles.end();
      ++particle_idxs, ++r_obj, ++cl_i, ++i){
    score_n_hess_dat dat = get_score_n_hess_dat(
      ran_vars, fixed_terms, is_event_in, event_times, tstart, tstop,
      *r_obj, i, fixed_params, family, F, Q);

    /* setup workers to perform computation */
    unsigned int n_parts = particle_idxs->size();
    std::vector<std::future<score_n_hess> > futures;
    futures.reserve(n_parts);
    const cloud &cl = *cl_i;

    for(auto j : *particle_idxs){
      PF_get_score_n_hess_worker wo(dat, cl[j], *cl[j].parent, only_score);
      futures.push_back(pool.submit(std::move(wo)));
    }

    std::map<arma::uword, score_n_hess> new_res;
    auto p_idx = particle_idxs->begin();
    for(unsigned int j = 0; j < n_parts; ++j, ++p_idx){
      new_res[*p_idx] = futures[j].get();

      if(is_first)
        continue;

      new_res[*p_idx] += old_res[cl[*p_idx].parent->get_cloud_idx()];
    }

    is_first = false;
    old_res = std::move(new_res);

    if(debug){
      /* means here are slightly odd as we condition on being an ancestor
       * in the end... */
      std::vector<const score_n_hess_base*> to_print_list;
      to_print_list.reserve(old_res.size());
      for(auto x : old_res)
        to_print_list.emplace_back(&x.second);

      print_means_score_n_hess(to_print_list, i + 1L);

    }
  }

  std::vector<std::unique_ptr<score_n_hess_base> > out;
  out.reserve(old_res.size());
  for(auto x : old_res)
    out.emplace_back(new score_n_hess(std::move(x.second)));

  return out;
}



/* functions and classes to make approximations as in
 * Poyiadjis, George, Arnaud Doucet, and Sumeetpal S. Singh. "Particle approximations of the score and observed information matrix in state space models with application to parameter estimation." Biometrika 98, no. 1 (2011): 65-80. */

class score_n_hess_O_N_sq : public score_n_hess_base {
  arma::vec a_state;
  arma::vec a_obs;
  arma::mat B_state;
  arma::mat B_obs;
  double weight;
public:

  const arma::mat &get_a_state() const override;
  const arma::mat &get_a_obs() const override;
  const arma::mat &get_B_state() const override;
  const arma::mat &get_B_obs() const override;
  const double get_weight() const override;

  score_n_hess_O_N_sq() = default;
  score_n_hess_O_N_sq(
    const score_n_hess_dat&, const particle&,
    const cloud&, const std::vector<double>&,
    const std::vector<score_n_hess_O_N_sq>, const bool);
};

const double neg_one = -1;

score_n_hess_O_N_sq::score_n_hess_O_N_sq(
  const score_n_hess_dat &dat, const particle &p,
  const cloud &old_cl, const std::vector<double> &ws,
  const std::vector<score_n_hess_O_N_sq> old_score_n_hess,
  const bool only_score): weight(exp(p.log_weight))
{
  /* is it first iteration? */
  const bool is_first_it = old_score_n_hess.size() < 1L;

  {
    obs_derivs_output o =
      get_obs_derivs_output(dat, p, only_score);
    a_obs = std::move(o.a_obs);
    B_obs = std::move(o.B_obs);
  }

  const arma::vec a_obs_org = a_obs;

  const int dim_state = p.get_state().n_elem * p.get_state().n_elem * 2L,
    dim_obs = a_obs.n_elem;
  a_state.zeros(dim_state);
  if(!only_score)
    B_state.zeros(dim_state, dim_state);

  auto w = ws.begin();
  auto old_res = old_score_n_hess.begin();
  for(const auto parent : old_cl){
    double w_i = exp(*(w++));

    state_derivs_output o =
      get_state_derivs_output(dat, p, parent, only_score);
    if(!is_first_it){
      arma::vec term_obs   = a_obs_org + old_res->get_a_obs(),
                term_state = o.a_state + old_res->get_a_state();

      a_state += w_i * term_state;
      a_obs   += w_i * old_res->get_a_obs();

      if(!only_score){
        B_state += w_i * (old_res->get_B_state() + o.B_state);
        B_obs   += w_i *  old_res->get_B_obs();

        /* add extra outer cross-product terms */
        R_BLAS_LAPACK::sym_mat_rank_one_update
          (&dim_obs  , &w_i, term_obs.memptr()  , B_obs.memptr());
        R_BLAS_LAPACK::sym_mat_rank_one_update
          (&dim_state, &w_i, term_state.memptr(), B_state.memptr());

      }

      ++old_res;
      continue;

    }

    /* first iteration */
    a_state += w_i * o.a_state;

    if(!only_score){
      B_state += w_i * o.B_state;

      /* add extra outer cross-product terms */
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&dim_state, &w_i, o.a_state.memptr(), B_state.memptr());

    }
  }

  if(!only_score){
    if(!is_first_it)
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&dim_obs, &neg_one, a_obs.memptr(), B_obs.memptr());
    R_BLAS_LAPACK::sym_mat_rank_one_update
      (&dim_state, &neg_one, a_state.memptr(), B_state.memptr());

  }
}

const arma::mat &score_n_hess_O_N_sq::get_a_state() const {
  return a_state;
}
const arma::mat &score_n_hess_O_N_sq::get_a_obs() const {
  return a_obs;
}
const arma::mat &score_n_hess_O_N_sq::get_B_state() const {
  return B_state;
}
const arma::mat &score_n_hess_O_N_sq::get_B_obs() const {
  return B_obs;
}
const double score_n_hess_O_N_sq::get_weight() const {
  return weight;
}

#define PF_GET_SCORE_N_HESS_O_N_SQ_ARGS                                   \
arma::mat &Q, const arma::mat &F,                                         \
const std::vector<arma::uvec> &risk_obj, const Rcpp::List &risk_obj_R,    \
const arma::ivec &is_event_in, const arma::vec &event_times,              \
arma::mat &ran_vars, arma::mat &fixed_terms,                              \
const arma::vec &tstart, const arma::vec &tstop,                          \
const arma::vec &fixed_params, const std::string family,                  \
const int max_threads, const bool debug, const bool only_score,           \
const arma::vec &a_0, const arma::mat &R, arma::mat &Q_0,                 \
const arma::mat &Q_tilde, const arma::uword N_fw_n_bw,                    \
const arma::uword N_first, const double nu, const double covar_fac,       \
const double ftol_rel,                                                    \
Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold

#define PF_GET_SCORE_N_HESS_O_N_SQ_INPUT                                  \
  Q, F, risk_obj, risk_obj_R, is_event_in, event_times, ran_vars,         \
  fixed_terms, tstart, tstop, fixed_params, family, max_threads,          \
  debug, only_score, a_0, R, Q_0, Q_tilde, N_fw_n_bw, N_first,            \
  nu, covar_fac, ftol_rel, forward_backward_ESS_threshold

template<template <bool> class T_resampler>
std::vector<std::unique_ptr<score_n_hess_base> > PF_get_score_n_hess_O_N_sq_comp
  (PF_GET_SCORE_N_HESS_O_N_SQ_ARGS)
{
#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#endif

  cloud old_cl, new_cl;
  const PF_data data
    (0L /* n_fixed_terms_in_state_vec */, ran_vars, fixed_terms, tstart,
     tstop, is_event_in, a_0, R, R.t(), Q_0, Q, risk_obj_R, F, 1L /* n_max*/,
     max_threads, fixed_params, Q_tilde, N_fw_n_bw, 1L /* N_smooth */,
     1L /* N_smooth_final */, forward_backward_ESS_threshold,
     debug ? 2L : 0L, N_first, nu, covar_fac, ftol_rel);
  std::string family_use = get_family(family);
  pf_dens dens_calc(data, family_use);

  /* sample at time zero */
  old_cl =
    importance_dens_normal_approx_w_cloud_mean_independent
    <true>::sample_first_state_n_set_weights(dens_calc, data);
  auto r_obj = risk_obj.begin();
  std::vector<score_n_hess_O_N_sq> old_n_hess_O_N_sqs, new_n_hess_O_N_sqs;

  if(data.debug > 0)
    data.log(1) << "Cloud means at time zero are" << std::endl
                << old_cl.get_weigthed_mean().t();

  /* run particle filter and compute smoothed functionals recursively */
  for(unsigned int t = 1;  t <= (unsigned int)data.d; ++t, ++r_obj){
    if((t + 1L) % 3L == 0L)
      Rcpp::checkUserInterrupt();

    std::shared_ptr<PF_cdist> y_dist = dens_calc.get_y_dist(t);

    if(data.debug > 0)
      data.log(1) << "Starting iteration " << t;

    arma::uvec resample_idx;
    bool did_resample;
    T_resampler<true>::resampler(
      dens_calc, data, old_cl, y_dist, t, resample_idx, did_resample);

    if(data.debug > 0){
      if(did_resample)
        data.log(1) << "Did resample";
      else
        data.log(1) << "Did not re-sample";

      data.log(1) << "Sampling";

    }

    if(did_resample){
      /* get new weights and indices to keep */
      std::map<arma::uword, double> new_idx_n_weights =
        get_resample_idx_n_log_weight
      <cloud, get_weight_from_particle, get_resample_weight_from_particle>
       (old_cl, resample_idx);

      /* remove non-sampled particles and update weights */
      cloud old_cl_resampled;
      std::vector<score_n_hess_O_N_sq> old_n_hess_O_N_sqs_resampled;
      old_cl_resampled.reserve(new_idx_n_weights.size());
      old_n_hess_O_N_sqs_resampled.reserve(new_idx_n_weights.size());

      for(auto it : new_idx_n_weights)
      {
        old_cl_resampled.new_particle(std::move(old_cl[it.first]));
        old_cl_resampled.back().log_weight = it.second;

        if(old_n_hess_O_N_sqs.size() > 0L)
          /* move only if it is not the first iteration */
          old_n_hess_O_N_sqs_resampled.emplace_back(
            std::move(old_n_hess_O_N_sqs[it.first]));
      }

      old_cl = std::move(old_cl_resampled);
      old_n_hess_O_N_sqs = std::move(old_n_hess_O_N_sqs_resampled);
    }

    new_cl = importance_dens_normal_approx_w_cloud_mean_independent
      <true>::sample(
        y_dist, dens_calc, data, old_cl, resample_idx, t, nothing());

    if(data.debug > 0)
      data.log(1) << "Setting weights";

    /* N x N log weights matrix */
    unsigned int n_elem = new_cl.size();
    std::vector<std::vector<double> > ws(n_elem);
    {
      double max_weight = -std::numeric_limits<double>::max();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
      for(unsigned int j = 0; j < n_elem; ++j){ // loop over new particles
        auto it = new_cl.begin() + j;
        it->log_weight =
          y_dist->log_dens(it->get_state()) - it->log_importance_dens;

        /* account for state transition */
        std::vector<double> &part_ws = ws[j];
        part_ws.resize(old_cl.size());
        auto part_ws_j = part_ws.begin();
        double max_weight_inner = -std::numeric_limits<double>::max();
        for(const auto parent : old_cl){ /* loop over parents */
          *(part_ws_j) =
            parent.log_weight +
            dens_calc.log_prob_state_given_parent(
              it->get_state(), parent.get_state());
          max_weight_inner = MAX(max_weight_inner, *(part_ws_j));
          ++part_ws_j;

        }

        /* normalize "inner" weights */
        {
          normalize_weights_output tmp = normalize_log_weights<false, true>
          (part_ws, max_weight_inner);
          it->log_weight += tmp.log_sum_logs;
        }

        max_weight = MAX(max_weight, it->log_weight);
      }

      // normalize weights
      if(data.debug > 0){
        normalize_weights_output tmp =
          normalize_log_weights<true, true>(new_cl, max_weight);

        data.log(1) << "ESS is at time " << t << " is " << tmp.ESS
                    << " and cloud mean is" << std::endl
                    << new_cl.get_weigthed_mean().t();

      } else
        normalize_log_weights<false, true>(new_cl, max_weight);

    }

    /* compute smoothed functionals */
    score_n_hess_dat dat = get_score_n_hess_dat(
      ran_vars, fixed_terms, is_event_in, event_times, tstart, tstop,
      *r_obj, t - 1L, fixed_params, family, F, Q);

    if(data.debug > 0)
      data.log(1) << "Computing smoothed functionals";

    /* TODO: use an openMP reduction */
    new_n_hess_O_N_sqs.clear();
    new_n_hess_O_N_sqs.resize(n_elem);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(unsigned int j = 0; j < n_elem; ++j){ // loop over new particles
      new_n_hess_O_N_sqs[j] = score_n_hess_O_N_sq(
        dat, *(new_cl.begin() + j), old_cl, ws[j], old_n_hess_O_N_sqs,
        only_score);

    }

    old_n_hess_O_N_sqs = std::move(new_n_hess_O_N_sqs);
    old_cl = std::move(new_cl);

    if(data.debug > 0){
      std::vector<const score_n_hess_base*> to_print_list;
      to_print_list.reserve(old_n_hess_O_N_sqs.size());
      for(auto x : old_n_hess_O_N_sqs)
        to_print_list.emplace_back(&x);

      print_means_score_n_hess(to_print_list, t);

    }
  }

  std::vector<std::unique_ptr<score_n_hess_base> > out;
  out.reserve(old_n_hess_O_N_sqs.size());
  for(auto x : old_n_hess_O_N_sqs)
    out.emplace_back(new score_n_hess_O_N_sq(std::move(x)));

  return out;
}



std::vector<std::unique_ptr<score_n_hess_base> > PF_get_score_n_hess_O_N_sq
  (PF_GET_SCORE_N_HESS_O_N_SQ_ARGS, const std::string method)
{
  if(method == BOOT_FILTER or method == PF_APPROX_CLOUD_MEAN or
       method == PF_APPROX_PARTICLE)
    return PF_get_score_n_hess_O_N_sq_comp
    <None_AUX_resampler>(PF_GET_SCORE_N_HESS_O_N_SQ_INPUT);
  if(method == AUX_APPROX_CLOUD_MEAN)
    return PF_get_score_n_hess_O_N_sq_comp
    <AUX_resampler_normal_approx_w_cloud_mean>(PF_GET_SCORE_N_HESS_O_N_SQ_INPUT);
  if(method == AUX_APPROX_PARTICLE)
    return PF_get_score_n_hess_O_N_sq_comp
    <AUX_resampler_normal_approx_w_particles>(PF_GET_SCORE_N_HESS_O_N_SQ_INPUT);

  Rcpp::stop("PF_get_score_n_hess_O_N_sq: Method not implemented");
}
