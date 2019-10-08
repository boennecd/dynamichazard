#include "PFs.h"
#include "../thread_pool.h"
#include "importance_samplers.h"
#include "resamplers.h"

constexpr static char C_U = 'U';
constexpr static int I_ONE = 1L;
static const double D_NEG_ONE = -1;

#define MAX(a,b) (((a)>(b))?(a):(b))

inline void print_means_score_n_hess
  (const std::vector<const score_n_hess_base*> &dat, const unsigned int t,
   const bool has_hess){
  arma::vec score;
  arma::mat hess_terms;

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
      score   = *w * i->get_score();
      if(has_hess)
        hess_terms  = *w * i->get_hess_terms();
      first = false;
      continue;

    }

    score    += *w * i->get_score();
    if(has_hess)
      hess_terms   += *w * i->get_hess_terms();

  }

  Rcpp::Rcout << "Score at " << t << '\n' << score.t();

  if(has_hess)
    Rcpp::Rcout
              << "Hessian terms at " << t << '\n' << hess_terms;


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

/* object to compute the score and Hessian for a given particle up to a
 * given point in time */
struct score_n_hess_dat;
class score_n_hess : public score_n_hess_base {
  arma::vec score;
  arma::mat hess_terms;
  bool is_set;
  double weight;
public:
  const arma::vec &get_score() const override {
    return score;
  }
  const arma::mat &get_hess_terms() const override {
    return hess_terms;
  }
  double get_weight() const override {
    return weight;
  }

  score_n_hess();
  score_n_hess(const score_n_hess_dat&, const particle&, const particle&,
               const bool);

  score_n_hess& operator+=(const score_n_hess&);
};

/* object used to store the data needed in computation of the gradient and
 * Hessian terms */
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
  const arma::mat K_1_2;

  score_n_hess_dat(arma::mat &&X, arma::vec &&y, arma::vec &&dts,
                   arma::mat &&ran_vars, const arma::vec &fixed_params,
                   const std::string family, const arma::mat &F,
                   const arma::mat &Q):
    X(X), eta(X.t() * fixed_params), y(y), dts(dts), ran_vars(ran_vars),
    family(get_fam<family_base>(family)), F(F), Q(Q), Q_chol(arma::chol(Q)),
    K(Q.i()), K_1_2(K * .5)
  { }
};

/* function to create a type of the above for a given point in time */
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

struct derivs_output {
  const arma::vec score;
  const arma::mat hess_terms;
};

/* function to compute the score and Hessian terms */
derivs_output get_derivs_output
  (const score_n_hess_dat &dat, const particle &p, const particle &parent,
   const bool only_score) {
  const int
    dfixd     = dat.X.n_rows,
    dsate     = p.get_state().n_elem, dsatesq = dsate * dsate,
    score_dim = dfixd + 2L * dsatesq;
  const arma::uword n = dat.X.n_cols;

  /* setup score and Hessian terms */
  arma::vec score(score_dim, arma::fill::zeros);
  arma::mat hess_terms =
    only_score ?
    arma::mat() : arma::mat(score_dim, score_dim, arma::fill::zeros);

  /* we will only compute the upper part of the Hessian. We copy it to the
   * lower part later. Start with terms from observation equation */
  {
    const arma::span obs_span(0L, dfixd - 1L);
    for(arma::uword i = 0; i < n; ++i){
      double eta = dat.eta[i] + arma::dot(dat.ran_vars.col(i), p.get_state());

      /* terms from observation equation */
      auto trunc_eta = dat.family->truncate_eta(
        dat.y[i], eta, exp(eta), dat.dts[i]);
      const double d  =
        dat.family->d_log_like(dat.y[i], trunc_eta, dat.dts[i]);
      score(obs_span) += d * dat.X.col(i);

      if(!only_score){
        const double dd =
          dat.family->dd_log_like(dat.y[i], trunc_eta, dat.dts[i]);
        R_BLAS_LAPACK::dsyr(
            &C_U, &dfixd, &dd, dat.X.colptr(i), &I_ONE, hess_terms.memptr(),
            &score_dim);
      }
    }
  }

  /* then terms from state equation */
  {
    double * const score_state = score.memptr() + dfixd;

    /* handle dL/dF */
    arma::vec innovation = p.get_state() - dat.F * parent.get_state();
    arma::vec innovation_Qinv =
      solve_w_precomputed_chol(dat.Q_chol, innovation);
    arma::mat Fd(score_state, dsate, dsate, false);
    Fd = parent.get_state() * innovation_Qinv.t();

    /* handle dL/dQ */
    arma::mat Qd(score_state + dsatesq, dsate, dsate, false);
    Qd = (innovation / 2.) * innovation_Qinv.t();
    Qd.diag() -= .5;
    Qd = solve_w_precomputed_chol(dat.Q_chol, Qd);

    if(!only_score){
      const arma::span Fspan(dfixd          , dfixd + dsatesq - 1L),
                       Qspan(dfixd + dsatesq, dfixd + 2L * dsatesq - 1L);
      hess_terms(Fspan, Fspan) =
        arma::kron(
          dat.K, (-parent.get_state()) * parent.get_state().t());
      hess_terms(Fspan, Qspan) =
        arma::kron(
          dat.K, (-parent.get_state()) * innovation_Qinv.t());
      hess_terms(Qspan, Qspan) =
          arma::kron(
            dat.K, (-innovation_Qinv) * innovation_Qinv.t() + dat.K_1_2);
    }
  }

  /* not needed as we later copy the upper part to the lower part */
  // if(!only_score)
  //   hess_terms = arma::symmatu(hess_terms);

  return { std::move(score), std::move(hess_terms) };
}

score_n_hess::score_n_hess():
  is_set(false), weight(std::numeric_limits<double>::quiet_NaN()) { }

score_n_hess::score_n_hess
  (const score_n_hess_dat &dat, const particle &p, const particle &parent,
   const bool only_score): is_set(true), weight(exp(p.log_weight))
{
  auto o = get_derivs_output(dat, p, parent, only_score);
  score = std::move(o.score);
  hess_terms = std::move(o.hess_terms);
}

score_n_hess& score_n_hess::operator+=(const score_n_hess& rhs){
  if(is_set){
    score      += rhs.score;
    hess_terms += rhs.hess_terms;

  } else {
    score      = rhs.score;
    hess_terms = rhs.hess_terms;

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

  unsigned i = 0;
  for(; particle_idxs != needed_particles.end();
      ++particle_idxs, ++r_obj, ++cl_i, ++i){
    score_n_hess_dat dat = get_score_n_hess_dat(
      ran_vars, fixed_terms, is_event_in, event_times, tstart, tstop,
      *r_obj, i, fixed_params, family, F, Q);

    /* setup workers to perform computation */
    const unsigned n_parts = particle_idxs->size();
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

      print_means_score_n_hess(to_print_list, i + 1L, !only_score);

    }
  }

  std::vector<std::unique_ptr<score_n_hess_base> > out;
  out.reserve(old_res.size());
  for(auto x : old_res)
    out.emplace_back(new score_n_hess(std::move(x.second)));

  return out;
}



/* functions and classes to make approximations suggested in
 * Poyiadjis, George, Arnaud Doucet, and Sumeetpal S. Singh. "Particle approximations of the score and observed information matrix in state space models with application to parameter estimation." Biometrika 98, no. 1 (2011): 65-80. */

namespace {
  /* struct used to save some computations in the member function in next
   * class */
  class extended_particle {
    const particle &p;
  public:
    /* stores Q^{-1}Fp */
    const arma::vec QiFp;

    extended_particle(const particle &p, const score_n_hess_dat &dat):
      p(p), QiFp(([&]{
        arma::vec out = dat.F * p.get_state();
        out = solve_w_precomputed_chol(dat.Q_chol, out);
        return out;
      })()) { }

    const arma::vec& get_state() const {
      return p.get_state();
    }
  };

  /* class to perform computation of the aforementioned algorithm */
  class score_n_hess_O_N_sq : public score_n_hess_base {
    arma::vec score;
    arma::mat hess_terms;
    double weight;
  public:
    const arma::vec &get_score() const override {
      return score;
    }
    const arma::mat &get_hess_terms() const override {
      return hess_terms;
    }
    double get_weight() const override {
      return weight;
    }

    score_n_hess_O_N_sq() = default;
    score_n_hess_O_N_sq(
      const score_n_hess_dat&, const particle&,
      const std::vector<extended_particle>&, const std::vector<double>&,
      const std::vector<score_n_hess_O_N_sq>&, const bool);
  };
}

score_n_hess_O_N_sq::score_n_hess_O_N_sq(
  const score_n_hess_dat &dat, const particle &p,
  const std::vector<extended_particle> &old_cl, const std::vector<double> &ws,
  const std::vector<score_n_hess_O_N_sq> &old_score_n_hess,
  const bool only_score): weight(exp(p.log_weight))
{
  /* is it first iteration? */
  const bool is_first_it = old_score_n_hess.size() < 1L;

  const int
      dfixd     = dat.X.n_rows,
      dsate     = p.get_state().n_elem, dsatesq = dsate * dsate,
      score_dim = dfixd + 2L * dsatesq;
  const arma::uword n = dat.X.n_cols;

  /* setup score and Hessian terms */
  score.zeros(score_dim);
  hess_terms =
    only_score ?
    arma::mat() : arma::mat(score_dim, score_dim, arma::fill::zeros);

  /* compute terms from observation equation */
  const arma::span obs_span(0L, dfixd - 1L), sta_span(dfixd, score_dim - 1L);
  {
    for(arma::uword i = 0; i < n; ++i){
      double eta = dat.eta[i] + arma::dot(dat.ran_vars.col(i), p.get_state());

      /* terms from observation equation */
      auto trunc_eta = dat.family->truncate_eta(
        dat.y[i], eta, exp(eta), dat.dts[i]);
      const double d =
        dat.family->d_log_like (dat.y[i], trunc_eta, dat.dts[i]);
      score(obs_span) += d * dat.X.col(i);

      if(!only_score){
        const double dd =
          dat.family->dd_log_like(dat.y[i], trunc_eta, dat.dts[i]);
        R_BLAS_LAPACK::dsyr(
            &C_U, &dfixd, &dd, dat.X.colptr(i), &I_ONE, hess_terms.memptr(),
            &score_dim);
      }
    }
  }

  /* take a copy of the gradient w.r.t. the terms from the observation equation */
  const arma::vec obs_score_term =
    only_score ? arma::vec() : score(obs_span);

  /* compute terms from state equation and handle additional outer product
   * terms with this algorithm */
  const arma::span Fspan(dfixd          , dfixd + dsatesq - 1L),
                   Qspan(dfixd + dsatesq, dfixd + 2L * dsatesq - 1L);
  {
    /* setup objects which we will use to store score terms */
    arma::vec
      score_terms   (score_dim, arma::fill::zeros);
    arma::mat
      score_terms_Fd(score_terms.memptr() + dfixd          , dsate, dsate, false),
      score_terms_Qd(score_terms.memptr() + dfixd + dsatesq, dsate, dsate, false);

    auto w = ws.begin();
    auto old_res = old_score_n_hess.begin();
    const arma::vec Qip = solve_w_precomputed_chol(dat.Q_chol, p.get_state());

    for(auto parent = old_cl.cbegin(); parent !=  old_cl.cend();
        ++w, ++old_res, ++parent){
      double w_i = exp(*w);

      /* handle dL/dF */
      /* Q^{-1}(y - F x) */
      const arma::vec innovation_Qinv = Qip - parent->QiFp;
      score_terms_Fd = parent->get_state() * innovation_Qinv.t();

      /* handle dL/dQ */
      score_terms_Qd =
        (innovation_Qinv / 2) * innovation_Qinv.t() - dat.K_1_2;

      /* add old score terms */
      if(!is_first_it){
        score_terms(obs_span)  = old_res->get_score()(obs_span);
        score_terms(sta_span) += old_res->get_score()(sta_span);
      } else
        score_terms(obs_span).zeros();

      /* add terms to score */
      score += w_i * score_terms;

      if(!only_score){
        /* add terms Hessian given state variables */
        hess_terms(Fspan, Fspan) +=
          arma::kron(
            dat.K, ((-w_i) * parent->get_state()) * parent->get_state().t());
        hess_terms(Fspan, Qspan) +=
          arma::kron(
            dat.K, ((-w_i) * parent->get_state()) * innovation_Qinv.t());
        hess_terms(Qspan, Qspan) +=
          arma::kron(
            dat.K, w_i * ((- innovation_Qinv) * innovation_Qinv.t() + dat.K_1_2));

        /* add outer product of score terms from this pair */
        score_terms(obs_span) += obs_score_term;
        R_BLAS_LAPACK::dsyr(
            &C_U, &score_dim, &w_i, score_terms.memptr(), &I_ONE,
            hess_terms.memptr(), &score_dim);

        if(!is_first_it)
          hess_terms += w_i * old_res->get_hess_terms();
      }
    }
  }

  if(!only_score){
   /* subtract outer product of score */
   R_BLAS_LAPACK::dsyr(
     &C_U, &score_dim, &D_NEG_ONE, score.memptr(), &I_ONE, hess_terms.memptr(),
     &score_dim);

   /* not needed as we later copy the upper part to the lower part */
   // hess_terms = arma::symmatu(hess_terms);
  }
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
  for(unsigned t = 1;  t <= (unsigned)data.d; ++t, ++r_obj){
    if(t % 5L == 0L)
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

      for(auto &it : new_idx_n_weights)
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
    const unsigned n_elem = new_cl.size();
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
        double max_weight_inner = -std::numeric_limits<double>::infinity();
        for(const auto &parent : old_cl){ /* loop over parents */
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


    /* setup vector of extended particle class that we will use to save some
     * computation in the next calls */
    const std::vector<extended_particle> ex_old_cl = ([&]{
      std::vector<extended_particle> out;
      out.reserve(old_cl.size());
      for(auto &x : old_cl)
        out.emplace_back(x, dat);

      return out;
    })();

    new_n_hess_O_N_sqs.clear();
    new_n_hess_O_N_sqs.resize(n_elem);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(unsigned j = 0; j < n_elem; ++j){ // loop over new particles
      new_n_hess_O_N_sqs[j] = score_n_hess_O_N_sq(
        dat, *(new_cl.begin() + j), ex_old_cl, ws[j], old_n_hess_O_N_sqs,
        only_score);

    }

    old_n_hess_O_N_sqs = std::move(new_n_hess_O_N_sqs);
    old_cl = std::move(new_cl);

    if(data.debug > 0){
      std::vector<const score_n_hess_base*> to_print_list;
      to_print_list.reserve(old_n_hess_O_N_sqs.size());
      for(auto &x : old_n_hess_O_N_sqs)
        to_print_list.emplace_back(&x);

      print_means_score_n_hess(to_print_list, t, !only_score);

    }
  }

  std::vector<std::unique_ptr<score_n_hess_base> > out;
  out.reserve(old_n_hess_O_N_sqs.size());
  for(auto &x : old_n_hess_O_N_sqs)
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
