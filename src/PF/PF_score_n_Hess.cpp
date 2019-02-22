#include "PF_utils.h"
#include "../thread_pool.h"

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
  const arma::mat &K;
  const arma::mat &K_3_2;

  score_n_hess_dat(arma::mat &&X, arma::vec &&y, arma::vec &&dts,
                   arma::mat &&ran_vars, const arma::vec &fixed_params,
                   const std::string family, const arma::mat &F,
                   const arma::mat &Q):
    X(X), eta(X.t() * fixed_params), y(y), ran_vars(ran_vars), dts(dts),
    family(get_fam<family_base>(family)), F(F), Q(Q), Q_chol(arma::chol(Q)),
    K(Q.i()), K_3_2(K * 1.5)
    { }
};

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

score_n_hess::score_n_hess(): is_set(false) { }

score_n_hess::score_n_hess
  (const score_n_hess_dat &dat, const particle &p, const particle &parent,
   const bool only_score)
{
  arma::uword k = p.get_state().n_elem, q = dat.X.n_rows, n = dat.X.n_cols,
    B_sta_end_1 = k * k - 1L, B_sta_end_2 = 2L * k * k - 1L;
  int qi = q;
  a_state.zeros(k * k * 2L);
  a_obs.zeros(q);
  if(!only_score){
    B_obs.zeros(q, q);
    B_state.zeros(k * k * 2L, k * k * 2L);

  } else {
    B_obs = arma::mat();
    B_state = arma::mat();

  }

  /* terms from state equation */
  arma::vec innovation = p.get_state() - dat.F * parent.get_state();
  arma::vec innovation_std = solve_w_precomputed_chol(dat.Q_chol, innovation);
  arma::mat Fd = innovation_std * parent.get_state().t(),
            Qd = (innovation / 2) * innovation_std.t();
  Qd.diag() -= .5;
  Qd = solve_w_precomputed_chol(dat.Q_chol, Qd);

  double *a_state_it = a_state.begin();
  for(auto f = Fd.begin(); f != Fd.end(); ++f, ++a_state_it)
    *a_state_it = *f;
  for(auto qd = Qd.begin(); qd != Qd.end(); ++qd, ++a_state_it)
    *a_state_it = *qd;

  if(!only_score){
    /* TODO: cam be done smart? */
    B_state.submat(0L, 0L, B_sta_end_1, B_sta_end_1) =
      arma::kron(dat.K, (-parent.get_state()) * parent.get_state().t());

    B_state.submat(0L, B_sta_end_1 + 1L, B_sta_end_1, B_sta_end_2).zeros();
    B_state.submat
      (B_sta_end_1 + 1L, B_sta_end_1 + 1L, B_sta_end_2, B_sta_end_2) =
        arma::kron(innovation_std * innovation_std.t() - dat.K_3_2, dat.K);
  }

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

  is_set = true;
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

std::vector<score_n_hess> PF_get_score_n_hess
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
    /* prepare data */
    arma::mat ran_vars_i = ran_vars   .cols(*r_obj);
    arma::mat X_i        = fixed_terms.cols(*r_obj);
    arma::vec y_i;
    {
      arma::ivec tmp = is_event_in.elem(*r_obj);
      tmp.for_each([i](arma::ivec::elem_type &val) { val = val == (int)i; } );
      y_i = arma::conv_to<arma::vec>::from(tmp);
    }

    arma::vec dts;
    {
      double int_stop  = event_times[i + 1L],
             int_start = event_times[i     ];
      arma::vec sta = tstart.elem(*r_obj), sto = tstop.elem(*r_obj);
      sta.for_each([int_start](arma::vec::elem_type &val) {
        val = std::max(val, int_start); } );
      sto.for_each([int_stop ](arma::vec::elem_type &val) {
        val = std::min(val, int_stop ); } );
      dts = sto - sta;
    }

    score_n_hess_dat dat(std::move(X_i), std::move(y_i), std::move(dts),
                         std::move(ran_vars_i), fixed_params, family, F, Q);

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
  }

  std::vector<score_n_hess> out;
  out.reserve(old_res.size());
  for(auto x : old_res)
    out.push_back(std::move(x.second));

  return out;
}
