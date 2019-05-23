#include "cond_approx.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"

namespace {
  struct mode_object_result {
    const double objective;
    const arma::vec grad;
    const arma::mat hess;
  };
}

inline mode_object_result mode_objective
  (const arma::vec &state, const std::vector<PF_cdist*> &dists,
   const bool do_grad_n_hess){
  arma::vec grad;
  arma::mat hess;
  const arma::uword dim = state.n_elem;

  if(do_grad_n_hess){
    grad.zeros(dim);
    hess.zeros(dim, dim);

  }

  double o = 0.;
  for(auto &d : dists){
    o += d->log_dens(state);

    if(do_grad_n_hess){
      grad += d->gradient(state);
      hess -= d->neg_Hessian(state);
    }
  }

  return { o, std::move(grad), std::move(hess) };
}

cdist_comb_generator::cdist_comb_generator
  (std::vector<PF_cdist*> &cdists, const int nu, const arma::mat *xtra_covar,
   const double covar_fac, const double ftol_rel):
  /* all should have the same dimension */
  cdist_comb_generator(
    cdists, arma::vec(cdists[0]->dim(), arma::fill::zeros), nu, xtra_covar,
    covar_fac, ftol_rel) { }

cdist_comb_generator::cdist_comb_generator(
  std::vector<PF_cdist*> &cdists, const arma::vec &start, const int nu,
  const arma::mat *xtra_covar, const double covar_fac, const double ftol_rel):
  cdists(cdists), nu(nu)
{
  std::vector<bool> is_mvn(cdists.size());
  std::transform(cdists.begin(), cdists.end(), is_mvn.begin(),
                 [](PF_cdist *p) { return p->is_mvn(); });

  bool all_mvn =
    std::all_of(is_mvn.begin(), is_mvn.end(), [](bool i){ return i; });

  /* all should have the same dimension */
  unsigned int n = cdists[0]->dim();
  arma::vec val = start;

  mode_object_result derivs_info = ([&]{
    if(!all_mvn){
      /* find mode. TODO: replace with penalized iteratively reweighted least
      *                  squares */
      static constexpr unsigned i_max = 500L;
      for(unsigned i = 0; i < 500L; ++i){
        /* get original objective value, gradient, and Hessian estimate */
        const mode_object_result start_info = mode_objective(val, cdists, true);

        /* find direction */
        const arma::vec direction = arma::solve(
          start_info.hess, -start_info.grad);

        /* step-half to find solution */
        double step_size = 1.;
        static constexpr unsigned j_max = 20L;
        unsigned j = 0;
        for(; j < j_max; ++j, step_size /= 2.){
          const arma::vec new_val = val + step_size * direction;
          const mode_object_result new_info =
            mode_objective(new_val, cdists, false);

          const double obj_diff = new_info.objective - start_info.objective;
          if(obj_diff > 0.){
            /* found better value. Check if have converged and update final
             * parameter */
            val = std::move(new_val);

            if(obj_diff / (std::abs(new_info.objective) + 1e-16) < ftol_rel)
              return mode_objective(val, cdists, true);

            break;
          }
        }

        if(j >= j_max or i + 1L >= i_max)
          /* should not reach this point */
          return start_info;
      }
    }

    /* we do not need to do any optimization */
    return mode_objective(start, cdists, true);
  })();

  /* compute negative Hessian */
  std::vector<arma::mat> neg_Ks(cdists.size());
  std::transform(cdists.begin(), cdists.end(), neg_Ks.begin(),
                 [&](PF_cdist* p) { return p->neg_Hessian(val); });

  neg_K.set_size(n, n);
  neg_K.zeros();
  for(auto K = neg_Ks.begin(); K != neg_Ks.end(); ++K)
    neg_K += *K;

  arma::mat Sig_mat_use = (xtra_covar) ?
    arma::mat(neg_K.i() + *xtra_covar) : arma::mat(neg_K.i());
  if(covar_fac > 0)
    Sig_mat_use *= covar_fac;

  covarmat Sig_obj = (nu > 2L) ?
    covarmat(Sig_mat_use * (nu - 2.) / nu) :
    covarmat(Sig_mat_use);

  Sig = std::make_shared<covarmat>(std::move(Sig_obj));

  /* mean part from non-Gaussian conditional distributions */
  k.set_size(n);
  k.zeros();
  for(unsigned int i = 0; i < cdists.size(); ++i){
    if(cdists[i]->is_mvn())
      continue;
    k += neg_Ks[i] * val + cdists[i]->gradient(val);
  }
}

class cdist_comb final : public dist_comb {
  std::shared_ptr<covarmat> Sig;
  const int nu;
  arma::vec mu;

public:
  cdist_comb(const std::initializer_list<arma::vec*>&,
             std::vector<PF_cdist*>&, const arma::mat&,
             const arma::vec&, std::shared_ptr<covarmat>, const int);
  ~cdist_comb() = default;

  arma::vec sample() const override;
  double log_density(const arma::vec&) const override;
  const arma::vec& get_mean() const override;
  const arma::mat& get_covar() const override;
};

std::unique_ptr<dist_comb> cdist_comb_generator::get_dist_comb(
    const std::initializer_list<arma::vec*> &states)
{
  return std::unique_ptr<dist_comb>(new cdist_comb(
      states, cdists, neg_K, k, Sig, nu));
}

cdist_comb::cdist_comb(
  const std::initializer_list<arma::vec*> &states,
  std::vector<PF_cdist*> &cdists, const arma::mat &K, const arma::vec &k,
  std::shared_ptr<covarmat> Sig, const int nu): Sig(Sig), nu(nu)
{
  auto s = states.begin();
  mu = k;
  for(auto p = cdists.begin(); p != cdists.end(); ++p){
    if(!(*p)->is_mvn())
      continue;
    if((*p)->is_grad_z_hes_const()){
      mu += (*p)->gradient_zero(nullptr);
      continue;
    }

    mu += (*p)->gradient_zero(*(s++));
  }

  /* TODO: may be expensive if called multiple times. Make decompostion
   *       or use Sig, Do not do the latter as we change sig if we use
   *       a multivariate t-distribution... */
  mu = arma::solve(K, mu);
}

arma::vec cdist_comb::sample() const {
  if(nu >= 2L)
    return mvtrnorm(mu, Sig->chol(), nu);
  return mvrnorm(mu, Sig->chol());
}

double cdist_comb::log_density(const arma::vec &state) const {
  if(nu >= 2L)
    return dmvtrm_log(state, mu, Sig->chol_inv(), nu);
  return dmvnrm_log(state, mu, Sig->chol_inv());
}

const arma::vec& cdist_comb::get_mean() const {
  return mu;
}

const arma::mat& cdist_comb::get_covar() const {
  return Sig->mat();
}
