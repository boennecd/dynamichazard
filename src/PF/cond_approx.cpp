#include "cond_approx.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"
#include <nloptrAPI.h>

double mode_objective(
    unsigned int n, const double *x, double *grad, void *data_in)
{
  std::vector<PF_cdist*> *data = (std::vector<PF_cdist*> *) data_in;
  unsigned int dim = (*data->begin())->dim();
  arma::vec state(x, n);

  if (grad) {
    arma::vec gr(grad, n, false);
    gr.zeros();

    for(auto d = data->begin(); d != data->end(); ++d)
      gr += (*d)->gradient(state);
    gr *= -1;
  }

  double o = 0.;
  for(auto d = data->begin(); d != data->end(); ++d)
    o += (*d)->log_dens(state);

  return -o;
}

cdist_comb_generator::cdist_comb_generator(std::vector<PF_cdist*> &cdists,
                                           const int nu):
  /* all should have the same dimension */
  cdist_comb_generator(
    cdists, arma::vec(cdists[0]->dim(), arma::fill::zeros), nu) { }

cdist_comb_generator::cdist_comb_generator(
  std::vector<PF_cdist*> &cdists, const arma::vec &start, const int nu):
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

  if(!all_mvn){
    /* find mode */
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_LBFGS, n);
    nlopt_set_min_objective(opt, mode_objective, &cdists);
    nlopt_set_xtol_rel(opt, 1e-5);

    double minf;
    if (nlopt_optimize(opt, val.memptr(), &minf) < 0){
      nlopt_destroy(opt);
      Rcpp::stop("'nlopt' failed to the mode");
    }

    nlopt_destroy(opt);
  }

  /* compute negative Hessian */
  std::vector<arma::mat> neg_Ks(cdists.size());
  std::transform(cdists.begin(), cdists.end(), neg_Ks.begin(),
                 [&](PF_cdist* p) { return p->neg_Hessian(val); });

  neg_K.set_size(n, n);
  neg_K.zeros();
  for(auto K = neg_Ks.begin(); K != neg_Ks.end(); ++K)
    neg_K += *K;

  covarmat Sig_obj = (nu > 2L) ?
    covarmat(neg_K.i() * (nu - 2.) / nu) :
    covarmat(neg_K.i());

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

class cdist_comb : public dist_comb {
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
