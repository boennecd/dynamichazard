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

cdist_comb_generator::cdist_comb_generator
  (std::vector<PF_cdist*> &cdists, const int nu, const arma::mat *xtra_covar,
   const double covar_fac, const double ftol_rel):
  /* all should have the same dimension */
  cdist_comb_generator(
    cdists, arma::vec(cdists[0]->dim(), arma::fill::zeros), nu, xtra_covar,
    covar_fac, ftol_rel) { }

#define NLOPT_RESULT_CODE_NOT_SET -100L

cdist_comb_generator::cdist_comb_generator(
  std::vector<PF_cdist*> &cdists, const arma::vec &start, const int nu,
  const arma::mat *xtra_covar, const double covar_fac, const double ftol_rel):
  cdists(cdists), nu(nu),
  /* just in case we do not set it */
  nlopt_result_code(NLOPT_RESULT_CODE_NOT_SET)
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
    nlopt_set_ftol_rel(opt, ftol_rel);
    nlopt_set_vector_storage(opt, 30L);

    double minf;
    nlopt_result_code = nlopt_optimize(opt, val.memptr(), &minf);
    nlopt_destroy(opt);

    if(nlopt_result_code < 1L)
      /* fall back to start value */
      val = start;

  } else
    nlopt_result_code = 1L;

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

nlopt_return_value_msg::nlopt_return_value_msg():
  nlopt_result_code(NLOPT_RESULT_CODE_NOT_SET), is_error(true), message() { }

nlopt_return_value_msg::nlopt_return_value_msg(const int nlopt_result_code):
  nlopt_result_code(nlopt_result_code),
  is_error(nlopt_result_code < 1L or nlopt_result_code > 4L),
  message(
    is_error ?
    "'nlopt' returned with codes" + std::to_string(nlopt_result_code) :
    "") { }


void nlopt_return_value_msgs::insert
  (const nlopt_return_value_msgs &other){
  for(auto o : other.msgs){
    if (msgs.find(o.first) != msgs.end())
      continue;

    msgs[o.first] = o.second;
    any_errors = any_errors or o.second.is_error;
  }
}

void nlopt_return_value_msgs::insert(const nlopt_return_value_msg &&val){
  if (msgs.find(val.nlopt_result_code) == msgs.end()){
    msgs[val.nlopt_result_code] = val;
    any_errors = any_errors or val.is_error;
  }
}

bool nlopt_return_value_msgs::has_any_errors() const {
  return any_errors;
}

std::string nlopt_return_value_msgs::message() const {
  if(!any_errors)
    return "";

  std::string out = "'nlopt' returned with codes";
  for(auto i : msgs)
    out += " " + std::to_string(i.second.nlopt_result_code);

  return out;
}

nlopt_return_value_msg cdist_comb_generator::get_result_code(){
  return nlopt_return_value_msg(nlopt_result_code);
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
