#include "dynamichazard.h"

//logit_fam
static double logit_fam::link_func(const double &mu){
  return (log(mu / (1 - mu)));
}

double logit_fam::link_func_inv(const double &eta) const{
  double tmp = (eta < MTHRESH) ? DOUBLE_EPS :
  ((eta > THRESH) ? INVEPS : exp(eta));

  return( tmp / (1.0 + tmp));
};

double logit_fam::variance(const double &mu) const{
  return(mu * (1 - mu));
};

double logit_fam::d_mu_d_eta(const double& eta) const{
  double exp_eta = exp(eta);
  double opexp = 1 + exp_eta;

  return((eta > THRESH || eta < MTHRESH) ?  DOUBLE_EPS : exp_eta / (opexp * opexp));
}

double logit_fam::dev_resids(const double &y, const double &mu, const double &w) const{
  return (y > 0.0) ? - log(mu) : - log(1 - mu);
};



//poisson_fam
double poisson_fam::link_func(const double &mu) const{
  return log(mu);
}

double poisson_fam::link_func_inv(const double &eta) const{
  return std::max(exp(eta), DOUBLE_EPS);
}


double poisson_fam::variance(const double &mu) const{
  return mu;
}

double poisson_fam::d_mu_d_eta(const double &eta) const{
  return std::max(exp(eta), DOUBLE_EPS);
}

double poisson_fam::dev_resids(const double &y, const double &mu, const double &w) const{
  return (y > 0.0) ? w * (y * log(y / mu) - (y - mu)) : mu * w;
}
