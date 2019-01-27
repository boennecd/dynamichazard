#ifndef COND_APPROX_H
#define COND_APPROX_H
#include "dists.h"

/* This class will either yield an exact conditional distribution if
 * all the passed conditional distributions are multivariate Gaussian
 * and otherwise it will yield an approximation
 */
class cdist_comb {
  std::vector<PF_cdist*> &cdists;
  const arma::mat neg_K;
  const arma::mat Sig;
  const arma::vec k;

public:
  cdist_comb(std::vector<PF_cdist*>&, const arma::vec&);

  arma::vec sample(std::vector<arma::vec&>) const;
  double proposal_density(std::vector<arma::vec&>) const;
  arma::vec get_mean(std::vector<arma::vec&>) const;
  arma::mat get_covar(std::vector<arma::vec&>) const;
};

#endif
