#ifndef COND_APPROX_H
#define COND_APPROX_H
#include "dists.h"
#include <memory>

class dist_comb {
public:
  virtual ~dist_comb() = default;

  virtual arma::vec sample() const = 0;
  virtual double log_density(const arma::vec&) const = 0;
  virtual const arma::vec& get_mean() const = 0;
  virtual const arma::mat& get_covar() const = 0;
};

/* This class will either yield an exact conditional distribution if
 * all the passed conditional distributions are multivariate Gaussian
 * and otherwise it will yield an approximation
 */

class cdist_comb_generator {
  std::vector<PF_cdist*> &cdists;
  arma::mat neg_K;
  std::shared_ptr<covarmat> Sig;
  arma::vec k;
  /* less than 2 implies multivariate t-distribution */
  const int nu;

public:
  cdist_comb_generator
    (std::vector<PF_cdist*>&, const arma::vec&,
     const int nu = -1L, const arma::mat *xtra_covar = nullptr,
     const double covar_fac = -1, const double ftol_rel = 1e-6);
  cdist_comb_generator
    (std::vector<PF_cdist*>&, const int nu = -1L,
     const arma::mat *xtra_covar = nullptr, const double covar_fac = -1,
     const double ftol_rel = 1e-6);

  std::unique_ptr<dist_comb> get_dist_comb(
      const std::initializer_list<arma::vec*>&);
};

#endif
