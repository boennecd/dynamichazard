#ifndef COND_APPROX_H
#define COND_APPROX_H
#include "dists.h"
#include <memory>

class dist_comb {
public:
  virtual ~dist_comb() {}

  virtual arma::vec sample() const = 0;
  virtual double log_density(const arma::vec&) const = 0;
  virtual const arma::vec& get_mean() const = 0;
  virtual const arma::mat& get_covar() const = 0;
};

/* This class will either yield an exact conditional distribution if
 * all the passed conditional distributions are multivariate Gaussian
 * and otherwise it will yield an approximation
 */

struct nlopt_return_value_msg {
  nlopt_return_value_msg();
  nlopt_return_value_msg(const int);

  int nlopt_result_code;
  bool is_error;
  std::string message;
};

struct nlopt_return_value_msgs {
  std::map<int, nlopt_return_value_msg> msgs;
  bool any_errors = false;

public:
  void insert(const nlopt_return_value_msgs&);
  void insert(const nlopt_return_value_msg&&);
  bool has_any_errors() const;
  std::string message() const;
};

class cdist_comb_generator {
  std::vector<PF_cdist*> &cdists;
  arma::mat neg_K;
  std::shared_ptr<covarmat> Sig;
  arma::vec k;
  /* less than 2 implies multivariate t-distribution */
  const int nu;
  int nlopt_result_code;

public:
  cdist_comb_generator
    (std::vector<PF_cdist*>&, const arma::vec&,
     const int nu = -1L, const arma::mat *xtra_covar = nullptr);
  cdist_comb_generator
    (std::vector<PF_cdist*>&, const int nu = -1L,
     const arma::mat *xtra_covar = nullptr);

  nlopt_return_value_msg get_result_code();

  std::unique_ptr<dist_comb> get_dist_comb(
      const std::initializer_list<arma::vec*>&);
};

#endif
