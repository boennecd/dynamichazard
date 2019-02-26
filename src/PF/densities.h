#ifndef DENSITIES
#define DENSITIES

#include "dists.h"

// TODO: check what is needed

#include "dmvnrm.h"
#include "PF_data.h"
#include "get_work_blocks.h"
#include "../family.h"
#include "particles.h"

class pf_dens {
  const std::string family;
  const PF_data &data;
  artificial_prior_generator art_gen;

public:
  pf_dens(const PF_data&, const std::string&);

  double log_prob_state_given_parent(const arma::vec&, const arma::vec&);
  double log_prob_state_given_child(const arma::vec&, const arma::vec&);
  std::unique_ptr<PF_cdist> get_fw_dist(const arma::vec&);
  std::unique_ptr<PF_cdist> get_bw_dist(const arma::vec&);
  std::shared_ptr<PF_cdist> get_prior(const arma::uword);
  std::shared_ptr<PF_cdist> get_y_dist(
      const int, const bool multithreaded = false);
};

std::string get_family(const std::string&);

#endif
