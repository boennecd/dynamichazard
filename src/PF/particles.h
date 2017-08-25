#ifndef PARTICLES
#define PARTICLES

#include "../arma_n_rcpp.h"

struct particle {
  const arma::vec state;
  const particle *const parent;
  const arma::uword cloud_idx; /* zero based */
  const particle *const child;
  /* log-scale */
  double log_importance_dens;
  double log_weight;
  double log_unnormalized_weight;
  double log_resampling_weight;

  particle(
    const arma::vec state, const particle *parent,
    const arma::uword cloud_idx, const particle *child = nullptr);
};

class cloud : private std::vector<particle> {
  using Base = std::vector<particle>;

public:
  using Base::size_type;
  using Base::iterator;
  using Base::const_iterator;
  using Base::reverse_iterator;

  particle& New_particle(arma::vec state, const particle *parent, const particle *child);
  particle& New_particle(arma::vec state, const particle *parent);
  arma::vec get_weigthed_mean();

  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::rend;
  using Base::operator[];
  using Base::reserve;
  using Base::size;
};

#endif
