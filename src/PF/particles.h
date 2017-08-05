#ifndef PARTICLES
#define PARTICLES

#include "../arma_n_rcpp.h"

struct particle {
  const arma::vec state;
  const particle *const parent;
  const int cloud_idx;
  const particle *const child;
  /* log-scale */
  double log_weight;
  double log_resampling_weight;

  particle(
    const arma::vec state, const particle *parent,
    const int cloud_idx, const particle *child = nullptr);
};

class cloud {
  using size_type = std::vector<particle>::size_type;
  using iterator = std::vector<particle>::iterator;

  std::vector<particle> particles;

public:
  particle New_particle(
      arma::vec state, const particle *parent, const particle *child = nullptr);
  iterator begin();
  iterator end();
  particle operator[](size_type pos);
};

#endif
