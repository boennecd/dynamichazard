#ifndef PARTICLES
#define PARTICLES

#include "../arma_n_rcpp.h"

struct particle {
  particle *const parent;
  const int cloud_idx;
  particle *const child;
  arma::vec state;
  /* log-scale */
  double log_weight;
  double log_resampling_weight;

  particle(particle *parent, const int cloud_idx, particle *child = nullptr);
};

class cloud {
  using size_type = std::vector<particle>::size_type;
  using iterator = std::vector<particle>::iterator;

  std::vector<particle> particles;

public:
  particle New_particle(particle *parent, particle *child = nullptr);
  iterator begin();
  iterator end();
  particle operator[](size_type pos);
};

#endif
