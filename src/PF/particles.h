#ifndef PARTICLES
#define PARTICLES

#include "../arma_n_rcpp.h"

struct particle {
  const arma::vec state;
  const particle *const parent;
  const arma::uword cloud_idx; /* zero based */
  const particle *const child;
  /* log-scale */
  double log_weight;
  double log_resampling_weight;

  particle(
    const arma::vec state, const particle *parent,
    const arma::uword cloud_idx, const particle *child = nullptr);
};

class cloud {
  using size_type = std::vector<particle>::size_type;
  using iterator = std::vector<particle>::iterator;
  using reverse_iterator = std::vector<particle>::reverse_iterator;

  std::vector<particle> particles;

public:
  particle New_particle(arma::vec state, const particle *parent, const particle *child);
  particle New_particle(arma::vec state, const particle *parent);
  iterator begin();
  iterator end();
  reverse_iterator rbegin();
  reverse_iterator rend();
  particle& operator[](size_type pos);
  void reserve (size_type n);
  size_type size();
};

#endif
