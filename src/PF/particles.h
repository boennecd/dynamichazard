#ifndef PARTICLES
#define PARTICLES

#include "../arma_n_rcpp.h"

class particle {
  arma::vec state;
  arma::uword cloud_idx; /* zero based */

public:
  const particle *parent;
  const particle *child;

  const arma::vec& get_state() const {
    return state;
  }

  arma::uword get_cloud_idx() const {
    return cloud_idx;
  }

  /* log-scale */
  double log_importance_dens;
  double log_weight;
  double log_likelihood_term;
  double log_resampling_weight;

  particle(
    const arma::vec state, const particle *parent,
    const arma::uword cloud_idx, const particle *child = nullptr);
  particle();
};

class cloud : private std::vector<particle> {
  using Base = std::vector<particle>;

public:
  using Base::size_type;
  using Base::iterator;
  using Base::const_iterator;
  using Base::reverse_iterator;

  cloud(): Base() {}
  cloud(size_type count): Base(count) {}

  void new_particle(
      arma::vec state, const particle *parent = nullptr, const particle *child = nullptr);
  particle& set_particle(
      arma::uword idx, arma::vec state,
      const particle *parent = nullptr, const particle *child = nullptr);
  arma::vec get_weigthed_mean();

  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::back;
  using Base::rend;
  using Base::operator[];
  using Base::reserve;
  using Base::size;
};

#endif
