#include "particles.h"

particle::particle(
  const arma::vec state, const particle *parent,
  const arma::uword cloud_idx, const particle *child):
  state(state), parent(parent), cloud_idx(cloud_idx), child(child)
  {
    log_unnormalized_weight = log_importance_dens = log_weight = log_resampling_weight =
      std::numeric_limits<double>::quiet_NaN();
  }

particle& cloud::New_particle(
    arma::vec state, const particle *parent, const particle *child){
  emplace_back(state, parent, size(), child);
  return back();
}

particle& cloud::New_particle(
    arma::vec state, const particle *parent){
  emplace_back(state, parent, size(), nullptr);
  return back();
}

/*
  Returns the weighted mean using cloud weights. Assumes that weights are
  scaled
*/
arma::vec cloud::get_weigthed_mean(){
  if(size() == 0)
    return arma::vec();

  arma::vec out(front().state.n_elem, arma::fill::zeros);
  for(auto it = begin(); it != end(); ++it){
    out += (it->state * exp(it->log_weight));
  }

  return out;
}

