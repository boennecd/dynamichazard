#include "particles.h"

particle::particle(
  const arma::vec state, const particle *parent,
  const arma::uword cloud_idx, const particle *child):
  state(state), cloud_idx(cloud_idx), parent(parent), child(child)
  {
    log_likelihood_term = log_importance_dens = log_weight =
      log_resampling_weight = std::numeric_limits<double>::quiet_NaN();
  }

particle::particle():
  particle::particle(arma::vec(), nullptr, 0, nullptr) {}

void cloud::new_particle(
    arma::vec state, const particle *parent, const particle *child){
  emplace_back(state, parent, size(), child);
}

void cloud::new_particle(particle &&other){
  emplace_back(other);
  back().cloud_idx = size() - 1L;
}

particle& cloud::set_particle(
    arma::uword idx, arma::vec state,
    const particle *parent, const particle *child){
  // use copy assigment
  particle &elem = this->at(idx);
  elem = particle(state, parent, idx, child);
  return elem;
}

/*
  Returns the weighted mean using cloud weights. Assumes that weights are
  scaled
*/
arma::vec cloud::get_weigthed_mean(){
  if(size() == 0)
    return arma::vec();

  arma::vec out(front().get_state().n_elem, arma::fill::zeros);
  for(auto it = begin(); it != end(); ++it){
    out += (it->get_state() * exp(it->log_weight));
  }

  return out;
}

