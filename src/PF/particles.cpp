#include "particles.h"

particle::particle(
  const arma::vec state, const particle *parent,
  const int cloud_idx, const particle *child):
  state(state), parent(parent), cloud_idx(cloud_idx), child(child)
  {
    log_weight = log_resampling_weight = std::numeric_limits<double>::quiet_NaN();
  }


particle cloud::New_particle(
    arma::vec state, const particle *parent, const particle *child){
  particles.emplace_back(state, parent, particles.size(), child);
  return particles.back();
}

std::vector<particle>::iterator cloud::begin(){
  return particles.begin();
}

std::vector<particle>::iterator cloud::end(){
  return particles.end();
}

particle cloud::operator[](size_type pos){
  return particles[pos];
}

