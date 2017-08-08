#include "particles.h"

particle::particle(
  const arma::vec state, const particle *parent,
  const arma::uword cloud_idx, const particle *child):
  state(state), parent(parent), cloud_idx(cloud_idx), child(child)
  {
    log_weight = log_resampling_weight = std::numeric_limits<double>::quiet_NaN();
  }

particle cloud::New_particle(
    arma::vec state, const particle *parent, const particle *child){
  particles.emplace_back(state, parent, particles.size(), child);
  return particles.back();
}

particle cloud::New_particle(
    arma::vec state, const particle *parent){
  particles.emplace_back(state, parent, particles.size(), nullptr);
  return particles.back();
}

std::vector<particle>::iterator cloud::begin(){
  return particles.begin();
}

std::vector<particle>::iterator cloud::end(){
  return particles.end();
}

std::vector<particle>::reverse_iterator cloud::rbegin(){
  return particles.rbegin();
}

std::vector<particle>::reverse_iterator cloud::rend(){
  return particles.rend();
}

particle& cloud::operator[](size_type pos){
  return particles[pos];
}
void cloud::reserve (size_type n){
  particles.reserve(n);
}

std::vector<particle>::size_type cloud::size(){
  return particles.size();
}

