#include "particles.h"

particle::particle(particle *parent, const int cloud_idx, particle *child):
  parent(parent), cloud_idx(cloud_idx), child(child) {
  log_weight = log_resampling_weight = std::numeric_limits<double>::quiet_NaN();
}

particle cloud::New_particle(particle *parent, particle *child){
  particles.emplace_back(parent, particles.size(), child);
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

