#ifndef PF_UTILS
#define PF_UTILS

#include "particles.h"

inline void normalize_log_weights(
    cloud &cl, const double max_weight){
  double norm_constant = 0;
  for(auto it = cl.begin(); it != cl.end(); ++it){
    /* back transform weights */
    it->log_weight = exp(it->log_weight - max_weight);
    norm_constant += it->log_weight;
  }

  for(auto it = cl.begin(); it != cl.end(); ++it){
    /* Re-scale and take log */
    it->log_weight = log(it->log_weight / norm_constant);
  }
}

#endif
