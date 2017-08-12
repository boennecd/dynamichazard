#ifndef PF_UTILS
#define PF_UTILS

#define MAX(a,b) (((a)>(b))?(a):(b))

#include "particles.h"

struct normalize_weights_output {
  double ESS = 0.;
  arma::vec weights;
};

template<double& (*F)(particle&), bool compute_ESS, bool update_particles>
inline normalize_weights_output normalize_weights(cloud &cl, const double max_weight){
  normalize_weights_output ans;
  double &ESS = ans.ESS;
  arma::vec &weights = ans.weights;
  weights.set_size(cl.size());

  auto w = weights.begin();
  double norm_constant = 0;
  for(auto it = cl.begin(); it != cl.end(); ++it, ++w){
    /* back transform weights */
    *w = exp(F(*it) - max_weight);
    norm_constant += *w;
  }

  w = weights.begin();
  for(auto it = cl.begin(); it != cl.end(); ++it, ++w){
    *w /= norm_constant;

    if(compute_ESS){
      ESS += *w * *w;
    }

    if(update_particles){
      /* Re-scale and take log */
      F(*it) = log(*w);
    }
  }

  if(compute_ESS){
    ESS = 1/ESS;
  }

  return ans;
}

template<double& (*F)(particle&), bool compute_ESS, bool update_particles>
inline normalize_weights_output normalize_weights(cloud &cl){
  double max_weight =  -std::numeric_limits<double>::max();
  for(auto it = cl.begin(); it != cl.end(); ++it){
    max_weight = MAX(F(*it), max_weight);
  }

  return(normalize_weights<F, compute_ESS, update_particles>(cl, max_weight));
}


inline double& normalize_log_weights_F(particle &p){
  return p.log_weight;
}
template<bool compute_ESS, bool update_particles>
inline normalize_weights_output normalize_log_weights(
    cloud &cl, const double max_weight){
  return normalize_weights<normalize_log_weights_F, compute_ESS, update_particles>(cl, max_weight);
}


inline double& normalize_log_resampling_weight_F(particle &p){
  return p.log_resampling_weight;
}
template<bool compute_ESS, bool update_particles>
inline normalize_weights_output normalize_log_resampling_weight(
    cloud &cl, const double max_weight){
  return normalize_weights<normalize_log_resampling_weight_F, compute_ESS, update_particles>(cl, max_weight);
}

struct nothing {};



#undef MAX
#endif
