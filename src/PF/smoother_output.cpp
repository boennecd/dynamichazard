#include "PF_utils.h"

smoother_output::pair::pair(const particle *p, double log_weight):
  p(p), log_weight(log_weight) {}

smoother_output::particle_pairs::particle_pairs
  (const particle *p, const double log_weight, std::vector<pair> &&pairs):
  p(p), log_weight(log_weight), transition_pairs(pairs) { }

smoother_output::particle_pairs::particle_pairs
  (const particle *p, const double log_weight):
  particle_pairs(p, log_weight, std::vector<pair>()) {}

smoother_output::particle_pairs::particle_pairs():
  particle_pairs(nullptr, std::numeric_limits<double>::quiet_NaN()) { }

smoother_output::smoother_output():
  transition_likelihoods(std::make_shared<trans_like_obj>()) { }

std::shared_ptr<smoother_output::trans_like_obj>
  smoother_output::get_transition_likelihoods
  (const bool do_make_if_len_0) const
  {
    if(!do_make_if_len_0 || transition_likelihoods->size() > 0)
      return transition_likelihoods;

    const unsigned n_clouds = smoothed_clouds.size();
    std::shared_ptr<trans_like_obj> out =
      std::make_shared<trans_like_obj>(n_clouds);
    auto out_begin = out->begin();
    auto clo_begin = smoothed_clouds.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1000) if(n_clouds > 10000)
#endif
    for(unsigned int i = 0; i < n_clouds; ++i){
      auto clo = clo_begin + i;

      std::vector<particle_pairs> &new_elem = *(i + out_begin);
      new_elem.reserve(clo->size());
      for(cloud::const_iterator pa = clo->begin(); pa != clo->end(); ++pa){
        std::vector<pair> pairs(1);
        pairs[0].p = pa->parent;
        pairs[0].log_weight = 0;
        new_elem.emplace_back(&(*pa), pa->log_weight, std::move(pairs));
      }
    }

    return out;
  }
