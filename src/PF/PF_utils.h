#ifndef PF_UTILS
#define PF_UTILS

#define MAX(a,b) (((a)>(b))?(a):(b))

#include "particles.h"
#include "../BLAS_and_LAPACK/arma_utils.h"

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

/* ------------------------------------------- */

struct nothing {};

/* ------------------------------------------- */

template<typename iterator>
class work_block {
public:
  iterator start;
  iterator end;
  const unsigned long block_size;

  work_block(iterator start, iterator end, const unsigned long block_size):
    start(start), end(end), block_size(block_size) {}
};

template<typename iterator>
std::vector<work_block<iterator>> get_work_blocks(
    iterator begin, iterator end, const unsigned long block_size){
  unsigned long const length = std::distance(begin, end);
  unsigned long const num_blocks= (length + block_size - 1) / block_size;

  std::vector<work_block<iterator>> ans;
  ans.reserve(num_blocks);

  iterator block_start = begin;
  for(unsigned long i = 0; i < num_blocks - 1; ++i){
    iterator block_end = block_start;
    std::advance(block_end, block_size);

    ans.emplace_back(block_start, block_end, block_size);
    block_start = block_end;
  }

  ans.emplace_back(block_start, end, std::distance(block_start, end));

  return ans;
}

/* ------------------------------------------- */

struct input_for_normal_approximation {
  arma::vec mu;
  arma::mat Sigma_inv_chol;
  arma::mat Sigma_chol;
  arma::mat sigma_chol_inv;

  /* the conditional means for each of the particles given their parents */
  std::vector<arma::vec> mu_js;
};

template<typename densities>
static input_for_normal_approximation compute_mu_n_Sigma_from_normal_approximation(
    const PF_data &data, const unsigned int t, const PF_data::covarmat &Q, const arma::vec &alpha_bar){
  if(data.debug > 2){
    data.log(3) << "Computing normal approximation with mean vector:";
    data.log(3) << alpha_bar.t();
    data.log(3) << "and covaraince matrix:";
    data.log(3) << Q.mat;
  }

  input_for_normal_approximation ans;

  /* Compute the terms that does not depend on the outcome */
  /* Sigma^-1 = (Q + \tilde{Q})^{-1} */
  arma::uword p = alpha_bar.n_elem;
  arma::mat Sigma_inv = Q.inv;
  ans.mu = arma::vec(p, arma::fill::zeros);
  arma::vec &mu = ans.mu;

  /* Add the terms that does depend on the outcome */
  arma::uvec r_set = Rcpp::as<arma::uvec>(data.risk_sets[t - 1]) - 1;
  auto jobs = get_work_blocks(r_set.begin(), r_set.end(), data.work_block_size);
  unsigned int n_jobs = jobs.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(std::min(data.n_threads, (int)std::ceil(n_jobs / 2.)))
#endif
  for(unsigned int i = 0; i < n_jobs; ++i){
    auto &job = jobs[i];
    arma::uvec my_r_set(job.start, job.block_size, false /* don't copy */);

    arma::vec eta =  alpha_bar.t() * data.X.cols(my_r_set);
    const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

    auto it_eta = eta.begin();
    auto it_is_event = is_event.begin();
    auto it_r = my_r_set.begin();
    arma::uword n_elem = eta.n_elem;
    /*
      Update with:
        Signa = ... + X^T (-G) X
        mu = X^T (-G) X \bar{alpha} + X^T (-g)
    */
    arma::mat my_Sigma_inv(p, p, arma::fill::zeros);
    arma::vec my_mu(p, arma::fill::zeros);
    for(arma::uword i = 0; i < n_elem; ++i, ++it_eta, ++it_is_event, ++it_r){
      double g = densities::log_p_prime(*it_is_event, *it_eta, t);
      double neg_G = - densities::log_p_2prime(*it_is_event, *it_eta, t);

      sym_mat_rank_one_update(neg_G, data.X.col(*it_r), my_Sigma_inv);

      my_mu += data.X.col(*it_r) * ((*it_eta * neg_G) + g);
    }

#ifdef _OPENMP
#pragma omp critical
{
#endif
    Sigma_inv += my_Sigma_inv;
    mu += my_mu;
#ifdef _OPENMP
}
#endif
  }

  /* copy to lower */
  Sigma_inv = arma::symmatu(Sigma_inv);

  /* Compute needed factorizations */
  ans.Sigma_inv_chol = arma::chol(Sigma_inv);
  ans.Sigma_chol = arma::chol(arma::inv(Sigma_inv)); // TODO: do something smarter
  ans.sigma_chol_inv = arma::inv(arma::trimatu(ans.Sigma_chol));

  return ans;
}

template<typename densities>
static input_for_normal_approximation compute_mu_n_Sigma_from_normal_approximation(
    const PF_data &data, const unsigned int t, const PF_data::covarmat &Q, const arma::vec &alpha_bar,
    cloud &cl /* set mu_js when cloud is passed to */){
  auto ans = compute_mu_n_Sigma_from_normal_approximation<densities>(
    data, t, Q, alpha_bar);

  auto n_elem = cl.size();
  ans.mu_js = std::vector<arma::vec>(n_elem);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(std::min(data.n_threads, (int)std::ceil(n_elem / 10.)))
#endif
  for(unsigned int i = 0; i < n_elem; ++i){
    arma::vec mu_j = solve_w_precomputed_chol(Q.chol, cl[i].state) + ans.mu;
    mu_j = solve_w_precomputed_chol(ans.Sigma_inv_chol, mu_j);

    ans.mu_js[i] = std::move(mu_j);
  }

  return ans;
}


#undef MAX
#endif
