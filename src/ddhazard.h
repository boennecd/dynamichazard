// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <algorithm>
#include "ddhazard_problem_data.h"
#include "arma_n_rcpp.h"

#ifndef DDHAZARD_MISC
#define DDHAZARD_MISC

inline bool is_exponential_model(std::string model){
  return(model == "exp_combined" ||
         model == "exp_bin" ||
         model == "exp_clip_time" ||
         model == "exp_clip_time_w_jump");
}
#endif

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;
};

// Classes for EKF method
class EKF_filter_worker{
protected:
  virtual void do_comps(const arma::uvec::const_iterator it, int &i,
                        const arma::vec &i_a_t, const bool &compute_z_and_H,
                        const int &bin_number,
                        const double &bin_tstart, const double &bin_tstop) = 0; // abstact method to be implemented

  bool is_first_call;
  problem_data_EKF &dat;

  // local variables to compute temporary result
  arma::colvec u_;
  arma::mat U_;

public:
  EKF_filter_worker(problem_data_EKF &p_data);

  void operator()(arma::uvec::const_iterator first, const arma::uvec::const_iterator &last,
                const arma::vec &i_a_t, const bool &compute_z_and_H,
                const int &i_start, const int &bin_number,
                const double &bin_tstart, const double &bin_tstop);
};



class EKF_helper{
  unsigned long const max_threads;
  problem_data_EKF &p_data;
  std::vector<std::shared_ptr<EKF_filter_worker> > workers;
  const std::string model;

public:
  EKF_helper(problem_data_EKF &p_data_, const std::string model_);

  void parallel_filter_step(arma::uvec::const_iterator first, arma::uvec::const_iterator last,
                            const arma::vec &i_a_t,
                            const bool &compute_H_and_z,
                            const int &bin_number,
                            const double &bin_tstart, const double &bin_tstop);
};


class EKF_solver : public Solver{
  problem_data_EKF &p_dat;
  EKF_helper filter_helper;


public:
  EKF_solver(problem_data_EKF &p_, const std::string model);

  void solve();
};
