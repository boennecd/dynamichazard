#include "ddhazard.h"
#include "bigglm_wrapper.h"

using uword = arma::uword;
using bigglm_updateQR_logit   = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;

using Posterior_approx_logit =  Posterior_approx<Posterior_approx_hepler_logit>;
using Posterior_approx_exp =  Posterior_approx<Posterior_approx_hepler_exp>;

// Define convergence criteria used later
inline double relative_norm_change(const arma::mat &prev_est, const arma::mat &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::mat&, const arma::mat&) = relative_norm_change;


extern std::vector<double> logLike_cpp(const arma::mat&, const Rcpp::List&,
                                       const arma::mat&, const arma::mat&,
                                       arma::mat Q, const arma::mat&,
                                       const arma::vec&, const arma::vec&,
                                       const arma::vec &,
                                       const int, const std::string);



// Method to estimate fixed effects like in biglm::bigglm
template<typename T>
void estimate_fixed_effects(problem_data * const p_data, const int chunk_size,
                            bigglm_updateQR<T> &updater){
  int it_outer = 0;
  arma::vec old_beta;
  do{
    int cursor_risk_set = 0;
    int n_elements = 0;

    // Set up look variables
    int t = 1; // start looping at one to be consistent with other implementations
    arma::mat fixed_terms(p_data->fixed_parems.n_elem, chunk_size, arma::fill::none);
    arma::vec offsets(chunk_size, arma::fill::zeros);
    arma::vec y(chunk_size);
    arma::vec eta;
    arma::vec w(chunk_size);
    qr_obj qr(p_data->fixed_parems.n_elem);
    auto it = p_data->risk_sets.begin();
    double bin_stop = p_data->min_start;

    for(; it != p_data->risk_sets.end(); ++it, ++t){

      double bin_start = bin_stop;
      double delta_t = p_data->I_len[t - 1];
      bin_stop += delta_t;

      // Find the risk set and the number of elements to take
      arma::uvec r_set = Rcpp::as<arma::uvec>(*it) - 1;
      int r_set_size = r_set.n_elem;
      int n_elements_to_take = std::min(chunk_size - n_elements, r_set_size - cursor_risk_set);
      r_set = r_set.subvec(cursor_risk_set, cursor_risk_set + n_elements_to_take - 1);

      // Find the outcomes, fixed terms and compute the offsets
      y.subvec(n_elements, n_elements + n_elements_to_take - 1) =
        arma::conv_to<arma::vec>::from(p_data->is_event_in_bin.elem(r_set) == (t - 1));

      w.subvec(n_elements, n_elements + n_elements_to_take - 1) =
        p_data->weights(r_set);

      fixed_terms.cols(n_elements, n_elements + n_elements_to_take - 1) =
        p_data->fixed_terms.cols(r_set);

      if(p_data->any_dynamic){
        offsets.subvec(n_elements, n_elements + n_elements_to_take - 1) =
          p_data->X.cols(r_set).t() * p_data->a_t_t_s.col(t).head(p_data->n_params_state_vec);
      } else {
        offsets.subvec(n_elements, n_elements + n_elements_to_take - 1).fill(0.);
      }

      for(arma::uword i = 0; i < r_set.n_elem; ++i){
        offsets(n_elements + i) +=
          T().time_offset(std::min(p_data->tstop(r_set(i)), bin_stop)
                            - std::max(p_data->tstart(r_set(i)), bin_start));
      }

      n_elements += n_elements_to_take;

      if(n_elements == chunk_size){ // we have reached the chunk_size

        arma::vec eta = fixed_terms.t() * p_data->fixed_parems;
        updater.update(qr, fixed_terms, eta, offsets, y, w);

        n_elements = 0;
      } else if(it == --p_data->risk_sets.end()){ // there is no more bins to process

        y = y.subvec(0, n_elements - 1);
        w = w.subvec(0, n_elements - 1);
        fixed_terms = fixed_terms.cols(0, n_elements - 1);
        offsets = offsets.subvec(0, n_elements - 1);

        arma::vec eta =  fixed_terms.t() * p_data->fixed_parems;
        updater.update(qr, fixed_terms, eta, offsets, y, w);
      }

      if(cursor_risk_set + n_elements_to_take < r_set_size){ // there are still elements left in the bin
        cursor_risk_set = cursor_risk_set + n_elements_to_take;
        --it;
        --t;
        bin_stop -= delta_t;
      } else
        cursor_risk_set = 0;
    }

    old_beta = p_data->fixed_parems;
    p_data->fixed_parems = bigglm_regcf(qr);

  } while(++it_outer < p_data->max_it_fixed_params && // Key that this the first condition we check when we use &&
    arma::norm(p_data->fixed_parems - old_beta, 2) / (arma::norm(old_beta, 2) + 1e-8) > p_data->eps_fixed_parems);

  static bool failed_to_converge_once = false;
  if(it_outer == p_data->max_it_fixed_params && !failed_to_converge_once){
    failed_to_converge_once = true;
    std::stringstream msg;
    msg << "Failed to estimate fixed effects in " << p_data->max_it_fixed_params << " iterations at least once" << std::endl;
    Rcpp::warning(msg.str());
  }
}



// [[Rcpp::export]]
Rcpp::List ddhazard_fit_cpp(arma::mat &X, arma::mat &fixed_terms, // Key: assumed to have observations in the columns for performance due to column-major storage
                            const arma::vec &weights,
                            const arma::vec &tstart, const arma::vec &tstop,
                            const arma::colvec &a_0,
                            const arma::vec &fixed_parems_start,
                            arma::mat Q_0, // by value copy. This  is key cuz we will change it if est_Q_0 = T
                            arma::mat Q, // similary this is a copy
                            const Rcpp::List &risk_obj,
                            const arma::mat &F_,
                            const double eps_fixed_parems, const int max_it_fixed_params,
                            const arma::uword n_max = 100, const double eps = 0.001,
                            const arma::uword verbose = 0,
                            const int order_ = 1, const bool est_Q_0 = true,
                            const std::string method = "EKF",
                            Rcpp::Nullable<Rcpp::NumericVector> kappa = R_NilValue, // see this link for nullable example http://blogs.candoerz.com/question/164706/rcpp-function-for-adding-elements-of-a-vector.aspx
                            Rcpp::Nullable<Rcpp::NumericVector> alpha = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> beta = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> NR_eps = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> LR = R_NilValue,
                            const std::string model = "logit",
                            const std::string M_step_formulation = "Fahrmier94",
                            const int fixed_effect_chunk_size = 2e4,
                            const bool debug = false,
                            const unsigned int NR_it_max = 100,
                            const int n_threads = -1,
                            const double ridge_eps = .0001,
                            const int n_fixed_terms_in_state_vec = 0,
                            const bool use_pinv = false,
                            const std::string criteria = "delta_coef",
                            const std::string posterior_version = "cholesky"){
  if(Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) &&
     is_exponential_model(model)){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = true which should be false for model '" + model  +"'");
  } else if(!Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) && model == "logit"){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = false which should be true for model '" + model  +"'");
  }

#ifdef USE_OPEN_BLAS
  openblas_set_num_threads(std::max(n_threads, 1));
#endif

  // Declare non constants and intialize some of them
  double delta_t, test_max_diff;
  const double Q_warn_eps = sqrt(std::numeric_limits<double>::epsilon());

  Rcpp::NumericVector conv_values;

  uword it = 0;

  // M-stp pointers for convenience
  arma::mat *B, *V_less, *V;
  arma::vec a_less, a;

  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  // Intialize the solver for the E-step
  std::unique_ptr<problem_data> p_data;
  std::unique_ptr<Solver> solver;

  if(method == "EKF"){
    p_data.reset(new problem_data_EKF(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      NR_eps, LR,
      eps_fixed_parems, max_it_fixed_params, weights,
      n_max, eps, verbose,
      order_, est_Q_0, model != "logit", NR_it_max, debug, n_threads,
      ridge_eps, use_pinv, criteria));
    solver.reset(new EKF_solver(static_cast<problem_data_EKF &>(*p_data.get()), model));

  } else if (method == "UKF"){
    if(model != "logit" &&
       !is_exponential_model(model))
      Rcpp::stop("UKF is not implemented for model '" + model  +"'");
    p_data.reset(new problem_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params, weights,
      n_max, eps, verbose,
      order_, est_Q_0, debug, LR, n_threads, ridge_eps, use_pinv,
      criteria));

    if(model == "logit"){
      solver.reset(new UKF_solver_New_logit(*p_data.get(), kappa, alpha, beta));

    } else if (model == "exp_combined"){
      solver.reset(new UKF_solver_New_exponential(*p_data.get(), kappa, alpha, beta));

    } else if (model == "exp_bin"){
      solver.reset(new UKF_solver_New_exp_bin(*p_data.get(), kappa, alpha, beta));
    } else if (model == "exp_clip_time"){
      solver.reset(new UKF_solver_New_exp_clip_time(*p_data.get(), kappa, alpha, beta));
    } else if (model == "exp_clip_time_w_jump"){
      solver.reset(new UKF_solver_New_exp_clip_time_w_jump(*p_data.get(), kappa, alpha, beta));
    } else
      Rcpp::stop("Model '", model ,"' is not implemented with UKF");

  } else if (method == "UKF_org"){
    if(model != "logit")
      Rcpp::stop("UKF is not implemented for model '" + model  +"'");

    p_data.reset(new problem_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params,
      weights,
      n_max, eps, verbose,
      order_, est_Q_0, debug, LR, n_threads, ridge_eps, use_pinv,
      criteria));

    if(p_data->any_fixed_in_M_step)
      Rcpp::stop("Fixed effects is not implemented with UKF");

    solver.reset(new UKF_solver_Org(*p_data.get(), kappa));

  } else if (method == "post_approx"){
    p_data.reset(new problem_data(
        n_fixed_terms_in_state_vec,
        X, fixed_terms, tstart, tstop, is_event_in_bin,
        a_0, fixed_parems_start, Q_0, Q,
        risk_obj, F_,
        eps_fixed_parems, max_it_fixed_params,
        weights,
        n_max, eps, verbose,
        order_, est_Q_0, debug, LR, n_threads, ridge_eps, use_pinv,
        criteria));

    if(model == "logit"){
      solver.reset(new Posterior_approx_logit(*p_data.get(), posterior_version));
    } else if(is_exponential_model(model)){
      solver.reset(new Posterior_approx_exp(*p_data.get(), posterior_version));
    } else
      Rcpp::stop("Model '", model ,"' is not implemented with rank one posterior approximation");

  }else{
    Rcpp::stop("method '" + method  + "'is not implemented");
  }

  arma::mat a_prev;
  double old_log_like = 0.0;
  if(p_data->criteria == "delta_coef"){
    a_prev.copy_size(p_data->a_t_t_s);
    a_prev.zeros();
  }

  do
  {
    if(p_data->debug){
      if(it > 0)
        Rcpp::Rcout << "\n\n\n";
      Rcpp::Rcout << "##########################################\nStarting iteration " << it
                  << " with the following values" << std::endl;
      my_print(p_data->a_t_t_s.col(0), "a_0");
      my_print(p_data->Q.diag(), "diag(Q)");
    }


    if((it + 1) % 25 == 0)
      Rcpp::checkUserInterrupt(); // this is expensive (on Windows) - you do not want to check too often


    if(p_data->any_dynamic){
      p_data->V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

      // E-step: filtering
      solver->solve();

      // E-step: smoothing
      if(p_data->debug){
        Rcpp::Rcout << "Started smoothing" << std::endl;
      }

      for (int t = p_data->d - 1; t > -1; t--){
        // we need to compute the correlation matrix first
        if(t > 0){
          p_data->lag_one_cov.slice(t - 1) = p_data->V_t_t_s.slice(t) * p_data->B_s.slice(t - 1).t() +
            p_data->B_s.slice(t) * (
                p_data->lag_one_cov.slice(t) - F_ * p_data->V_t_t_s.slice(t)) * p_data->B_s.slice(t - 1).t();
        }

        p_data->a_t_t_s.col(t) = p_data->a_t_t_s.unsafe_col(t) + p_data->B_s.slice(t) *
          (p_data->a_t_t_s.unsafe_col(t + 1) - p_data->a_t_less_s.unsafe_col(t));
        p_data->V_t_t_s.slice(t) = p_data->V_t_t_s.slice(t) + p_data->B_s.slice(t) *
          (p_data->V_t_t_s.slice(t + 1) - p_data->V_t_less_s.slice(t)) * p_data->B_s.slice(t).t();

        if(p_data->debug){
          std::stringstream ss;
          ss << t << "|" <<  p_data->d;
          my_print(p_data->a_t_t_s.col(t), "a(" + ss.str() + ")");
          my_print(p_data->V_t_t_s.slice(t).diag(), "diag(V(" + ss.str() + "))");
        }
      }

      // M-step
      if(est_Q_0){
        Q_0 = p_data->V_t_t_s.slice(0);
      }
      Q.zeros();
      for (int t = 1; t < p_data->d + 1; t++){
        delta_t = p_data->I_len[t - 1];

        V_less = &p_data->V_t_t_s.slice(t - 1);
        V = &p_data->V_t_t_s.slice(t);
        a_less = p_data->a_t_t_s.unsafe_col(t - 1);
        a = p_data->a_t_t_s.unsafe_col(t);

        if(M_step_formulation == "Fahrmier94"){
          B = &p_data->B_s.slice(t - 1);

          Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
                  - F_ * *B * *V
                  - (F_ * *B * *V).t()
                  + F_ * *V_less * p_data->T_F_) / delta_t;

        } else if (M_step_formulation == "SmoothedCov"){
          B = &p_data->lag_one_cov.slice(t - 1); // this is not B but the lagged one smooth correlation. Do not mind the variable name

          Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
                  - F_ * *B
                  - (F_ * *B).t()
                  + F_ * *V_less * p_data->T_F_) / delta_t;
        } else
          Rcpp::stop("'M_step_formulation' of type '" + M_step_formulation + "' is not implemented");

      }
      Q /= p_data->d;

      if(p_data->debug){
        my_print(p_data->Q.diag(), "Diag(Q) before changes in M-step");
      }


      if(p_data->any_fixed_terms_in_state_vec){
        Q.rows(p_data->span_fixed_params).zeros();
        Q.cols(p_data->span_fixed_params).zeros();
      }

      if((test_max_diff = static_cast<arma::mat>(Q - Q.t()).max()) > Q_warn_eps){
        std::ostringstream warning;
        warning << "Q - Q.t() maximal element difference was " << test_max_diff <<
          " in iteration " << it + 1;
        Rcpp::warning(warning.str());
      }

      if((test_max_diff = static_cast<arma::mat>(Q_0 - Q_0.t()).max()) > Q_warn_eps){
        std::ostringstream warning;
        warning << "Q_0 - Q_0.t() maximal element difference was " << test_max_diff <<
          " in iteration " << it + 1;
        Rcpp::warning(warning.str());
      }

      // Ensure that Q and Q_0 are symmetric
      Q = (Q + Q.t()) / 2.0;
      Q_0 = (Q_0 + Q_0.t()) / 2.0;

      if(order_ > 1){
        arma::mat tmp_Q = Q(p_data->span_current_cov, p_data->span_current_cov);
        Q.zeros();
        Q(p_data->span_current_cov , p_data->span_current_cov) = tmp_Q;
      }
    }

    if(p_data->debug){
      my_print(p_data->Q.diag(), "Diag(Q) after changes in M-step");
    }

    if(p_data->criteria == "delta_coef"){
      if(p_data->any_fixed_terms_in_state_vec ||
         p_data->any_dynamic){
        conv_values.push_back(conv_criteria(a_prev(p_data->span_current_cov, arma::span::all),
                                            p_data->a_t_t_s(p_data->span_current_cov, arma::span::all)));
      } else
        conv_values.push_back(0.0);
    }

    if(p_data->any_fixed_in_M_step){
      arma::vec old = p_data->fixed_parems;

      if(model == "logit"){
        bigglm_updateQR_logit  updater;
        estimate_fixed_effects(p_data.get(), fixed_effect_chunk_size, updater);

      } else if(is_exponential_model(model)){
        bigglm_updateQR_poisson updater;
        estimate_fixed_effects(p_data.get(), fixed_effect_chunk_size, updater);

      } else
        Rcpp::stop("Fixed effects is not implemented for '" + model  +"'");

      if(p_data->criteria == "delta_coef"){
        *(conv_values.end() -1) += conv_criteria(old, p_data->fixed_parems);
      }
    }

    double log_like = 0.0;
    if(p_data->criteria == "delta_likeli" || (verbose && it % 5 < verbose)){
      arma::mat varying_only_F = p_data->F_; // take copy. TODO: only do this once
      arma::mat varying_only_a = p_data->a_t_t_s; // take copy
      arma::vec fixed_effects_offsets;

      if(p_data->any_fixed_in_M_step){
        fixed_effects_offsets = p_data->fixed_terms.t() * p_data->fixed_parems;

      } else if(p_data->any_fixed_terms_in_state_vec){
        fixed_effects_offsets =
          p_data->X(p_data->span_fixed_params, arma::span::all).t() *
          p_data->a_t_t_s(p_data->span_fixed_params, arma::span::all).col(0);

        varying_only_a.shed_rows(p_data->span_fixed_params.a, p_data->span_fixed_params.b);
        varying_only_F.shed_rows(p_data->span_fixed_params.a, p_data->span_fixed_params.b);
        varying_only_F.shed_cols(p_data->span_fixed_params.a, p_data->span_fixed_params.b);


      } else{
        fixed_effects_offsets = arma::vec(p_data->X.n_cols, arma::fill::zeros);

      }

      log_like =
        logLike_cpp(p_data->X(p_data->span_current_cov_varying, arma::span::all),
                    risk_obj,
                    varying_only_F,
                    Q_0(p_data->span_current_cov_varying, p_data->span_current_cov_varying),
                    Q(p_data->span_current_cov_varying, p_data->span_current_cov_varying),
                    varying_only_a,
                    p_data->tstart, p_data->tstop,
                    fixed_effects_offsets, order_, model)[0];

      if(p_data->criteria == "delta_likeli"){
        if(it == 0){
          conv_values.push_back(1e6); // something large
        } else{
          conv_values.push_back(std::abs((log_like - old_log_like) / (old_log_like - 1e-8)));
        }
      }
    }

    if(!p_data->any_dynamic) // No reason to take further iterations
      break;

    if(verbose && it % 5 < verbose){
      auto rcout_width = Rcpp::Rcout.width();


      Rcpp::Rcout << "Iteration " <<  std::setw(5)<< it + 1 <<
        " ended with conv criteria " << std::setw(15) << *(conv_values.end() -1) <<
          "\t" << "The log likelihood is " << log_like <<
            std::setw(rcout_width) << std::endl;
    }

    if(*(conv_values.end() -1) < eps)
      break;

    if(p_data->criteria == "delta_coef"){
      a_prev = p_data->a_t_t_s;
    } else if(p_data->criteria == "delta_likeli"){
      old_log_like = log_like;
    }
  }while(++it < n_max);

  if(it == n_max)
    Rcpp::warning("EM algorithm did not converge within the n_max number of iterations");

  return(Rcpp::List::create(Rcpp::Named("V_t_d_s") = Rcpp::wrap(p_data->V_t_t_s),
                            Rcpp::Named("a_t_d_s") = Rcpp::wrap(p_data->a_t_t_s.t()),
                            Rcpp::Named("B_s") = Rcpp::wrap(p_data->B_s),
                            Rcpp::Named("lag_one_cov") = Rcpp::wrap(p_data->lag_one_cov),
                            Rcpp::Named("fixed_effects") = Rcpp::wrap(p_data->fixed_parems),

                            Rcpp::Named("n_iter") = it + 1,
                            Rcpp::Named("conv_values") = conv_values,
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("Q_0") = Rcpp::wrap(Q_0)));
}


// Exported for test only
// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset, arma::vec &y,
                          const arma::vec &w){
  qr_obj qr;
  qr.D = std::shared_ptr<arma::vec>(&D, [](arma::vec*x) -> void { });
  qr.rbar = std::shared_ptr<arma::vec>(&rbar, [](arma::vec*x) -> void { });
  qr.thetab = std::shared_ptr<arma::vec>(&thetab, [](arma::vec*x) -> void { });
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = std::shared_ptr<arma::vec>(&tol, [](arma::vec*x) -> void { });

  if(model == "logit"){
    return(bigglm_updateQR_logit().update(qr, X, eta, offset, y, w));
  } else if (is_exponential_model(model)){
    return(bigglm_updateQR_poisson().update(qr, X, eta, offset, y, w));
  }
}
