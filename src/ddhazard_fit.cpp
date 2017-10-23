#include "ddhazard.h"
#include "family.h"
#include "estimate_fixed_effects_M_step.h"

using uword = arma::uword;
using bigglm_updateQR_logit   = bigglm_updateQR<logistic>;
using bigglm_updateQR_poisson = bigglm_updateQR<exponential>;

using SMA_logit =  SMA<SMA_hepler_logit>;
using SMA_exp =  SMA<SMA_hepler_exp>;

using GMA_logit =  GMA<GMA_hepler_logit>;
using GMA_exp = GMA<GMA_hepler_exp>;

// Define convergence criteria used later
inline double relative_norm_change(const arma::mat &prev_est, const arma::mat &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::mat&, const arma::mat&) = relative_norm_change;

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
                            const double denom_term = .0001,
                            const int n_fixed_terms_in_state_vec = 0,
                            const bool use_pinv = false,
                            const std::string criteria = "delta_coef",
                            const std::string posterior_version = "cholesky",
                            const signed int GMA_max_rep = 10,
                            const double GMA_NR_eps = 0.1,
                            const int EKF_batch_size = 5000L){
  if(Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) &&
     is_exponential_model(model)){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = true which should be false for model '" + model  +"'");
  } else if(!Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) && model == "logit"){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = false which should be true for model '" + model  +"'");
  }

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
  std::unique_ptr<ddhazard_data> p_data;
  std::unique_ptr<Solver> solver;

  if(method == "EKF"){
    p_data.reset(new ddhazard_data_EKF(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      NR_eps, LR,
      eps_fixed_parems, max_it_fixed_params, weights,
      n_max, eps, verbose,
      order_, est_Q_0, model != "logit", NR_it_max, debug, n_threads,
      denom_term, use_pinv, criteria, EKF_batch_size));
    if(model == "logit"){
      solver.reset(new EKF_solver<logistic>(
          static_cast<ddhazard_data_EKF &>(*p_data.get()), model));

    } else if (is_exponential_model(model)){
      solver.reset(new EKF_solver<exponential>(
          static_cast<ddhazard_data_EKF &>(*p_data.get()), model));

    } else
      Rcpp::stop("EKF is not implemented for model '" + model  +"'");

  } else if (method == "UKF"){
    if(model != "logit" &&
       !is_exponential_model(model))
      Rcpp::stop("UKF is not implemented for model '" + model  +"'");
    p_data.reset(new ddhazard_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params, weights,
      n_max, eps, verbose,
      order_, est_Q_0, debug, LR, n_threads, denom_term, use_pinv,
      criteria));

    if(model == "logit"){
      solver.reset(
        new UKF_solver_New<logistic>(*p_data.get(), kappa, alpha, beta));

    } else if (is_exponential_model(model)){
      solver.reset(
        new UKF_solver_New<exponential>(*p_data.get(), kappa, alpha, beta));

    } else
      Rcpp::stop("Model '", model ,"' is not implemented with UKF");

  } else if (method == "UKF_org"){
    if(model != "logit")
      Rcpp::stop("UKF_org is not implemented for model '" + model  +"'");

    p_data.reset(new ddhazard_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params,
      weights,
      n_max, eps, verbose,
      order_, est_Q_0, debug, LR, n_threads, denom_term, use_pinv,
      criteria));

    if(p_data->any_fixed_in_M_step)
      Rcpp::stop("Fixed effects is not implemented with UKF_org");

    solver.reset(new UKF_solver_Org(*p_data.get(), kappa));

  } else if (method == "SMA"){
    p_data.reset(new ddhazard_data(
        n_fixed_terms_in_state_vec,
        X, fixed_terms, tstart, tstop, is_event_in_bin,
        a_0, fixed_parems_start, Q_0, Q,
        risk_obj, F_,
        eps_fixed_parems, max_it_fixed_params,
        weights,
        n_max, eps, verbose,
        order_, est_Q_0, debug, LR, n_threads, denom_term, use_pinv,
        criteria));

    if(model == "logit"){
      solver.reset(new SMA_logit(*p_data.get(), posterior_version));
    } else if(is_exponential_model(model)){
      solver.reset(new SMA_exp(*p_data.get(), posterior_version));
    } else
      Rcpp::stop("Model '", model ,"' is not implemented with rank one posterior approximation");

  } else if (method == "GMA"){
    p_data.reset(new ddhazard_data(
        n_fixed_terms_in_state_vec,
        X, fixed_terms, tstart, tstop, is_event_in_bin,
        a_0, fixed_parems_start, Q_0, Q,
        risk_obj, F_,
        eps_fixed_parems, max_it_fixed_params,
        weights,
        n_max, eps, verbose,
        order_, est_Q_0, debug, LR, n_threads, denom_term, use_pinv,
        criteria));

    if(model == "logit"){
      solver.reset(new GMA_logit(*p_data.get(), GMA_max_rep, GMA_NR_eps));
    }  else if(is_exponential_model(model)){
      solver.reset(new GMA_exp(*p_data.get(), GMA_max_rep, GMA_NR_eps));
    }else
      Rcpp::stop("Model '", model ,"' is not implemented with rank one posterior approximation");

  } else{
    Rcpp::stop("method '" + method  + "'is not implemented");
  }

  arma::mat a_prev;
  double old_log_like = 0.0;
  if(p_data->criteria == "delta_coef"){
    a_prev.copy_size(p_data->a_t_t_s);
    a_prev.zeros();
  }

  arma::mat varying_only_F = p_data->F_;
  if(p_data->any_fixed_in_E_step){
    if(p_data->any_fixed_in_E_step){
      varying_only_F.shed_rows(p_data->span_fixed_params->a, p_data->span_fixed_params->b);
      varying_only_F.shed_cols(p_data->span_fixed_params->a, p_data->span_fixed_params->b);
    } else{
      varying_only_F = arma::mat();
    }
  }

  do
  {
    p_data->em_iteration = it;
    p_data->computation_stage = "Starting EM";

    if(p_data->debug){
      my_debug_logger(*p_data)
        << "##########################################";
      my_debug_logger(*p_data)
        << "Starting iteration " << it << " with the following values";
      my_print(*p_data, p_data->a_t_t_s.col(0), "a_0");
      my_print(*p_data, p_data->Q, "Q");
    }

    if((it + 1) % 25 == 0)
      Rcpp::checkUserInterrupt(); // this is expensive (on Windows) - you do not want to check too often


    if(p_data->any_dynamic){
      p_data->V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

      // E-step: filtering
      p_data->computation_stage = "E-step filtering";
      solver->solve();

      // E-step: smoothing
      p_data->computation_stage = "E-step smoothing";
      if(p_data->debug){
        my_debug_logger(*p_data) << "Started smoothing";
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
          my_print(*p_data, p_data->a_t_t_s.col(t), "a_(" + ss.str() + ")");
          my_print(*p_data, p_data->V_t_t_s.slice(t), "V_(" + ss.str() + ")");
          my_debug_logger(*p_data)
            << "Condition number of V_(" + ss.str() + ") is "
            << arma::cond(p_data->V_t_t_s.slice(t));
        }
      }

      // M-step
      p_data->computation_stage = "M-step";
      if(est_Q_0){
        Q_0 = p_data->V_t_t_s.slice(0);
      }

      if(p_data->debug){
        my_print(*p_data, p_data->Q, "Q before changes in M-step");
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


      if(p_data->any_fixed_in_E_step){
        Q.rows(*p_data->span_fixed_params).zeros();
        Q.cols(*p_data->span_fixed_params).zeros();
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
        arma::mat tmp_Q = Q(*p_data->span_current_cov, *p_data->span_current_cov);
        Q.zeros();
        Q(*p_data->span_current_cov, *p_data->span_current_cov) = tmp_Q;
      }
    }

    if(p_data->debug){
      my_print(*p_data, p_data->Q, "Q after changes in M-step");
    }

    if(p_data->criteria == "delta_coef"){
      if(p_data->any_fixed_in_E_step ||
         p_data->any_dynamic){
        conv_values.push_back(conv_criteria(a_prev(*p_data->span_current_cov, arma::span::all),
                                            p_data->a_t_t_s(*p_data->span_current_cov, arma::span::all)));
      } else
        conv_values.push_back(0.0);
    }

    if(p_data->any_fixed_in_M_step){
      arma::vec old = p_data->fixed_parems;

      if(model == "logit"){
        estimate_fixed_effects_M_step<bigglm_updateQR_logit>(
          p_data.get(), fixed_effect_chunk_size);

      } else if(is_exponential_model(model)){
        estimate_fixed_effects_M_step<bigglm_updateQR_poisson>(
          p_data.get(), fixed_effect_chunk_size);

      } else
        Rcpp::stop("Fixed effects is not implemented for '" + model  +"'");

      if(p_data->criteria == "delta_coef"){
        *(conv_values.end() -1) += conv_criteria(old, p_data->fixed_parems);
      }
    }

    double log_like = 0.0;
    if(p_data->criteria == "delta_likeli" || (verbose && it % 5 < verbose)){
      arma::mat varying_only_a = p_data->a_t_t_s; // take copy
      arma::vec fixed_effects_offsets;

      if(p_data->any_fixed_in_M_step){
        fixed_effects_offsets = p_data->fixed_terms.t() * p_data->fixed_parems;

      } else if(p_data->any_fixed_in_E_step){
        fixed_effects_offsets =
          p_data->X(*p_data->span_fixed_params, arma::span::all).t() *
          p_data->a_t_t_s(*p_data->span_fixed_params, arma::span::all).col(0);

        varying_only_a.shed_rows(p_data->span_fixed_params->a, p_data->span_fixed_params->b);

      } else{
        fixed_effects_offsets = arma::vec(p_data->X.n_cols, arma::fill::zeros);

      }

      log_like =
        logLike_cpp(p_data->X(*p_data->span_current_cov_varying, arma::span::all),
                    risk_obj,
                    varying_only_F,
                    Q_0(*p_data->span_current_cov_varying, *p_data->span_current_cov_varying),
                    Q(*p_data->span_current_cov_varying, *p_data->span_current_cov_varying),
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


      my_debug_logger(*p_data)
        << "Iteration " <<  std::setw(5)<< it + 1
        << " ended with conv criteria " << std::setw(15) << *(conv_values.end() -1)
        << "\t" << "The log likelihood of the mean path is " << log_like << std::setw(rcout_width);
    }

    if(*(conv_values.end() -1) < eps)
      break;

    if(p_data->criteria == "delta_coef"){
      a_prev = p_data->a_t_t_s;
    } else if(p_data->criteria == "delta_likeli"){
      old_log_like = log_like;
    }
  } while(++it < n_max);

  if(it == n_max)
    Rcpp::warning("EM algorithm did not converge within the n_max number of iterations");

  return(Rcpp::List::create(Rcpp::Named("V_t_d_s") = Rcpp::wrap(p_data->V_t_t_s),
                            Rcpp::Named("a_t_d_s") = Rcpp::wrap(p_data->a_t_t_s.t()),
                            Rcpp::Named("B_s") = Rcpp::wrap(p_data->B_s),
                            Rcpp::Named("lag_one_cov") = Rcpp::wrap(p_data->lag_one_cov),
                            Rcpp::Named("fixed_effects") = Rcpp::wrap(p_data->fixed_parems),

                            Rcpp::Named("n_iter") = std::min(it + 1, n_max),
                            Rcpp::Named("conv_values") = conv_values,
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("Q_0") = Rcpp::wrap(Q_0)));
}
