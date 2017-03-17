#include "dynamichazard.h"
#include "exp_model_funcs.h"

using uword = arma::uword;

bool is_exponential_model(std::string model){
  return(model == "exp_combined" ||
         model == "exp_bin" ||
         model == "exp_clip_time" ||
         model == "exp_clip_time_w_jump");
}

// Define convergence criteria
inline double relative_norm_change(const arma::mat &prev_est, const arma::mat &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::mat&, const arma::mat&) = relative_norm_change;




// This is the orginal UKF formulation from:
// Julier, Simon J., and Jeffrey K. Uhlmann. "New extension of the Kalman filter
// to nonlinear systems." AeroSense'97. International Society for Optics and
// Photonics, 1997

// Sigma points are though re-generated as suggested in:
// Menegaz, Henrique Marra Taira. "Unscented kalman filtering on euclidean and
// riemannian manifolds." (2016).

// The methods requires inversion of matrix with dimension equal to the dim of
// observational equation. Hence, it does not scale well in the number of
// observation per bin. The code is kept to test against
class UKF_solver_Org : public Solver{
  problem_data &p_dat;
  const uword m;
  const double k;
  const double w_0;
  const double w_i;
  const double sqrt_m_k;
  arma::mat sigma_points;

  inline void compute_sigma_points(const arma::vec &a_t,
                                   arma::mat &s_points,
                                   const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_Org(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_):
  p_dat(p_),
  m(p_.a_t_t_s.n_rows),
  k(!k_.isNull() ?
      Rcpp::as< Rcpp::NumericVector >(k_)[0] :
      m * (1 + 1 * (.1 - 1)) / 1 * (1 - .1)),

  w_0(k / (m + k)),
  w_i(1 / (2 * (m + k))),
  sqrt_m_k(std::sqrt(m + k)),
  sigma_points(arma::mat(m, 2 * m + 1))
  {}

  void solve(){
#ifdef USE_OPEN_BLAS
    const int prev_n_thread = openblas_get_num_threads();
    openblas_set_num_threads(p_dat.n_threads);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

    double event_time = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){

      double delta_t = p_dat.I_len[t - 1];
      event_time += delta_t;

      // Update sigma pooints
      compute_sigma_points(p_dat.a_t_t_s.unsafe_col(t - 1),
                           sigma_points, p_dat.V_t_t_s.slice(t - 1));

      // E-step: Filter step
      // First we compute the mean
      p_dat.a_t_less_s.col(t - 1) = w_0 * sigma_points.unsafe_col(0) +
        w_i * arma::sum(sigma_points.cols(1, sigma_points.n_cols - 1), 1);

      // Then the variance
      p_dat.V_t_less_s.slice(t - 1) = delta_t * p_dat.Q;
      for(uword i = 0; i < sigma_points.n_cols; ++i){
        const double &w = i == 0 ? w_0 : w_i;

        p_dat.V_t_less_s.slice(t - 1) +=
          (w * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
          (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1)).t();
      }

      // Re-generate
      compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                           sigma_points, p_dat.V_t_less_s.slice(t - 1));

      // E-step: correction-step
      // First, we compute the mean of the outcome
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

      arma::mat Z_t = (sigma_points.t() * p_dat.X.cols(r_set)).t(); // we transpose due to the column-major
      Z_t = arma::trunc_exp(Z_t);
      Z_t.for_each([](arma::mat::elem_type &val) { val = val / (1 + val); });

      // Compute y_bar, P_a_v and P_v_v
      const arma::vec y_bar = w_0 * Z_t.unsafe_col(0) +
        w_i * arma::sum(Z_t.cols(1, Z_t.n_cols - 1), 1);

      arma::mat P_a_v = (w_0 * (sigma_points.unsafe_col(0) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
        (Z_t.unsafe_col(0) - y_bar).t();

      arma::mat P_v_v = (w_0 * (Z_t.unsafe_col(0) - y_bar)) * (Z_t.unsafe_col(0) - y_bar).t() +
        arma::diagmat(w_0 * Z_t.unsafe_col(0) % (1 - Z_t.unsafe_col(0)));

      for(uword i = 1; i < sigma_points.n_cols; ++i){
        P_a_v += (w_i * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
          (Z_t.unsafe_col(i) - y_bar).t();
        P_v_v += (w_i * (Z_t.unsafe_col(i) - y_bar)) * (Z_t.unsafe_col(i) - y_bar).t() +
          arma::diagmat(w_i * Z_t.unsafe_col(i) % (1 - Z_t.unsafe_col(i)));
      }

      P_v_v.diag() += p_dat.ridge_eps;

      // Compute new estimates
      inv_sympd(P_v_v, P_v_v, p_dat.use_pinv, "ddhazard_fit_cpp estimation error: Failed to invert P_v_v");

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) +
        P_a_v * (P_v_v * ((p_dat.is_event_in_bin(r_set) == t - 1) - y_bar));

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) - P_a_v * P_v_v * P_a_v.t();

      p_dat.B_s.slice(t - 1) = arma::solve(p_dat.V_t_less_s.slice(t - 1), p_dat.F_ * p_dat.V_t_t_s.slice(t - 1)).t();
    }

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(prev_n_thread);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif
  }
};


// New method that use the Woodbury matrix identity to make the UKF algorithm
// scale linearly with the dimension of the observationally.

// beta = 0.0 and alpha = 1.0 yields the same sigma points as in:
// Julier, Simon J., and Jeffrey K. Uhlmann. "New extension of the Kalman filter
// to nonlinear systems." AeroSense'97. International Society for Optics and
// Photonics, 1997

// Altering these will yields parameter estimate similar to:
// Wan, Eric A., and Rudolph Van Der Merwe. "The unscented Kalman filter for
// nonlinear estimation." Adaptive Systems for Signal Processing, Communications
// and Control Symposium 2000. AS-SPCC. The IEEE 2000. Ieee, 2000.

// We re-generate the sigma points as suggested in:
// Menegaz, Henrique Marra Taira. "Unscented kalman filtering on euclidean and
// riemannian manifolds." (2016).
class UKF_solver_New : public Solver{
protected:
  problem_data &p_dat;
  const uword m;
  const double a;
  const double k;
  const double b;
  const double lambda;
  const double w_0;
  const double w_0_c;
  const double w_0_cc;
  const double w_i;
  const double sqrt_m_lambda;
  arma::mat sigma_points;

  arma::vec weights_vec;
  arma::vec weights_vec_inv;
  arma::vec weights_vec_c;
  arma::vec weights_vec_c_inv;
  arma::vec weights_vec_cc;

  virtual void Compute_intermediates(const arma::uvec &r_set,
                                     const arma::vec offsets,
                                     const int t,
                                     const double bin_tstart, const double bin_tstop,
                                     arma::vec &c_vec, arma::mat &O) = 0;

  void compute_sigma_points(const arma::vec &a_t,
                            arma::mat &s_points,
                            const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_New(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                 Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                 Rcpp::Nullable<Rcpp::NumericVector> &beta):
  p_dat(p_),
  m(p_.a_t_t_s.n_rows),

  a(!alpha.isNull() ? Rcpp::as< Rcpp::NumericVector >(alpha)[0] : 1e-1),
  //k(!kappa.isNull() ? Rcpp::as< Rcpp::NumericVector >(kappa)[0] : m / pow(a, 2.0) - m),
  k(!kappa.isNull() ?
      Rcpp::as< Rcpp::NumericVector >(kappa)[0] :
      m * (1 + pow(a, 2) * (.1 - 1)) / (pow(a, 2) * (1 - .1))),
  b(!beta.isNull() ? Rcpp::as< Rcpp::NumericVector >(beta)[0] : 2.0),
  lambda(pow(a, 2) * (m + k) - m),

  w_0(lambda / (m + lambda)),
  w_0_c(w_0 + 1 - pow(a, 2) + b),
  w_0_cc(w_0 + 1 - a),
  w_i(1 / (2 * (m + lambda))),
  sqrt_m_lambda(std::sqrt(m + lambda)),

  sigma_points(arma::mat(m, 2 * m + 1))
  {
    if(w_0 == 0)
      Rcpp::stop("UKF not implemented for hyperparameters that yield zero weight on first sigma point");

    weights_vec = arma::vec(std::vector<double>(2 * m + 1, w_i));
    weights_vec[0] = w_0;
    weights_vec_inv = arma::pow(weights_vec, -1);

    weights_vec_c = weights_vec;
    weights_vec_c[0] = w_0_c;
    weights_vec_c_inv = arma::pow(weights_vec_c, -1);

    weights_vec_cc = weights_vec;
    weights_vec_cc[0] = w_0_cc;

    if(p_dat.debug){
      Rcpp::Rcout << "alpha, beta, kappa, n_dim = "
                  << a << ", "
                  << b << ", "
                  << k << ", "
                  << m << std::endl;
      Rcpp::Rcout << "w_0, w_0_c, w_0_cc, w_i, lambda = "
                  << w_0 << ", "
                  << w_0_c << ", "
                  << w_0_cc << ", "
                  << w_i << ", "
                  << lambda << std::endl;
    }
  }

  void solve(){
#ifdef USE_OPEN_BLAS
    const int prev_n_thread = openblas_get_num_threads();
    openblas_set_num_threads(p_dat.n_threads);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

    const arma::vec offsets = p_dat.any_fixed_in_M_step ?
      p_dat.fixed_terms.t() * p_dat.fixed_parems : arma::vec(p_dat.X.n_cols, arma::fill::zeros);
    double bin_stop = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){
      double bin_start = bin_stop;
      double delta_t = p_dat.I_len[t - 1];
      bin_stop += delta_t;

      // E-step: Filter step
      p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
      p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

      // Re-generate
      if(p_dat.debug){
        my_print(p_dat.V_t_less_s.slice(t - 1), "Chol decomposing for regenerations:");
      }


      compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                           sigma_points, p_dat.V_t_less_s.slice(t - 1));

      if(p_dat.debug){
        my_print(sigma_points, "new sigma points");
      }

      // E-step: correction-step
      arma::mat O;
      arma::vec c_vec;
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

      Compute_intermediates(r_set, offsets(r_set), t, bin_start, bin_stop, c_vec, O);

      // Substract mean to get delta sigma points
      arma::mat delta_sigma_points = sigma_points.each_col() - p_dat.a_t_less_s.unsafe_col(t - 1);

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) + p_dat.LR * delta_sigma_points * (weights_vec_cc % c_vec);

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) -
        (delta_sigma_points.each_row() % weights_vec_cc.t()) * O * (delta_sigma_points.each_row() % weights_vec_cc.t()).t();

      if(p_dat.debug){
        std::stringstream str, str_less;
        str_less << t << "|" << t - 1;
        str << t << "|" << t;

        my_print(p_dat.V_t_less_s.slice(t - 1).diag(), "diag(V_(" + str_less.str() + "))");
        my_print(p_dat.V_t_t_s.slice(t).diag(), "diag(V_(" + str.str()  + "))");
        my_print(p_dat.a_t_less_s.col(t - 1), "a_(" + str_less.str() + ")");
        my_print(p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
      }

      // Solve should be faster than inv_sympd http://arma.sourceforge.net/docs.html#inv_sympd
      // Solves yields a solution X to A * X = B <=> X = A^-1 B
      // We are looking at:
      //  X = B A^-1
      // X^T = A^-1 B^T <=> A X^T = B^T
      p_dat.B_s.slice(t - 1) = arma::solve(
        p_dat.V_t_less_s.slice(t - 1), p_dat.F_ * p_dat.V_t_t_s.slice(t - 1)).t();
    }

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(prev_n_thread);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif
  }
};


class UKF_solver_New_logit : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1: Compute expected outcomes given sigma points **
    O = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // we transpose due to the column-major

    O.each_col() += offsets;

    O = arma::trunc_exp(O);

    O.for_each([](arma::mat::elem_type &val) { val = val / (1 + val); });

    // ** 2: Compute mean observation sing sigma points **
    const arma::vec y_bar = w_0 * O.unsafe_col(0) +
      w_i * arma::sum(O.cols(1, O.n_cols - 1), 1);

    // ** 3: Compute variances and expected variance **
    arma::vec vars = w_0_c * (O.unsafe_col(0) % (1.0 - O.unsafe_col(0)));
    for(uword i = 1; i < O.n_cols; ++i){
      vars += w_i * (O.unsafe_col(i) % (1.0 - O.unsafe_col(i)));
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() % (p_dat.weights(r_set) / vars)).t() * ((p_dat.is_event_in_bin(r_set) == t - 1) - y_bar);
    O = O.t() * (O.each_col() % (p_dat.weights(r_set) / vars));

    // Compute intermediate matrix
    arma::mat tmp_mat;
    inv(tmp_mat, arma::diagmat(weights_vec_inv) + O, p_dat.use_pinv,
       "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    //arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O);
    inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O, p_dat.use_pinv,
       "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_logit(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                       Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                       Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exp_bin : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      at_risk_length(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_chance_die(v, inv_exp_v);
        vars(j) += w_c * exp_model_funcs::var_chance_die(v, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() % (p_dat.weights(r_set) / vars)).t() * (do_die - y_bar);
    O = O.t() * (O.each_col() % (p_dat.weights(r_set) / vars));

    // Compute intermediate matrix
    arma::mat tmp_mat;
    inv(tmp_mat, arma::diagmat(weights_vec_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_bin(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};

class UKF_solver_New_exp_clip_time : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
      bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_time(v, at_risk_length(j), inv_exp_v, exp_eta);
        vars(j) += w_c * exp_model_funcs::var_wait_time(v, at_risk_length(j),  exp_eta, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() % (p_dat.weights(r_set) / vars)).t() * (time_outcome - y_bar);
    O = O.t() * (O.each_col() % (p_dat.weights(r_set) / vars));

    // Compute intermediate matrix
    arma::mat tmp_mat;
    inv(tmp_mat, arma::diagmat(weights_vec_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_clip_time(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exp_clip_time_w_jump : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
        bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);

      if(do_die(i)){
        // we deduct the at risk lenght if the indvidual dies
        time_outcome(i) -= at_risk_length(i);
      }
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double inv_exp_eta = pow(exp_eta, -1);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(- v);

        O(j, i) = exp_model_funcs::expect_time_w_jump(exp_eta, inv_exp_eta, inv_exp_v, at_risk_length(j));
        vars(j) += w_c * exp_model_funcs::var_wait_time_w_jump(exp_eta, inv_exp_v, at_risk_length(j));
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() % (p_dat.weights(r_set) / vars)).t() * (time_outcome - y_bar);
    O = O.t() * (O.each_col() % (p_dat.weights(r_set) / vars));

    // Compute intermediate matrix
    arma::mat tmp_mat;
    inv(tmp_mat, arma::diagmat(weights_vec_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_clip_time_w_jump(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exponential : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O)
  {
    // See comments in UKF_solver_New_logit. The main difference here is that
    // we have tuples as outcomes. Thus, we have to deal with covariance terms

    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk * 2, sigma_points.n_cols);
    arma::vec vars(n_risk * 2, arma::fill::zeros);
    arma::vec covars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk * 2, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk and the outcome in case the
    // individual dies before the end of the bin
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
      bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_chance_die(v, inv_exp_v);
        vars(j) += w_c * exp_model_funcs::var_chance_die(v, inv_exp_v);

        covars(j) += w_c * exp_model_funcs::covar(v, inv_exp_v, exp_eta);

        O(j + n_risk, i) = exp_model_funcs::expect_time(
          v, at_risk_length(j), inv_exp_v, exp_eta);
        vars(j + n_risk) += w_c * exp_model_funcs::var_wait_time(
          v, at_risk_length(j), exp_eta, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Defines span to avoid 40-error
    arma::span span_binary(0, n_risk - 1);
    arma::span span_time(n_risk, 2 * n_risk -1);

    // Compute diagonal terms of inverse covariance matrix
    arma::vec inv_covmat_diag(2 * n_risk);
    inv_covmat_diag.subvec(span_binary) =
      (vars.subvec(span_binary) - covars % arma::pow(vars.subvec(span_time), -1) % covars);

    inv_covmat_diag.subvec(span_time) =
      (vars.subvec(span_time) - covars % arma::pow(vars.subvec(span_binary), -1) % covars);

    inv_covmat_diag.transform([](double val) { return(1 / val); });

    // Compute off diagonal terms of inverse covariance matrix
    arma::vec inv_covmat_off_diag(n_risk);
    inv_covmat_off_diag = -1 * inv_covmat_diag.subvec(span_time) % covars / vars.subvec(span_binary);

    // Compute inverse covariance matrix dot centered sigma points
    O.each_col() -= y_bar; // center

    O = O.t(); // transpose to ease the next computations

    arma::mat tmp_mat = O.each_row() % inv_covmat_diag.t();
    tmp_mat.cols(span_binary) +=
      O.cols(span_time).each_row() %  inv_covmat_off_diag.t();
    tmp_mat.cols(span_time) +=
      O.cols(span_binary).each_row() %  inv_covmat_off_diag.t();

    tmp_mat.cols(span_time).each_row() %= p_dat.weights(r_set).t();
    tmp_mat.cols(span_binary).each_row() %= p_dat.weights(r_set).t();

    O = O.t(); // transpose back

    {
      arma::vec outcome(n_risk * 2);
      outcome.subvec(span_binary) = arma::conv_to<arma::vec>::from(do_die);
      outcome.subvec(span_time) = time_outcome;

      c_vec = tmp_mat * (outcome - y_bar);
    }

    O = tmp_mat * O;
    inv(tmp_mat, arma::diagmat(weights_vec_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O, p_dat.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    tmp_mat = O * tmp_mat;

    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exponential(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                             Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                             Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {
  }
};















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
                            const std::string criteria = "delta_coef"){
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

      // E-step
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
