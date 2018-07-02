#include "../ddhazard.h"
#include "../utils.h"
#include "../family.h"

// This is the orginal UKF formulation from:
// Julier, Simon J., and Jeffrey K. Uhlmann. "New extension of the Kalman filter
// to nonlinear systems." AeroSense'97. International Society for Optics and
// Photonics, 1997

// Sigma points are re-generated as suggested in:
// Menegaz, Henrique Marra Taira. "Unscented kalman filtering on euclidean and
// riemannian manifolds." (2016).

// The methods requires inversion of matrix with dimension equal to the dim of
// observational equation. Hence, it does not scale well in the number of
// observation per bin. The code is kept to test against

UKF_solver_Org::UKF_solver_Org(ddhazard_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_):
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

void UKF_solver_Org::solve(){
    double event_time = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){

      double delta_t = p_dat.I_len[t - 1];
      event_time += delta_t;

      // Update sigma pooints
      compute_sigma_points(p_dat.a_t_t_s.col(t - 1),
                           sigma_points, p_dat.V_t_t_s.slice(t - 1));

      // E-step: Filter step
      // First we compute the mean
      p_dat.a_t_less_s.col(t - 1) = w_0 * sigma_points.col(0) +
        w_i * arma::sum(sigma_points.cols(1, sigma_points.n_cols - 1), 1);

      // Then the variance
      p_dat.V_t_less_s.slice(t - 1) = delta_t * p_dat.Q;
      for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
        const double &w = i == 0 ? w_0 : w_i;

        p_dat.V_t_less_s.slice(t - 1) +=
          (w * (sigma_points.col(i) - p_dat.a_t_less_s.col(t - 1))) *
          (sigma_points.col(i) - p_dat.a_t_less_s.col(t - 1)).t();
      }

      // Re-generate
      compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                           sigma_points, p_dat.V_t_less_s.slice(t - 1));

      // E-step: correction-step
      // First, we compute the mean of the outcome
      arma::uvec r_set = get_risk_set(p_dat, t);

      arma::mat Z_t = (sigma_points.t() * p_dat.X.cols(r_set)).t(); // we transpose due to the column-major
      Z_t = arma::trunc_exp(Z_t);
      Z_t.for_each([](arma::mat::elem_type &val) { val = val / (1 + val); });

      // Compute y_bar, P_a_v and P_v_v
      const arma::vec y_bar = w_0 * Z_t.col(0) +
        w_i * arma::sum(Z_t.cols(1, Z_t.n_cols - 1), 1);

      arma::mat P_a_v = (w_0 * (sigma_points.col(0) - p_dat.a_t_less_s.col(t - 1))) *
        (Z_t.col(0) - y_bar).t();

      arma::mat P_v_v = (w_0 * (Z_t.col(0) - y_bar)) * (Z_t.col(0) - y_bar).t() +
        arma::diagmat(w_0 * Z_t.col(0) % (1 - Z_t.col(0)));

      for(arma::uword i = 1; i < sigma_points.n_cols; ++i){
        P_a_v += (w_i * (sigma_points.col(i) - p_dat.a_t_less_s.col(t - 1))) *
          (Z_t.col(i) - y_bar).t();
        P_v_v += (w_i * (Z_t.col(i) - y_bar)) * (Z_t.col(i) - y_bar).t() +
          arma::diagmat(w_i * Z_t.col(i) % (1 - Z_t.col(i)));
      }

      P_v_v.diag() += p_dat.denom_term;

      // Compute new estimates
      inv_sympd(P_v_v, P_v_v, p_dat.use_pinv, "ddhazard_fit_cpp estimation error: Failed to invert P_v_v");

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.col(t - 1) +
        P_a_v * (P_v_v * ((p_dat.is_event_in_bin(r_set) == t - 1) - y_bar));

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) - P_a_v * P_v_v * P_a_v.t();

      p_dat.B_s.slice(t - 1) =
        arma::solve(p_dat.V_t_less_s.slice(t - 1),
                    p_dat.state_trans->map(
                      p_dat.V_t_t_s.slice(t - 1), right).sv).t();
    }
}



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
UKF_solver_New::UKF_solver_New(
  ddhazard_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
  Rcpp::Nullable<Rcpp::NumericVector> &alpha,
  Rcpp::Nullable<Rcpp::NumericVector> &beta,
  family_base &fam):
  p_dat(p_),
  m(p_.a_t_t_s.n_rows),

  a(!alpha.isNull() ? Rcpp::as< Rcpp::NumericVector >(alpha)[0] : 1e-1),
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

  sigma_points(arma::mat(m, 2 * m + 1)),
  fam(fam)
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
    my_debug_logger(p_dat)
      << "alpha, beta, kappa, n_dim = "
      << a << ", "
      << b << ", "
      << k << ", "
      << m;
    my_debug_logger(p_dat)
      << "w_0, w_0_c, w_0_cc, w_i, lambda = "
      << w_0 << ", "
      << w_0_c << ", "
      << w_0_cc << ", "
      << w_i << ", "
      << lambda;
  }
}

void UKF_solver_New::compute_sigma_points(
    const arma::vec &a_t, arma::mat &s_points, const arma::mat &P_x_x){
  arma::mat cholesky_decomp;
  if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
    Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
  }

  s_points.col(0) = a_t;
  for(arma::uword i = 1; i < s_points.n_cols; ++i)
    if(i % 2 == 0)
      s_points.col(i) = a_t + sqrt_m_lambda * cholesky_decomp.col((i - 1) / 2); else
        s_points.col(i) = a_t - sqrt_m_lambda * cholesky_decomp.col((i - 1) / 2);
}

void UKF_solver_New::solve(){
  double bin_stop = p_dat.min_start;
  const bool uses_at_risk_length = fam.uses_at_risk_length();

  for (int t = 1; t < p_dat.d + 1; t++){
    double bin_start = bin_stop;
    double delta_t = p_dat.I_len[t - 1];
    bin_stop += delta_t;

    // E-step: Prediction
    p_dat.a_t_less_s.col(t - 1) =
      p_dat.state_trans->map(p_dat.a_t_t_s.col(t - 1)).sv;
    p_dat.V_t_less_s.slice(t - 1) =
      p_dat.state_trans->map(p_dat.V_t_t_s.slice(t - 1)).sv +
      delta_t * p_dat.err_state->map(p_dat.Q).sv;

    // Re-generate
    if(p_dat.debug){
      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1),
               "Chol decomposing for regenerations:");
    }

    compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                         sigma_points, p_dat.V_t_less_s.slice(t - 1));

    if(p_dat.debug){
      my_print(p_dat, sigma_points, "new sigma points");
    }

    // E-step: correction-step
    arma::uvec r_set = get_risk_set(p_dat, t);

    // ** 1: compute means and variances **
    arma::uword n_risk = r_set.n_elem;
    arma::vec sqrt_weights_to_sds(n_risk, arma::fill::zeros);
    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat O(n_risk, sigma_points.n_cols);
    for(arma::uword i = 0; i < O.n_cols; ++i){
      O.col(i) = (p_dat.state_lp->map(sigma_points.col(i)).sv.t() *
        p_dat.X.cols(r_set)).t();
    }
    O.each_col() += p_dat.fixed_effects(r_set);

    arma::ivec do_die =
      arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);
    arma::vec at_risk_length(n_risk);
    arma::vec starts = p_dat.tstart(r_set);
    arma::vec stops = p_dat.tstop(r_set);

    if(uses_at_risk_length){
      for(arma::uword i = 0; i < n_risk; i++){

          at_risk_length[i] = get_at_risk_length(
            stops[i], bin_stop, starts[i], bin_start);
      }
    } else
      at_risk_length.zeros();

    for(arma::uword j = 0; j < sigma_points.n_cols; ++j){
      double w = (j == 0) ? w_0 : w_i;
      double w_c = (j == 0) ? w_0_c : w_i;

      const int *die_i = do_die.memptr();
      const double *at_i = at_risk_length.memptr();
      double *O_i_j = O.colptr(j);
      double *sqrt_w_i = sqrt_weights_to_sds.memptr();
      for(arma::uword i = 0; i < n_risk;
          i++, die_i++, at_i++, O_i_j++, sqrt_w_i++){
        auto trunc_eta = fam.truncate_eta(*die_i, *O_i_j, exp(*O_i_j), *at_i);

        *O_i_j = fam.linkinv(trunc_eta, *at_i);
        *sqrt_w_i += w_c * fam.var(trunc_eta, *at_i);
      }

      y_bar += w * O.col(j);
    }

    sqrt_weights_to_sds += p_dat.denom_term;
    sqrt_weights_to_sds = p_dat.weights(r_set) / sqrt_weights_to_sds;
    sqrt_weights_to_sds.for_each([](double &val) { val = std::sqrt(val); } );

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;
    O = (O.each_col() % sqrt_weights_to_sds).t();

    arma::vec c_vec = O * (sqrt_weights_to_sds % (do_die - y_bar));
    O = O * O.t();

    // Compute intermediate matrix
    arma::mat tmp_mat = arma::diagmat(weights_vec_c_inv) + O;

    // Compute vector for state space vector update
    c_vec = c_vec -  O * arma::solve(tmp_mat, c_vec);

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    O = O - O * arma::solve(tmp_mat, O);

    // Substract mean to get delta sigma points
    arma::mat delta_sigma_points = sigma_points.each_col()
      - p_dat.a_t_less_s.col(t - 1);

    p_dat.a_t_t_s.col(t) =
      p_dat.a_t_less_s.col(t - 1) +
      p_dat.LR * delta_sigma_points * (weights_vec_cc % c_vec);

    p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) -
      (delta_sigma_points.each_row() % weights_vec_cc.t()) * O * (delta_sigma_points.each_row() % weights_vec_cc.t()).t();

    if(p_dat.debug){
      std::stringstream str, str_less;
      str_less << t << "|" << t - 1;
      str << t << "|" << t;

      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1).diag(), "diag(V_(" + str_less.str() + "))");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str_less.str() + ") is "
        << arma::cond(p_dat.V_t_less_s.slice(t - 1));

      my_print(p_dat, p_dat.V_t_t_s.slice(t).diag(), "diag(V_(" + str.str()  + "))");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_t_s.slice(t));

      my_print(p_dat, p_dat.a_t_less_s.col(t - 1), "a_(" + str_less.str() + ")");
      my_print(p_dat, p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
    }

    // We are looking at:
    //  X = B A^-1
    // X^T = A^-1 B^T <=> A X^T = B^T
    p_dat.B_s.slice(t - 1) = arma::solve(
      p_dat.V_t_less_s.slice(t - 1),
      p_dat.state_trans->map(p_dat.V_t_t_s.slice(t - 1), left).sv).t();
  }
}
