#include "densities.h"

pf_dens::pf_dens(const PF_data &data, const std::string &family):
  family(family), data(data), art_gen(data.F, data.Q, data.a_0, data.Q_0)
{}

double pf_dens::log_prob_state_given_parent(
    const arma::vec &child, const arma::vec &parent){
  return state_fw::log_dens_func(child, parent, data.F, data.Q);
}

double pf_dens::log_prob_state_given_child(
    const arma::vec &parent, const arma::vec &child){
  return state_bw::log_dens_func(parent, child, data.F, data.Q);
}

std::unique_ptr<PF_cdist> pf_dens::get_fw_dist
  (const arma::vec &parent){
  return std::unique_ptr<state_fw>(new state_fw(
      parent, data.F, data.Q));
}

std::unique_ptr<PF_cdist> pf_dens::get_bw_dist
  (const arma::vec &child){
  return std::unique_ptr<state_bw>(new state_bw(child, data.F, data.Q));
}

std::shared_ptr<PF_cdist> pf_dens::get_prior(const arma::uword t){
  return std::make_shared<artificial_prior>(art_gen.get_artificial_prior(t));
}

std::shared_ptr<PF_cdist> pf_dens::get_y_dist(
    const int t, const bool multithreaded){
  arma::uvec r_set = get_risk_set(data.risk_sets, t);
  /* zero indexed while t is not */
  arma::uvec is_event = data.is_event_in_bin(r_set) == t - 1L;
  auto bin_start_stop = get_bin_times(data.min_start, data.I_len, t);

  return get_observational_cdist(
    family, data.X.cols(r_set), std::move(is_event),
    data.fixed_effects(r_set),  data.tstart(r_set), data.tstop(r_set),
    bin_start_stop.start, bin_start_stop.stop, multithreaded);
}

std::string get_family(const std::string &model){
  if(model == "logit" || model == "binomial"){
    return BINOMIAL;

  } else if (model == "exponential" || model == "poisson"){
    return POISSON;

  } else if (model == "cloglog") {
    return CLOGLOG;

  }

  std::stringstream stream;
  stream << "model '" << model << "' is not implemented";
  Rcpp::stop(stream.str());

  return "";
}
