#include "dists.h"
#include "dmvnrm.h"

state_fw::state_fw(
  const arma::vec &parent, const arma::mat &F,const covarmat &Q):
  parent(parent), F(F), Q(Q),
  /* TODO: maybe exploit that we have decomposition of Q */
  QiF(arma::solve(Q.mat(), F)), mu(F * parent) { }

double state_fw::log_dens_func(
    const arma::vec &child, const arma::vec &pa,
    const arma::mat &Fm, const covarmat &Qm){
  return dmvnrm_log(child, Fm * pa, Qm.chol_inv());
}

bool state_fw::is_mvn() const {
  return true;
}

double state_fw::log_dens(const arma::vec &child) const {
  return dmvnrm_log(child, mu, Q.chol_inv());
}

arma::vec state_fw::gradient(const arma::vec &child) const {
  /* TODO: maybe exploit that we have decomposition of Q */
  return arma::solve(Q.mat(), mu - child);
}

arma::vec state_fw::gradient_zero(const arma::vec &pa) const {
  return QiF * pa;
}

arma::mat state_fw::neg_Hessian(const arma::vec &child) const {
  return Q.inv();
}



state_bw::state_bw(
  const arma::vec &child, const arma::mat &F, const covarmat &Q):
  child(child), F(F), Q(Q),
  /* TODO: maybe exploit that we have decomposition of Q */
  FtQi(arma::solve(Q.mat(), F).t()), n_hes(FtQi * F) {}

double state_bw::log_dens_func(
    const arma::vec &parent, const arma::vec &ch,
    const arma::mat &Fm, const covarmat &Qm){
  return dmvnrm_log(ch, Fm * parent, Qm.chol_inv());
}

bool state_bw::is_mvn() const {
  return true;
}

double state_bw::log_dens(const arma::vec &parent) const {
  return dmvnrm_log(child, F * parent, Q.chol_inv());
}

arma::vec state_bw::gradient(const arma::vec &parent) const {
  return FtQi * (child - F * parent);
}

arma::vec state_bw::gradient_zero(const arma::vec &ci) const {
  return FtQi * ci;
}

arma::mat state_bw::neg_Hessian(const arma::vec &parent) const {
  return n_hes;
}



artificial_prior::artificial_prior(
  const arma::vec &mut, const covarmat &Qt): mut(mut), Qt(Qt),
  /* TODO: maybe exploit that we have decomposition of Q */
  dz(arma::solve(Qt.mat(), mut)) { }

bool artificial_prior::is_mvn() const {
  return true;
}

double artificial_prior::log_dens(const arma::vec &state) const {
  return dmvnrm_log(state, mut, Qt.chol_inv());
}

arma::vec artificial_prior::gradient(const arma::vec &state) const {
  /* TODO: maybe exploit that we have decomposition of Q */
  return arma::solve(Qt.mat(), mut - state);
}

arma::vec artificial_prior::gradient_zero(const arma::vec &state) const {
  return dz;
}

arma::mat artificial_prior::neg_Hessian(const arma::vec &state) const {
  return Qt.inv();
}



artificial_prior_generator::artificial_prior_generator(
  const arma::mat &F, const covarmat &Q, const arma::vec &mu_0,
  const covarmat &Q_0): F(F), Q(Q), mu_0(mu_0), Q_0(Q_0) {
  mt.insert(std::make_pair(0L, mu_0));
  Pt.insert(std::make_pair(0L, covarmat(Q_0)));
}

artificial_prior artificial_prior_generator::get_artificial_prior(
    const arma::uword time){
  if(Pt.find(time) == Pt.end()){
    /* need to add elements up to `time` */
    auto it_p = (--Pt.end()); /* must have a value due to constructor */
    auto it_m = (--mt.end());

    if(it_p->first != it_m->first)
      Rcpp::stop("indices does not match in 'artificial_prior_generator::get_artificial_prior'");

    for(arma::uword j = it_p->first + 1L; j <= time; ++j, ++it_p, it_m++){
      Pt.insert(std::make_pair(
          j, covarmat(F * it_p->second.mat() * F.t() + Q.mat())));
      mt.insert(std::make_pair(j, F * it_m->second));
    }
  }

  return artificial_prior(mt.find(time)->second, Pt.find(time)->second);
}
