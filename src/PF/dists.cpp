#include "dists.h"
#include "dmvnrm.h"
#include "../utils.h"
#include "../family.h"

state_fw::state_fw(
  const arma::vec &parent, const arma::mat &F,const covarmat &Q):
  parent(parent), Q(Q),
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

bool state_fw::is_grad_z_hes_const() const {
  return false;
}

const arma::vec& state_fw::get_mean() const {
  return mu;
}

arma::uword state_fw::dim() const {
  return parent.n_elem;
}

double state_fw::log_dens(const arma::vec &child) const {
  return dmvnrm_log(child, mu, Q.chol_inv());
}

arma::vec state_fw::gradient(const arma::vec &child) const {
  /* TODO: maybe exploit that we have decomposition of Q */
  return arma::solve(Q.mat(), mu - child);
}

arma::vec state_fw::gradient_zero(const arma::vec *pa) const {
  return QiF * *pa;
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

bool state_bw::is_grad_z_hes_const() const {
  return false;
}

const arma::vec& state_bw::get_mean() const {
  Rcpp::stop("'state_bw' is not a density in x");

  throw std::logic_error("this code will should not be reached");
}

arma::uword state_bw::dim() const {
  return child.n_elem;
}

double state_bw::log_dens(const arma::vec &parent) const {
  return dmvnrm_log(child, F * parent, Q.chol_inv());
}

arma::vec state_bw::gradient(const arma::vec &parent) const {
  return FtQi * (child - F * parent);
}

arma::vec state_bw::gradient_zero(const arma::vec *ci) const {
  return FtQi * *ci;
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

bool artificial_prior::is_grad_z_hes_const() const {
  return true;
}

const arma::vec& artificial_prior::get_mean() const {
  return mut;
}

arma::uword artificial_prior::dim() const {
  return mut.n_elem;
}

double artificial_prior::log_dens(const arma::vec &state) const {
  return dmvnrm_log(state, mut, Qt.chol_inv());
}

arma::vec artificial_prior::gradient(const arma::vec &state) const {
  /* TODO: maybe exploit that we have decomposition of Q */
  return arma::solve(Qt.mat(), mut - state);
}

arma::vec artificial_prior::gradient_zero(const arma::vec *state) const {
  return dz;
}

arma::mat artificial_prior::neg_Hessian(const arma::vec &state) const {
  return Qt.inv();
}



artificial_prior_generator::artificial_prior_generator(
  const arma::mat &F, const covarmat &Q, const arma::vec &mu_0,
  const covarmat &Q_0): F(F), Q(Q) {
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



template<class T>
class observational_cdist :
  public virtual PF_cdist,  public virtual T
{
  const arma::mat X;
  const arma::uvec is_event;
  const arma::vec offsets;
  const arma::vec tstart;
  const arma::vec tstop;
  const double bin_start;
  const double bin_stop;
  const bool multithreaded;
  const arma::vec at_risk_length;

public:
  observational_cdist(
    const arma::mat&, const arma::uvec&, const arma::vec&,
    const arma::vec&, const arma::vec&, const double, const double,
    const bool multithreaded = false);
  ~observational_cdist() = default;

  bool is_mvn() const override;
  bool is_grad_z_hes_const() const override;
  const arma::vec& get_mean() const override;
  arma::uword dim() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec*) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

arma::vec set_at_risk_length(
    const arma::vec &tstart, const arma::vec &stop,
    double bin_start, double bin_stop,
    const bool uses_at_risk_length) {
  if(!uses_at_risk_length)
    return arma::vec(tstart.n_elem, arma::fill::zeros);

  arma::vec out(tstart.n_elem);
  double *o = out.begin();
  const double *sta = tstart.begin(), *sto = stop.begin();
  for(arma::uword i = 0; i < tstart.n_elem; ++i)
    *(o++) = get_at_risk_length(*(sto++), bin_stop, *(sta++), bin_start);

  return out;
}

#ifdef _OPENMP
/* openMP reductions */
#pragma omp declare reduction(armaVP: arma::vec: omp_out += omp_in)
#pragma omp declare reduction(armaMP: arma::mat: omp_out += omp_in)
#endif

template<class T>
observational_cdist<T>::observational_cdist(
  const arma::mat &X, const arma::uvec &is_event,
  const arma::vec &offsets, const arma::vec &tstart, const arma::vec &tstop,
  const double bin_start, const double bin_stop, const bool multithreaded):
  X(X), is_event(is_event), offsets(offsets), tstart(tstart),
  tstop(tstop), bin_start(bin_start), bin_stop(bin_stop),
  multithreaded(multithreaded),
  at_risk_length(set_at_risk_length(tstart, tstop, bin_start, bin_stop,
                                    this->uses_at_risk_length())) { }

template<class T>
bool observational_cdist<T>::is_mvn() const {
  return false;
}

template<class T>
bool observational_cdist<T>::is_grad_z_hes_const() const {
  return false;
}

template<class T>
const arma::vec& observational_cdist<T>::get_mean() const {
  Rcpp::stop("'observational_cdist<T>' is not a density in x");

  throw std::logic_error("this code will should not be reached");
}

template<class T>
arma::uword observational_cdist<T>::dim() const {
  return X.n_rows;
}

template<class T>
double observational_cdist<T>::log_dens(const arma::vec &coefs) const {
  const arma::vec eta = X.t() * coefs + offsets;
  arma::uword n = eta.n_elem;

  if(n < 1L)
    return 0.;

  /* compute log likelihood */
  double result = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:result) \
  if(multithreaded)
#endif
  for(unsigned int i = 0; i < n; ++i){
    auto trunc_eta = this->truncate_eta(
      is_event[i], eta[i], exp(eta[i]), at_risk_length[i]);
    result += this->log_like(is_event[i], trunc_eta, at_risk_length[i]);
  }

  return result;
}

template<class T>
arma::vec observational_cdist<T>::gradient(const arma::vec &coefs) const {
  const arma::vec eta = X.t() * coefs + offsets;
  const arma::uword n = eta.n_elem;

  /* compute gradient */
  arma::vec result(coefs.n_elem, arma::fill::zeros);
  if(n < 1L)
    return result;

#ifdef _OPENMP
  bool first_it = true;
#pragma omp parallel for schedule(static) if(multithreaded)\
  reduction(armaVP:result) firstprivate(first_it)
#endif
  for(unsigned i = 0; i < n; ++i){
#ifdef _OPENMP
    if(first_it){
      result.zeros(coefs.n_elem);
      first_it = false;
    }
#endif
    auto trunc_eta = this->truncate_eta(
      is_event[i], eta[i], exp(eta[i]), at_risk_length[i]);
    double d_l = this->d_log_like(is_event[i], trunc_eta, at_risk_length[i]);
    result += X.col(i) * d_l;
  }

  return result;
}

template<class T>
arma::vec observational_cdist<T>::gradient_zero(const arma::vec *coefs) const {
  Rcpp::stop("'observational_cdist<T>::gradient' is not implemented");

  throw std::logic_error("this code will should not be reached");
}

template<class T>
arma::mat observational_cdist<T>::neg_Hessian(const arma::vec &coefs) const {
  const arma::vec eta = X.t() * coefs  + offsets;
  const arma::uword n = eta.n_elem;

  /* compute Hessian */
  arma::mat result(coefs.n_elem, coefs.n_elem, arma::fill::zeros);

  if(n < 1L)
    return result;

#ifdef _OPENMP
  bool first_it = true;
#pragma omp parallel for schedule(static) reduction(armaMP:result) \
  if(multithreaded) firstprivate(first_it)
#endif
  for(unsigned int i = 0; i < n; ++i){
#ifdef _OPENMP
    if(first_it){
      result.zeros(coefs.n_elem, coefs.n_elem);
      first_it = false;
    }
#endif
    auto trunc_eta = this->truncate_eta(
      is_event[i], eta[i], exp(eta[i]), at_risk_length[i]);
    double dd_l = this->dd_log_like(
      is_event[i], trunc_eta, at_risk_length[i]);
    sym_mat_rank_one_update(dd_l, X.col(i), result);
  }

  return -arma::symmatu(result);
}

std::shared_ptr<PF_cdist> get_observational_cdist(
    const std::string& fam, const arma::mat &X,
    const arma::uvec &is_event, const arma::vec &offsets,
    const arma::vec &tstart, const arma::vec &tstop, const double bin_start,
    const double bin_stop, const bool multithreaded){
  if(fam == BINOMIAL)
    return(std::shared_ptr<PF_cdist>((new observational_cdist<logistic>(
        X, is_event, offsets, tstart, tstop, bin_start, bin_stop,
        multithreaded))));
  else if(fam == CLOGLOG)
    return(std::shared_ptr<PF_cdist>((new observational_cdist<cloglog>(
        X, is_event, offsets, tstart, tstop, bin_start, bin_stop,
        multithreaded))));
  else if(fam == POISSON)
    return(std::shared_ptr<PF_cdist>((new observational_cdist<exponential>(
        X, is_event, offsets, tstart, tstop, bin_start, bin_stop,
        multithreaded))));

  Rcpp::stop("'fam' not implemented");
  throw std::logic_error("this code will should not be reached");
}
