#ifndef PF_DATA
#define PF_DATA

#include "../arma_n_rcpp.h"
#include "../problem_data.h"
#include <chrono>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Logger class for debug information */
class PF_logger {
public:
#ifdef _OPENMP
  static omp_lock_t lock; /* NOTICE: have to call omp_init_lock somewhere
                             before usage!                                */
#endif

  PF_logger(const bool log, const unsigned int level);

  PF_logger(PF_logger&& other);

  ~PF_logger();

  PF_logger(const PF_logger&) = delete;
  PF_logger& operator=(const PF_logger&) = delete;

  /*
    Class to prefix stream. See:
      https://stackoverflow.com/a/27336473
  */

  class prefixbuf : public std::streambuf
  {
    std::string     prefix;
    std::streambuf* sbuf;
    bool            need_prefix;

    int sync() {
      return this->sbuf->pubsync();
    }

    int overflow(int c) {
      if (c != std::char_traits<char>::eof()) {
        if (this->need_prefix
              && !this->prefix.empty()
              && this->prefix.size() != (unsigned int)this->sbuf->sputn(&this->prefix[0], this->prefix.size())) {
              return std::char_traits<char>::eof();
        }
        this->need_prefix = c == '\n';
      }
      return this->sbuf->sputc(c);
    }

  public:
    prefixbuf(std::string const& prefix, std::streambuf* sbuf)
      : prefix(prefix), sbuf(sbuf) , need_prefix(true) {}
  };

  class oprefixstream
    : private virtual prefixbuf, public std::ostream
  {
  public:
    oprefixstream(std::string const& prefix, std::ostream& out)
      : prefixbuf(prefix, out.rdbuf())
    , std::ios(static_cast<std::streambuf*>(this))
    , std::ostream(static_cast<std::streambuf*>(this)) {}
  };

  template<typename T>
  std::ostream& operator<<(const T &obj){
    if(log){
      if(!os_w_prefix){
        os_w_prefix.reset(new oprefixstream(get_prefix(level), os));
      }

      return *os_w_prefix << obj;
    }

    std::ostringstream tmp;
    if(!oprefixstream_dummy){
      oprefixstream_dummy.reset(new oprefixstream("", tmp));
    }
    return *oprefixstream_dummy;
  }

private:
  bool log;
  unsigned int level;
  std::ostringstream os;
  std::unique_ptr<oprefixstream> os_w_prefix;
  std::unique_ptr<oprefixstream> oprefixstream_dummy;

  const static unsigned int n_spaces = 3;
  static std::string get_prefix(const unsigned int level);

  using tp = std::chrono::time_point<std::chrono::system_clock>;
  using tp_pointer =std::unique_ptr<tp>;

  static tp_pointer last_message_time;

  static double get_elapsed_seconds_n_set_last_message_time(){
    double elapsed_seconds;

#ifdef _OPENMP
    omp_set_lock(&lock);
#endif
    tp_pointer now(new tp());
    *now = std::chrono::system_clock::now();

    if(last_message_time){
      std::chrono::duration<double> tmp = *now - *last_message_time;
      elapsed_seconds = tmp.count();
    } else {
      elapsed_seconds = std::numeric_limits<double>::quiet_NaN();

    }

    last_message_time.reset(now.release());

#ifdef _OPENMP
    omp_unset_lock(&lock);
#endif

    return elapsed_seconds;
  }
};

/* data class for pre-computed factorization and matrices which are only
 * computed once */
class covarmat {
private:
#ifdef _OPENMP
  std::unique_ptr<omp_lock_t> lock =
    std::unique_ptr<omp_lock_t>((new omp_lock_t()));
  class LockGuard {
  public:
    explicit LockGuard(omp_lock_t& lock) : m_lock(lock){
      omp_set_lock(&m_lock);
    }

    ~LockGuard() {
      omp_unset_lock(&m_lock);
    }

  private:
    omp_lock_t& m_lock;
  };
#endif

  enum output { e_mat, e_chol, e_chol_inv, e_inv };

  std::unique_ptr<const arma::mat> mat_;
  std::unique_ptr<bool> is_chol_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_;
  std::unique_ptr<bool> is_chol_inv_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_inv_;
  std::unique_ptr<bool> is_inv_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> inv_;

  const arma::mat& get_mat(output what) const {
    if(what == e_mat)
      return(*mat_.get());

    bool is_computed, *this_flag;

    this_flag = is_chol_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *chol_.get() += arma::chol(*mat_.get());
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *chol_.get() += arma::chol(*mat_.get());
      *this_flag = true;
    }
#endif
    if(what == e_chol)
      return(*chol_.get());


    this_flag = is_chol_inv_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *chol_inv_.get() += arma::inv(arma::trimatu(*chol_.get()));
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *chol_inv_.get() += arma::inv(arma::trimatu(*chol_.get()));
      *this_flag = true;
    }
#endif
    if(what == e_chol_inv)
      return(*chol_inv_.get());


    this_flag = is_inv_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *inv_.get() += mat_->i();
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *inv_.get() += mat_->i();
      *this_flag = true;
    }
#endif
    return(*inv_.get());
  }

public:
  /* Covar matrix: Q */
  const arma::mat& mat() const  {
    return get_mat(e_mat);
  }
  /* C in Q = C^\topC in  */
  const arma::mat& chol() const {
    return get_mat(e_chol);
  }
  /* C^{-1} */
  const arma::mat& chol_inv() const {
    return get_mat(e_chol_inv);
  }
  /* Q^{-1} */
  const arma::mat& inv() const {
    return get_mat(e_inv);
  }

  template<typename T>
  covarmat(T Q):
    mat_(new arma::mat(Q)),
    chol_    (new arma::mat(arma::size(Q), arma::fill::zeros)),
    chol_inv_(new arma::mat(arma::size(Q), arma::fill::zeros)),
    inv_     (new arma::mat(arma::size(Q), arma::fill::zeros)) {
#ifdef _OPENMP
    omp_init_lock(lock.get());
#endif
  }

  covarmat(const covarmat &other): covarmat(other.mat()) { }

  ~covarmat(){
#ifdef _OPENMP
    omp_destroy_lock(lock.get());
#endif
  }
};

// data holder for particle filtering
class PF_data : public problem_data {
  using uword = arma::uword;

  /* Objects used to compute P(\alpha_t \vert \alpha_{t - 1}) */
  std::map<const uword, const std::unique_ptr<linear_mapper>> bw_mean_maps;
  std::map<const uword, const arma::vec> bw_mean_const_term;
  std::map<const uword, const covarmat> bw_covar_map;

  /* Pre-computed objects for the artificial prior */
  std::map<const uword, const arma::vec> uncond_means;
  std::map<const uword, const covarmat>  uncond_covarmats;

  /* uncondtional mean and covariance matrix terms in backward filter's
   * and smoother's proposal distribution in the Taylor approximation */
  std::map<const uword, const arma::vec> uncond_mean_terms;
  std::map<const uword, const arma::mat> uncond_covs_inv;

protected:
  virtual void set_maps(){
    arma::vec m_t, m_t_p_1 = a_0;
    arma::mat Q_state = err_state->map(Q.mat()).sv;
    arma::mat P_t, P_t_p_1 = Q_0.mat();
    for(int t = 0; t <= d + 1; ++t){

      // move and update
      m_t = std::move(m_t_p_1);
      P_t = std::move(P_t_p_1);
      m_t_p_1 = state_trans->map(m_t).sv;
      P_t_p_1 = state_trans->map(P_t).sv + Q_state;

      // insert map elements
      uncond_means.insert(std::make_pair(t + 1L, m_t_p_1));
      uncond_covarmats.insert(std::make_pair(t + 1L, covarmat(P_t_p_1)));

      arma::mat S_t =
        arma::solve(P_t_p_1, Q_state);
      S_t = state_trans_inv->map(S_t, right).sv;
      S_t = state_trans->map(P_t, right).sv * S_t;

      /* 1. P_{t +1}^{-\top}F
       * 2. P_t(P_{t +1}^{-\top}F)^\top */
      arma::mat tmp = arma::solve(P_t_p_1.t(), state_trans->map());
      bw_mean_maps.insert(std::make_pair(
          t, std::unique_ptr<dens_mapper>(
              new dens_mapper(arma::mat(P_t * tmp.t())))));

      bw_mean_const_term.insert(std::make_pair(
          t, S_t * arma::solve(P_t, m_t)));

      bw_covar_map.insert(std::make_pair(t, std::move(S_t)));

      uncond_mean_terms.insert(std::make_pair(
          t, arma::vec(arma::solve(P_t, m_t))));

      uncond_covs_inv.insert(std::make_pair(
          t, arma::mat(P_t.i())));

      if(debug > 4)
        log(5) << "P_" << t << std::endl << P_t
               << "m_" << t << std::endl << m_t.t()
               << "S_" << t << std::endl << bw_covar_map.at(t).mat()
               << "P(a_" << t << " | a_" << t  + 1 << ") mean term"
               << std::endl << uncond_mean_terms.at(t).t();
    }
  }

public:
  /* Number of paprticles in forward and/or backward filter */
  const uword N_fw_n_bw;
  const uword N_smooth;
  const uword N_smooth_final;
  const double forward_backward_ESS_threshold;

  /* Inital state, number of particles to draw at time 0 and d + 1 and debug level */
  const arma::vec &a_0;
  const unsigned int debug; /* < 1 is no info and greater values yields more info */
  const uword N_first;
  const unsigned long work_block_size;

  /* pre-computed factorization */
  const covarmat Q;
  const covarmat Q_0;
  const covarmat Q_proposal;
  const covarmat Q_proposal_state;

  PF_data(const int n_fixed_terms_in_state_vec,
          arma::mat &X,
          arma::mat &fixed_terms,
          const arma::vec &tstart,
          const arma::vec &tstop, const arma::ivec &is_event_in_bin,
          const arma::colvec &a_0,
          const arma::mat &R,
          const arma::mat &L,
          arma::mat &Q_0,
          arma::mat &Q,
          const Rcpp::List &risk_obj,
          const arma::mat &F,
          const int n_max,
          const int n_threads,
          const arma::vec &fixed_parems,

          // new arguments
          const arma::mat Q_tilde,
          const uword N_fw_n_bw,
          const uword N_smooth,
          const uword N_smooth_final,
          Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
          const unsigned int debug,
          const uword N_first) :
    problem_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin, a_0, R, L, Q_0, Q,
      risk_obj, F, n_max, n_threads, fixed_parems),

      N_fw_n_bw(N_fw_n_bw), N_smooth(N_smooth), N_smooth_final(N_smooth_final),
      forward_backward_ESS_threshold(
        forward_backward_ESS_threshold.isNotNull() ?
          Rcpp::as<Rcpp::NumericVector>(forward_backward_ESS_threshold)[0] :
          N_fw_n_bw / 2.),

      a_0(a_0),
      debug(debug),
      N_first(N_first),
      work_block_size(500),

      Q(Q),
      Q_0(Q_0),
      Q_proposal(Q_tilde),
      Q_proposal_state(err_state->map(Q_tilde).sv)
    {
#ifdef _OPENMP
      omp_init_lock(&PF_logger::lock);
#endif
      set_maps();
    }

  ~PF_data(){
#ifdef _OPENMP
    omp_destroy_lock(&PF_logger::lock);
#endif
  }

  PF_logger log(const unsigned int level) const{
    return PF_logger(level <= debug, level);
  }

  PF_data & operator=(const PF_data&) = delete;
  PF_data(const PF_data&) = delete;
  PF_data() = delete;

  const arma::vec bw_mean(const uword t, const arma::vec &a) const {
    return bw_mean_maps.at(t)->map(a).sv + bw_mean_const_term.at(t);
  }

  const covarmat& bw_covar(const uword t) const {
    return bw_covar_map.at(t);
  }

  const arma::vec& uncond_mean_state(const uword t) const {
    return uncond_means.at(t);
  }

  const covarmat& uncond_covar_state(const uword t) const {
    return uncond_covarmats.at(t);
  }

  const arma::vec& uncond_mean_term(const int t) const {
    return uncond_mean_terms.at(t);
  }

  const arma::mat& uncond_covar_inv(const uword t) const {
    return uncond_covs_inv.at(t);
  }
};

#endif
