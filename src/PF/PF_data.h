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
    } else{
      elapsed_seconds = std::numeric_limits<double>::quiet_NaN();

    }

    last_message_time.reset(now.release());

#ifdef _OPENMP
    omp_unset_lock(&lock);
#endif

    return elapsed_seconds;
  }
};

/* data class for pre-computed factorization */
class covarmat {
public:
  const arma::mat mat;
  const arma::mat chol;
  const arma::mat chol_inv;
  const arma::mat inv;

  template<typename T>
  covarmat(T Q):
    mat(Q), chol(arma::chol(mat)), chol_inv(arma::inv(arma::trimatu(chol))), inv(arma::inv(Q))
  {}
};

// data holder for particle filtering
class PF_data : public problem_data {
public:
  /* Number of paprticles in forward and/or backward filter */
  const arma::uword N_fw_n_bw;
  const arma::uword N_smooth;
  const double forward_backward_ESS_threshold;

  /* Inital state, number of particles to draw at time 0 and d + 1 and debug level */
  const arma::vec &a_0;
  const unsigned int debug; /* < 1 is no info and greater values yields more info */
  const arma::uword N_first;
  const unsigned long work_block_size;

  /* pre-computed factorization */
  const covarmat Q;
  const covarmat Q_0;
  const covarmat Q_proposal;
  const covarmat Q_proposal_smooth;

  PF_data(const int n_fixed_terms_in_state_vec,
          arma::mat &X,
          arma::mat &fixed_terms,
          const arma::vec &tstart,
          const arma::vec &tstop, const arma::ivec &is_event_in_bin,
          const arma::colvec &a_0,
          arma::mat &Q_0,
          arma::mat &Q,
          const Rcpp::List &risk_obj,
          const arma::mat &F,
          const int n_max,
          const int n_threads,

          // new arguments
          const arma::mat Q_tilde,
          const arma::uword N_fw_n_bw,
          const arma::uword N_smooth,
          Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
          const unsigned int debug,
          const arma::uword N_first):
    problem_data(
      n_fixed_terms_in_state_vec,
      X,
      fixed_terms,
      tstart,
      tstop, is_event_in_bin,
      a_0,
      Q_0,
      Q,
      risk_obj,
      F,
      n_max,
      n_threads),

      N_fw_n_bw(N_fw_n_bw),
      N_smooth(N_smooth),
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
      Q_proposal(Q + Q_tilde),
      Q_proposal_smooth((Q + Q_tilde) * .5)
    {
#ifdef _OPENMP
    omp_init_lock(&PF_logger::lock);
#endif
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
};

#endif
