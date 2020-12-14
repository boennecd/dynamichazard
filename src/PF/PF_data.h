#ifndef PF_DATA
#define PF_DATA

#include "../arma_n_rcpp.h"
#include "covarmat.h"
#include <thread>

/* Logger class for debug information */
class tracker;
class PF_logger {
public:
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

    int sync();

    int overflow(int);

  public:
    prefixbuf(std::string const&, std::streambuf*);
  };

  class oprefixstream
    : private virtual prefixbuf, public std::ostream
  {
  public:
    oprefixstream(std::string const&, std::ostream&);
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

  static std::unique_ptr<tracker> last_message_time;

  static double get_elapsed_seconds_n_set_last_message_time();
};

// data holder for particle filtering
class PF_data {
  using uword = arma::uword;
public:
  const int d;
  const Rcpp::List risk_sets;

  /* these are not const due the arma::mat constructor which does not allow
  * copy_aux_mem = false with const pointer. They should not be changed
  * though...                                                              */
  arma::mat X;
  arma::mat fixed_terms;

  const std::vector<double> I_len;
  const arma::mat F;
  const arma::uword state_dim;

  const int n_threads;

  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;

  const double min_start;

  arma::vec fixed_params;
  arma::vec fixed_effects;

  /* Number of paprticles in forward and/or backward filter */
  const uword N_fw_n_bw;
  const uword N_smooth;
  const uword N_smooth_final;
  const double forward_backward_ESS_threshold;

  /* Inital state, number of particles to draw at time 0 and d + 1 and debug level */
  const arma::vec &a_0;
  const unsigned int debug; /* < 1 is no info and greater values yields more info */
  const uword N_first;
  const int nu;
  const unsigned long work_block_size;
  const double covar_fac, ftol_rel;

  /* pre-computed factorization */
  const covarmat Q;
  const covarmat Q_0;
  const arma::mat xtra_covar;

  PF_data(const int n_fixed_terms_in_state_vec,
          arma::mat &X,
          arma::mat &fixed_terms,
          const arma::vec &tstart,
          const arma::vec &tstop, const arma::ivec &is_event_in_bin,
          const arma::vec &a_0,
          const arma::mat &R,
          const arma::mat &L,
          arma::mat &Q_0,
          arma::mat &Q,
          const Rcpp::List &risk_obj,
          const arma::mat &F,
          const int n_max,
          const int n_threads,
          const arma::vec &fixed_params,

          // new arguments
          const arma::mat Q_tilde,
          const uword N_fw_n_bw,
          const uword N_smooth,
          const uword N_smooth_final,
          Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
          const unsigned int debug,
          const uword N_first, const int nu, const double covar_fac,
          const double ftol_rel) :
    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    X(X.begin(), X.n_rows, X.n_cols, false),
    fixed_terms(fixed_terms.begin(), fixed_terms.n_rows,
                fixed_terms.n_cols, false),
    I_len(Rcpp::as<std::vector<double> >(risk_obj["I_len"])), F(F),
    state_dim(a_0.size()),

    n_threads((n_threads > 0) ? n_threads : std::thread::hardware_concurrency()),

    tstart(tstart),
    tstop(tstop),
    is_event_in_bin(is_event_in_bin),
    min_start(Rcpp::as<double>(risk_obj["min_start"])),

    fixed_params(fixed_params),
    fixed_effects(
      (fixed_terms.n_elem > 0) ?
      fixed_terms.t() * fixed_params :
      arma::vec(X.n_cols, arma::fill::zeros)),

    N_fw_n_bw(N_fw_n_bw), N_smooth(N_smooth), N_smooth_final(N_smooth_final),
    forward_backward_ESS_threshold(
      forward_backward_ESS_threshold.isNotNull() ?
        Rcpp::as<Rcpp::NumericVector>(forward_backward_ESS_threshold)[0] :
        N_fw_n_bw / 2.),

    a_0(a_0),
    debug(debug),
    N_first(N_first), nu(nu),
    work_block_size(500),

    covar_fac(covar_fac), ftol_rel(ftol_rel), Q(Q), Q_0(Q_0),
    xtra_covar(Q_tilde)
    {
#ifdef _OPENMP
      omp_set_num_threads(n_threads);
      omp_set_max_active_levels(1);
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
