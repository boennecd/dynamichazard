#include "PF_data.h"
#include <chrono>
#include <ctime>

using tp = std::chrono::time_point<std::chrono::system_clock>;
using tp_pointer = std::unique_ptr<tp>;

PF_logger::PF_logger(const bool log, const unsigned int level):
  log(log), level(level) {}

PF_logger::PF_logger(PF_logger&& other):
  log(other.log), level(other.level)
{
  // does not work with some GCC versions https://stackoverflow.com/a/27152585
  // os = std::move(other.os);

  os << other.os.str();
  other.log = false;
}

PF_logger::~PF_logger(){
  if(log){
#ifdef _OPENMP
    omp_set_lock(&lock);
#endif
    Rcpp::Rcout << os.str() << std::endl;
#ifdef _OPENMP
    omp_unset_lock(&lock);
#endif
  }
}



int PF_logger::prefixbuf::sync() {
  return this->sbuf->pubsync();
}

int PF_logger::prefixbuf::overflow(int c) {
  if (c != std::char_traits<char>::eof()) {
    if (this->need_prefix
          && !this->prefix.empty()
          && this->prefix.size() !=
          (unsigned int)this->sbuf->sputn(
              &this->prefix[0], this->prefix.size())) {
          return std::char_traits<char>::eof();
    }
    this->need_prefix = c == '\n';
  }
  return this->sbuf->sputc(c);
}

PF_logger::prefixbuf::prefixbuf(std::string const &prefix, std::streambuf *sbuf)
  : prefix(prefix), sbuf(sbuf) , need_prefix(true) {}



PF_logger::oprefixstream::oprefixstream(std::string const& prefix, std::ostream& out)
  : prefixbuf(prefix, out.rdbuf())
  , std::ios(static_cast<std::streambuf*>(this))
  , std::ostream(static_cast<std::streambuf*>(this)) {}



std::string PF_logger::get_prefix(const unsigned int level){
  std::stringstream ss;

#ifdef _OPENMP
  unsigned int me;
  if(omp_get_level() > 1 && !omp_get_nested()){
    me = omp_get_ancestor_thread_num(omp_get_level() - 1);

  } else
    me = omp_get_thread_num();
  ss << "Thread:" << std::setw(3) << me  << "\t";
#endif

  ss << "delta T: " << std::setw(10) << std::setprecision(6)
     << get_elapsed_seconds_n_set_last_message_time() << "\t";

  ss << std::string(n_spaces * (level - 1), ' ');
  return ss.str();
}

class tracker {
public:
  template<class T>
  tracker(T val): time_pt(val) { }

  std::chrono::time_point<std::chrono::system_clock> time_pt;
};

std::unique_ptr<tracker> PF_logger::last_message_time;

#ifdef _OPENMP
omp_lock_t PF_logger::lock;
#endif

double PF_logger::get_elapsed_seconds_n_set_last_message_time(){
  double elapsed_seconds;

#ifdef _OPENMP
  omp_set_lock(&lock);
#endif
  tp_pointer now(new tp());
  *now = std::chrono::system_clock::now();

  if(last_message_time){
    std::chrono::duration<double> tmp = *now - last_message_time->time_pt;
    elapsed_seconds = tmp.count();
  } else {
    elapsed_seconds = std::numeric_limits<double>::quiet_NaN();

  }

  last_message_time.reset(new tracker(*now.release()));

#ifdef _OPENMP
  omp_unset_lock(&lock);
#endif

  return elapsed_seconds;
}
