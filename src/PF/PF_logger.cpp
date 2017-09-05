#include "PF_data.h"

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

    os << std::endl;
    Rcpp::Rcout << os.str();

#ifdef _OPENMP
    omp_unset_lock(&lock);
#endif
  }
}

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

std::unique_ptr
  <std::chrono::time_point<std::chrono::system_clock>>
    PF_logger::last_message_time;

#ifdef _OPENMP
omp_lock_t PF_logger::lock;
#endif
