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

std::unique_ptr
  <std::chrono::time_point<std::chrono::system_clock>>
    PF_logger::last_message_time;

#ifdef _OPENMP
omp_lock_t PF_logger::lock;
#endif
