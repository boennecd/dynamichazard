#include "PF_data.h"

PF_logger::PF_logger(const bool log, const unsigned int level):
  log(log), level(level),
  os_w_prefix(get_prefix(level), os) {}

PF_logger::PF_logger(PF_logger&& other):
  log(other.log), level(other.level),
  os_w_prefix(get_prefix(level), os)
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

#ifdef _OPENMP
omp_lock_t PF_logger::lock;
#endif
