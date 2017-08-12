#include "PF_data.h"

PF_logger::PF_logger(const bool log, const unsigned int level):
  log(log), level(level) {}

PF_logger::PF_logger(PF_logger&& other)
{
  // does not work with some GCC versions https://stackoverflow.com/a/27152585
  // os = std::move(other.os);
  os << other.os.str();
  log = other.log;
  level = other.level;
  other.log = false;
}

PF_logger::~PF_logger(){
  if(log){
    std::lock_guard<std::mutex> guard(mx);
    os << std::endl;
    Rcpp::Rcout << os.str();
  }
}

std::mutex PF_logger::mx;
