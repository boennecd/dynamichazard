## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.5.0
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.5.1
* win-builder (devel and release)
* Local Ubuntu 17.04 (64-bit) with R devel and with clang 6.0.0 with ASAN and 
  UBSAN checks
* The following rhub platforms:
  debian-gcc-devel
  fedora-clang-devel
  fedora-gcc-devel
  debian-gcc-patched
  debian-gcc-release
  linux-x86_64-rocker-gcc-san
  
I use the build with codename trusty for the GCC 4.8.4 which has C++11 support.

## R CMD check results
All platforms have a note about the package size except the win-builder with 
the devel version.

I still get this error with the UBSAN check: https://github.com/RcppCore/Rcpp/issues/874#issue-337282814

Though, it does not seem to trigger on CRAN and Dirk Eddelbuettel and Kevin 
Ushey indicate that it may not be a problem.
