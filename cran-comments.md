## Test environments
* Ubuntu 18.04 LTS
  R version 3.5.2
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.5.2
* win-builder (devel and release)
* Local Ubuntu 18.04 with R 3.5.2 and with clang 6.0.0 with ASAN and 
  UBSAN checks
* The following rhub platforms:
  fedora-clang-devel
  fedora-gcc-devel
  debian-gcc-patched
  debian-gcc-devel
  debian-gcc-release
  linux-x86_64-rocker-gcc-san
  
I use the build with codename trusty for the GCC 4.8.4 which has C++11 support.

## R CMD check results
All platforms have a note about the package size except the win-builder with 
the devel version.

I still get this error with the UBSAN check: https://github.com/RcppCore/Rcpp/issues/874#issue-337282814

Though, it does not seem to trigger on CRAN and Dirk Eddelbuettel and Kevin 
Ushey indicate that it may not be a problem.

## Resubmission
This is a resubmission. In this version I have:

* Added vignettes which are currently not included due to an error in 
  `.Rbuildignore`.
* Tried to fix the error on Solaris in the CRAN check results.
