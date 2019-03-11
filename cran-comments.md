## Test environments
* Ubuntu 18.04 LTS
  R version 3.5.2
* Ubuntu 16.04.5 LTS (on travis-ci with codename: xenial)
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

## R CMD check results
All platforms have a note about the package size except the win-builder with 
the devel version.

I still get this error with the UBSAN check: https://github.com/RcppCore/Rcpp/issues/874#issue-337282814

I have added `suppressWarnings(RNGversion("3.5.0"))` when I run the tests. 

The C++ code has been more or less been rewritten for the particle filters and
smoothers which may have fixed the error on r-patched-solaris-x86 which I 
cannot reproduce.

I have used Ghostscript this time so the tarball is less than 5 MB.
