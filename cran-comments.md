## Test environments
* Ubuntu 18.04 LTS
  R version 3.5.3
* Ubuntu 16.04.5 LTS (on travis-ci with codename: xenial)
  R version 3.5.2
* win-builder (devel and release)
* Local Ubuntu 18.04 with R 3.5.2 and with clang 6.0.0 with ASAN and 
  UBSAN checks
* The following rhub platforms:
  fedora-clang-devel
  fedora-gcc-devel
  debian-gcc-devel
  linux-x86_64-rocker-gcc-san

## R CMD check results
There is a note about the size on most platforms.

I cannot reproduce the errors on devel on CRAN related to `dtrtri`. I have 
tried with R devel from 2019-04-04 (r76316) on Ubuntu 18.04 with gcc 
7.3.0.

I could not reproduce the `free(): double free detected in tcache` on 
r-devel-linux-x86_64-fedora-clang but I did find some issues when using 
clang 8.0.0 and Valgrind. These are fixed.
