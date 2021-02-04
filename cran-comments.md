## Test environments
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3
* Ubuntu 18.04 LTS with clang-6.0.0 
  R devel 2020-12-15 r79637 with ASAN and UBSAN
* Ubuntu 16.04 LTS (on travis-ci)
  R version 4.0.0
* win-builder (devel and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("debian-gcc-devel", "debian-clang-devel", "fedora-clang-devel", "fedora-clang-devel"))`
  
## R CMD check results
The test which caused the package to be removed has been fixed.

There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.
