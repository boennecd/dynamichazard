## Test environments
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3  
* Ubuntu 18.04 LTS with gcc 8.3.0
  R devel 2020-11-30 r79529 with LTO checks
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3 with valgrind
* Ubuntu 18.04 LTS with clang-6.0.0 
  R devel 2020-11-30 r79529 with ASAN and UBSAN
* Ubuntu 16.04 LTS (on travis-ci)
  R version 4.0.0
* win-builder (devel and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("solaris-x86-patched", "debian-gcc-devel", "fedora-clang-devel"))`
  
## R CMD check results
The issues with CRAN's checks have been solved and the package has been 
updated to handle the new version of `all.equal`.

There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.
