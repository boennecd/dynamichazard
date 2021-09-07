## Test environments
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.1.1
* Ubuntu 20.04 LTS with gcc 10.1.0
  R devel 2021-09-06 r80861 with LTO checks
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.1.1 with valgrind
* Ubuntu 20.04 LTS with gcc 10.1.0
  R devel 2021-09-05 r80859 with ASAN and UBSAN
* Github actions on windows-latest (release), macOS-latest (release), 
  ubuntu-20.04 (release), and ubuntu-20.04 (devel)
* win-builder (devel, oldrelease, and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("fedora-clang-devel", "macos-highsierra-release-cran"))`
  
## R CMD check results
There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.

There is a NOTE about a 404 error for one of the DOIs in the description. This  
is for a Journal of Statistical Software paper that will be published soon. 
The paper is the reason for this update.

There is NOTE that Christoffersen is possibly misspelled. It is not.
