## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.5.0
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.5.0
* win-builder (devel and release)
* Local Ubuntu 17.04 (64-bit) with R devel and with clang 4.0.0 with ASAN and 
  UBSAN checks
* The following rhub platforms:
  debian-gcc-release
  fedora-clang-devel
  fedora-gcc-devel
  
  I failed to build with these platforms. It seems to be a bug with rhub https://github.com/r-hub/rhub/issues/141
  debian-gcc-devel
  debian-gcc-patched
  linux-x86_64-rocker-gcc-san
  
I use the build with codename trusty for the GCC 4.8.4 which has C++11 support.

## R CMD check results
All platforms have a note about the package size expect the win-builder with the 
devel version.

## Resubmission
This is a resubmission. In this version I have:

 * Skipped some of the tests on r-patched-linux-x86_64 and 
   r-release-linux-x86_64. The test failed due to errors in the biglm package 
   version 0.9-1. I could not reproduce the errors but now the tests do not fail 
   if the calls to `biglm::bigglm` fails.
 * Fixed the UBSAN error. I could re-produce it after I upgraded to 
   clang-6.0.
