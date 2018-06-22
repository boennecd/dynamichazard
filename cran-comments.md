## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.5.0
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.5.0
* win-builder (devel and release)
* Local Ubuntu 17.04 (64-bit)
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
All platforms have a note about the package size.
