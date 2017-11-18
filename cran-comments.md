## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.4.2
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.4.2
* win-builder (devel and release)
* The following rhub platforms:
  debian-gcc-devel
  debian-gcc-patched
  debian-gcc-release
  fedora-clang-devel
  fedora-gcc-devel
  
I use the build with codename trusty for the GCC 4.8.4 with C++11 support.

## R CMD check results
All platforms have a note about the package size.
