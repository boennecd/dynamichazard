## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.4.2
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.4.2
* win-builder (devel and release)
* Local Ubuntu 17.04 (64-bit)
* The following rhub platforms:
  debian-gcc-devel
  debian-gcc-patched
  debian-gcc-release
  fedora-clang-devel
  fedora-gcc-devel
  linux-x86_64-rocker-gcc-san
  
I use the build with codename trusty for the GCC 4.8.4 with C++11 support.

## R CMD check results
All platforms have a note about the package size.

## Resubmission
> Thanks, please add more small executable examples in your Rd-files.

I have added small examples to all the functions associated with the main 
estimation method, `ddhazard`.

> We also see that your version 0.4.0 had Additional issues: ...

The error ASAN error does not appear on

* Local Ubuntu 17.04 (64-bit) both with ASAN and UBSAN with 
  gcc version 6.3.0
  clang version 4.0.0
* The linux-x86_64-rocker-gcc-san platform from rhub

Sorry for not writing this in the previous submission.
