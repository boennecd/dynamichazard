## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R 3.4.0
* Ubuntu 14.04.4 LTS (on travis-ci with codename: trusty)
  R 3.4.0
* The win builder on release image and develop image
  
I use the build with codename trusty for the GCC 4.8.4 C++11

## R CMD check results
There were no errors or warnings. The are some notes all regarding the size

There is one note on travis-ci:
* checking installed package size ... NOTE
  installed size is 18.2Mb
  sub-directories of 1Mb or more:
    doc    3.4Mb
    libs  12.7Mb
    R      2.0Mb

The note from win-builder.r-project.org on the release image is:
* checking installed package size ... NOTE
  installed size is  8.3Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    3.4Mb
    libs   2.9Mb

The note from win-builder.r-project.org on the develop image is:
* installed size is  8.7Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    3.4Mb
    libs   3.2Mb

The note on my own laptop is:
* checking installed package size ... NOTE
  installed size is  8.4Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    3.4Mb
    libs   2.9Mb

## Resubmission
Should have fixed build errors on: 

* r-devel-linux-x86_64-fedora-clang
* r-patched-solaris-sparc
* r-patched-solaris-x86

## Resubmission
Changed the rule Makevars and Makevars.win that professor Kurt Hornik and professor Brian Ripley mentioned.

Changed dynamichazard/src/dchur.f to take in LOGICAL arguments instead of CHARACTER arguments.

Sorry for violating the CRAN policies with the rule and wasting your time. Sorry for sending a HTML formatted e-mail to professor Brian Ripley and for not using the mailing list.

## Resubmission
Removed the __thread to make the package build on r-patched-solaris-sparc.
Reduced the size of the vignettes.
Minor small changes to decrease the computation speed.

## Resubmission
Fixed issues with r-patched-solaris-x86. I figure it was due to the $(SHLIB_OPENMP_CXXFLAGS) flag in last submission which triggered armadillo to use OpenMP with C++11 code which caused the error. 
Decreased the size of the R code in the final package.
Fixed the clang-UBSAN errors in https://www.stats.ox.ac.uk/pub/bdr/memtests/clang-UBSAN/dynamichazard/tests/testthat.Rout. I tested the new version on Ubuntu 17.04 with clang-3.9. However, I got a lot of errors from Rcpp/include/Rcpp/routines.h as in https://www.stats.ox.ac.uk/pub/bdr/memtests/clang-UBSAN/Rcpp/build_vignettes.log (checked on 24/06/2017 at 2 pm). I gather these can be disregared from Kurt Hornik's mail on 	26. June 2017 at 14.53. 
