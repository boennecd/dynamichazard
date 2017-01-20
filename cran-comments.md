## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R 3.3.2
* Ubuntu 14.04.4 LTS (on travis-ci with codename: trusty)
  R 3.3.2
* The win builder on devel image
  
I used the build with codename trusty for the GCC 4.8.4. It is needed for the c++11 syntax I use

## R CMD check results
There were no errors or warnings. The where a some notes

There is one note on travis-ci:
checking installed package size ... NOTE
  installed size is 14.0Mb
  sub-directories of 1Mb or more:
    doc    3.0Mb
    libs   8.9Mb
    R      2.0Mb

The note from win-builder.r-project.org on the devel image is:
* checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    3.0Mb
    libs   2.8Mb

The note on my own laptop is:
* checking installed package size ... NOTE
  installed size is  7.5Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    2.9Mb
    libs   2.5Mb

## ASAN error
I have fixed the ASAN error that is present in the current 0.1.0 version
