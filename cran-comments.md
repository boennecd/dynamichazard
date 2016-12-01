## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.3.2
* Ubuntu 14.04.4 LTS (on travis-ci with codename: trusty)
  R 3.3.2
  
I used the build with codename trusty for the GCC 4.8.4. It is needed for the c++ 11 syntax I use

## R CMD check results
There were no ERRORs or WARNINGs. 

There are two notes travis-ci:
checking installed package size ... NOTE
installed size is 11.0Mb
  sub-directories of 1Mb or more:
    libs   8.0Mb
    R      2.0Mb
checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences

I do not get the size note on my Windows computer. The resulting package is only 4.3MB