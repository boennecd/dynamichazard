## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R 3.3.2
* Ubuntu 14.04.4 LTS (on travis-ci with codename: trusty)
  R 3.3.2
* The win builder on both devel and release image
  
I used the build with codename trusty for the GCC 4.8.4. It is needed for the c++11 syntax I use

## R CMD check results
There were no errors and only one warning which is on travis-ci. 

There are two notes travis-ci:
checking installed package size ... NOTE
  installed size is 11.4Mb
  sub-directories of 1Mb or more:
    doc    1.0Mb
    libs   8.3Mb
    R      2.0Mb
checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

The warning on travis-ci is:
checking sizes of PDF files under ‘inst/doc’ ... WARNING
  ‘gs+qpdf’ made some significant size reductions:
     compacted ‘ddhazard.pdf’ from 509Kb to 207Kb

The notes from win-builder.r-project.org on the devel image are:
* checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    1.0Mb
    libs   2.8Mb
* checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

The notes from win-builder.r-project.org on the release image are:
* checking installed package size ... NOTE
  installed size is  5.5Mb
  sub-directories of 1Mb or more:
    R      2.0Mb
    doc    1.0Mb
    libs   2.4Mb
* checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

I am not sure how to deal with pdf size warning on travis-ci. Would you recommend that a pre-build the vignette and add the pdfs? 

Further, I cannot get rid of the note in regard to the description
