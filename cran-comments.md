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
* checking installed package size ... NOTE
  installed size is 11.4Mb
  sub-directories of 1Mb or more:
    doc    1.0Mb
    libs   8.3Mb
    R      2.0Mb
* checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

The warning on travis-ci is:
* checking sizes of PDF files under ‘inst/doc’ ... WARNING
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

Further, I cannot get rid of the note regarding the description

## Resubmission
This is a resubmission. In this version, I have:

* Updated the Description such that I got rid of the previous note. Further, I have cited the papers that the package is based on. I added period at the end and corrected typos. Though, the description reads 'coxph()' where it read 'coxph' before. It refers to the survival package function. Is this the right way to do it?   

* Got rid of the compiler warnings that professor Hornik mentioned. Though, I could not get my compiler on Windows or Travis-ci to yield the warnings. What flags do you use to build with? I figure you should not get the warnings now as I have removed the tab characters and added END DO statements

* Build package on my site with --compact-vignettes=gs+qpdf

* Sorry for the issues with the authors pointed out by Dr. Uwe Ligges. I am not sure how to cite the source code I am using by other authors. The code from other authors are: 
  The code in src/biglm/boundedQRf.f. It was included as in the biglm package. That is, no citation in of it in the DESCRIPTION file. Is it correct that I have added Allan Miller (the author of the Fortran code) with role = c("ctb")?
  The code in src/thread_pool.cpp is a modified version of the code from the book "Williams, Anthony. C++ concurrency in action. London, 2012." The code is available from: https://manning-content.s3.amazonaws.com/download/0/78f6c43-a41b-4eb0-82f2-44c24eba51ad/CCiA_SourceCode.zip. It requires no login and I have not found any part of the book where Williams states if/how one can use the code. The code is slightly modified. I have copied the License_1_0.txt from the .zip file and included it in top of the code. Is it correct to add Anthony Williams with role = c("ctb")? Further, I have added a "Boost developers" with role = c("ctb", "cph") 
  src/stats/lm.c and src/stats/statsR.h is from the source code of R. The code is from R-3.3.2/src/library/stats/src/lm.c and R-3.3.2/src/library/stats/src/statsR.h. I have added person("R-core", email = "R-core@R-project.org", role = c("ctb", "cph"))
  I am sorry for the inconvenience and wasting your time
  
* Added DOI in the description
