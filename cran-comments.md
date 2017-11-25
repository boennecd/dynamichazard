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

The ASAN error does not appear on

* Local Ubuntu 17.04 (64-bit) both with ASAN and UBSAN with 
  gcc version 6.3.0
  clang version 4.0
* The linux-x86_64-rocker-gcc-san platform from rhub

Sorry for not writing this in the previous submission.

## Resubmission
I see the errors on 
https://www.stats.ox.ac.uk/pub/bdr/memtests/clang-ASAN/dynamichazard/00check.log
https://www.stats.ox.ac.uk/pub/bdr/memtests/valgrind/dynamichazard/dynamichazard-Ex.Rout

I cannot reproduce them on
* Local Ubuntu 17.04 (64-bit) both with ASAN and UBSAN with
  clang version 4.0
* The linux-x86_64-rocker-gcc-san platform from rhub

Though, I am quite sure I have fixed the error in this version.

I also see the error in
https://www.stats.ox.ac.uk/pub/bdr/memtests/gcc-ASAN/dynamichazard/00check.log

Again, I cannot re-produce it. Here I cannot figure out what is wrong. I might 
have fixed it in this version.

### Longer comment to the gcc log
Please, skip this if you have no interest. I added this comment as there seems
to be similar issues in other packages.

There seems to be similar issues in the `eDMA` package and the bug exactly  
match with the one in the `Morpho` package. The log from the `Morpho` package 
is partly given below

  ==43264==ERROR: AddressSanitizer: stack-use-after-scope on address 0x7fff6b2ce230 at pc 0x7f8e99e90ed5 bp 0x7fff6b2cde00 sp 0x7fff6b2cddf0
  WRITE of size 8 at 0x7fff6b2ce230 thread T0
      #0 0x7f8e99e90ed4 in createL._omp_fn.0 /data/gannet/ripley/R/test-3.5/RcppArmadillo/include/armadillo_bits/subview_meat.hpp:37
      #1 0x7f8eae8b8cde in GOMP_parallel (/lib64/libgomp.so.1+0xdcde)
      #2 0x7f8e99e92cea in createL /data/gannet/ripley/R/packages/tests-gcc-SAN/Morpho/src/createL.cpp:14
      ...
      
What seems to be common is that we both construct a subview inside a parallel 
region which uses the constructor in `subview_meat.hpp:37`. I do it in 
`GMA_solver.cpp` at

  sym_mat_rank_one_update(h_2d_neg, X_t.col(i), my_X_cross);
  
  
They do it in `Morpho` in `createL.cpp` at 

  mat diff = MatrixA.row(i)-MatrixA.row(j);
  ...
  K(i,j) = -sqrt(dot(diff,diff));
  ...
  K(i,j) = r2*log(r2);
  
The new version I have submitted does not create a subview but uses the 
`unsafe_col` method instead. I am not sure if this helps but now there are no 
calls to `subview` constructor.
