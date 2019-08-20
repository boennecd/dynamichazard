## Test environments
* Ubuntu 18.04 LTS
  R version 3.6.1
* Ubuntu 16.04 LTS (on travis-ci with codename: xenial)
  R version 3.6.0
* win-builder (devel and release)
* Local Ubuntu 18.04 with R 3.5.2 and with clang 6.0.0 with ASAN and 
  UBSAN checks
* `rhub::check_for_cran()`

## Resubmission
This is a resubmission. In this version I have:

* Used the `FCLEN` macro in the extern declaration of `dsyr` and `FCONE` 
  when calling `dsyr`.

There is a note about the size on most platforms.
