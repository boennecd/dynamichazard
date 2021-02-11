## Test environments
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3
* Ubuntu 18.04 LTS with clang-6.0.0 
  R devel 2021-02-03 r79933 with ASAN and UBSAN
* Ubuntu 16.04 LTS (on travis-ci)
  R version 4.0.0
* win-builder (devel and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("debian-gcc-devel", "debian-clang-devel", "fedora-clang-devel", "fedora-clang-devel"))`
  
## R CMD check results
The issue with test which caused the package to be removed form CRAN has been 
fixed.

There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.

## Resubmission
This is a resubmission. I have addressed the issues stated below.

> Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or
> 'F' as vector names.)

`T` and `F` is no longer used. Instead `TRUE` and `FALSE` is used. Also, T and F 
are no longer used for vector names.

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      ddFixed.Rd: \value
      ddhazard_app.Rd: \value
      logLik.ddhazard.Rd: \value
      plot.ddhazard_space_errors.Rd: \value
      plot.ddhazard.Rd: \value
      predict.ddhazard.Rd: \value
      print.ddhazard_boot.Rd: \value
      residuals.ddhazard.Rd: \value

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary in every case.

Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.


You are setting options(warn=-1) in your function. This is not allowed.
Please rather use suppressWarnings() if really needed.

Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...

e.g.: inst/doc.. files


Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.: ddsurvcurve.Rd
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

old <- options(digits = 3)
...
options(old)

Please do not modify the .GlobalEnv. This is not allowed by the CRAN
policies.

Please do not modify the global environment (e.g. by using <<-) in your
functions. This is not allowed by the CRAN policies.

Please ensure that you do not use more than 2 cores in your examples,
vignettes, etc.


Additionally:
Have the issues why your package was archived been fixed?
Please explain this in the submission comments.

