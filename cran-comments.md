## New version submission

CRAN archived `dynamichazard` on 2023-09-01 after the `ddhazard` vignette could not rebuild without `captioner`. Version 1.0.3 uses a LaTeX figure reference in that vignette and removes the dependency.

## Test environment

- Debian GNU/Linux 13 (x86_64), R 4.6.2, GCC and GNU Fortran 14.2.0, with OpenMP enabled.
- Built with `R CMD build --compact-vignettes=gs+qpdf .`.
- Checked with `R CMD check --as-cran dynamichazard_1.0.3.tar.gz`.

## R CMD check results

0 ERRORs, 0 WARNINGs, and 1 NOTE. The NOTE identifies the package as a new submission after archival. The examples, tests, vignette rebuilds, and PDF and HTML manuals passed.

The installed size was 51.4 MB. The shared library accounts for 45.2 MB, primarily from compiler debug information.

The current CRAN source package index lists no reverse dependencies.
