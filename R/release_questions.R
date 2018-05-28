release_questions <- function() {
  c(
    "Have you updated the 'Sim_study_with_logit' vignette?",
    "Have you updated the 'Bootstrap_illustration' vignette?",
    "Have you pushed to github before running release?",
    "Have you searched the project through for 'TODO'?",
    "Have you updated cran-comments.md?",
    "Have you updated NEWS.md?",
    "Have you run devtools::spell_check()?",
    "Have you checked reverse dependencies?",
    "Have you set args = c('--compact-vignettes=gs+qpdf')",
    "Have you search through for CHECK and TODO?",
    "Have you cleaned up the TODO file?",
    "Are 'ARMA_NO_DEBUG' and 'ARMA_DONT_PRINT_ERRORS' defined in arma_n_rcpp.h'?",
    "Have you run non-cran tests with ubsan and asan?"
  )
}
