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
    "Have you cleaned up the TODO file?",
    "Have you run non-cran tests with ubsan and asan?",
    "Is 'DDHAZ_DEBUG' defined?"
  )
}
