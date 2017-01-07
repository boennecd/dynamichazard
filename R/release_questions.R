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
    "Did you set args = c('--compact-vignettes=gs+qpdf')"
  )
}
