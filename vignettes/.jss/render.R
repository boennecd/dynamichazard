dev.off()
rm(list = ls())

options(ddhazard_max_threads = 7)

setwd(dir = paste0(
  stringr::str_extract(getwd(), "^.+/dynamichazard"), "/vignettes/.jss"))
knitr::knit("dynamichazard.Rnw", output = "jss.tex")

tools::texi2pdf("jss.tex")
tools::texi2pdf("jss.tex")


knitr::purl("dynamichazard.Rnw", output = "reproducible_script.R")
