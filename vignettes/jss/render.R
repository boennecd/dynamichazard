devtools::reload("knitr")
rm(list = ls())
setwd(dir = "C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/vignettes/jss")
knitr::knit("dynamichazard.Rnw", output = "jss.tex")

tools::texi2pdf("jss.tex")
tools::texi2pdf("jss.tex")
