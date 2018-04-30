#####
# Script from http://rstudio-pubs-static.s3.amazonaws.com/5411_d0603a8032954b98ac187bc20ec644e6.html

## Load libaries
library(survival)
library(plyr)
library(doBy)

## Set wd
cur_path <- str_extract(getwd(), ".+/dynamichazard")
setwd(paste0(cur_path, "/vignettes/Diagnostics"))

## Load data from online
Rossi <- read.table(file = "Rossi.txt", header = TRUE)
head(Rossi)

## Change to long format
Rossi.long <- reshape(data      = Rossi,
                      varying   = paste0("emp", 1:52),
                      v.names   = "employed",
                      timevar   = "time",
                      idvar     = "id",
                      direction = "long",
                      sep       = "")

## Sort by id and time
Rossi.long <- orderBy( ~  + id + time, data = Rossi.long)

## Drop rows where emp is NA (time after event/censoring)
Rossi.long <- Rossi.long[!is.na(Rossi.long$emp),]

## Create time variables and various forms of exposure variables
Rossi.long <- ddply(.data = Rossi.long,
                    .variables = c("id"),
                    .drop = TRUE, .parallel = TRUE,
                    .fun = function(DF) {

                      ## Start time is the start of each interval
                      DF$start <- c(0, head(DF$time, -1))

                      ## Stop time is the end of each interval
                      DF$stop <- DF$time

                      ## Event indicator for each interval
                      DF$event <- 0
                      ## Use arrest value for the last interval
                      DF[nrow(DF),"event"] <- DF[nrow(DF),"arrest"]

                      ## Initial employment status
                      DF$employed.initial <- DF$employed[1]

                      ## Lagged employment status (employment status from last interval matters)
                      DF$employed.lag1 <- c(rep(NA, 1), head(DF$employed, -1))

                      ## Cumulative number of weeks in employment
                      DF$employed.cumsum <- cumsum(DF$employed)

                      ## Ever employed status
                      DF$employed.ever <- as.numeric(DF$employed.cumsum > 0)

                      ## % of time in employment
                      DF$employed.percent <- with(DF, employed.cumsum / stop)*100

                      ## Return DF
                      DF
                    })

## Save output
Rossi <- Rossi.long
save(Rossi, file = "Rossi.RData")


#####
# whas
# Used in e.g.: http://stats.idre.ucla.edu/sas/seminars/sas-survival/

library(smoothHR)

data(whas500)
save(whas500, file = "whas500.RData")
