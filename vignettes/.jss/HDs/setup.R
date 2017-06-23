#####
# This script assumes that you are in the directory .../dynamichazard or
# .../dynamichazard/... where the exists a directory called
# .../dynamichazard/vignettes/jss/HDs/

library(RSQLite)
library(sqldf)
library(tcltk)
require(stringr)
library(testthat)
library(survival)
library(testthat)

#####
# Define object used through out
cur_ls <- ls()
tmp_files <- list()
new_folder <- NULL
pb <- NULL
conn <- NULL
log_file <- NULL
warn_old <- getOption("warn")
db <- "driver_stats.db"
final_tbl_name <- "HDs.RDS"

sqlFromFile <- function(file){
  sql <- readLines(file)
  sql <- gsub("(--.*$)|(\\t)"," " ,sql)
  sql <- unlist(str_split(paste(sql,collapse=" "),";"))
  sql <- sql[grep("^ *$", sql, invert=T)]
  sql
}

dbSendQueries <- function(con,sql){
  dummyfunction <- function(sql,con){
    cat("Running the following sql:\n", sql, sep = "\n")
    tm <- system.time(out <- dbSendQuery(con,sql))
    cat("Run time info is:\n")
    print(tm)
    out
  }
  lapply(sql, dummyfunction, con)
}

reset_dir <- bquote(setwd(.(getwd())))

set_dir <- bquote(
  setwd(.(paste0(
    stringr::str_match(getwd(), ".+dynamichazard"), "/vignettes/.jss/HDs/"))))

conn_exp <- expression({
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, db, cache_size = 8e5, synchronous = "off")
})

tryCatch({
  eval(set_dir)
  options(warn = 1)

  # Check if db should made
  make_db <- file.exists(db)
  if(make_db)
    make_db <- "YES" ==
    winDialog("yesno", "Do you want to remake the data base for the hard disks?")

  log_file <- file("setup.log", if(make_db) "wt" else "at")
  sink(log_file, append = F, split = T, type = "output")
  sink(log_file, append = T, split = F, type = "message")
  cat("\n#############\n ",  as.character(as.POSIXlt(Sys.time())), "\n#############\n\n")

  #####
  # Create database
  if(make_db){
    if(file.exists(db))
      unlink(db, force = T)

    eval(conn_exp)

    queries <- sqlFromFile("create_table.sql")
    dbSendQueries(conn, queries)

    dbListTables(conn)

    suffixes <- c("2013", "2014", "2015", "Q1_2016", "Q2_2016", "Q3_2016")
    for(suffix in suffixes){
      data_url <- paste0("https://f001.backblaze.com/file/Backblaze-Hard-Drive-Data/data_", suffix, ".zip")

      # Dl data files and unzip
      cat("Downloading files for", suffix, "\n")
      tmp_files$data <- tempfile()
      download.file(data_url, tmp_files$data)
      tmp_files$data_dir <- suffix
      if(substr(suffix, 0, 1) == "Q")
        tmp_files$data_dir <- paste0("data_", tmp_files$data_dir)
      unzip(tmp_files$data)
      unlink(tmp_files$data)
      tmp_files$data <- NULL

      # Import
      cat("Adding data to database for", suffix, "\n")
      fs <- dir(tmp_files$data_dir)
      pb <- txtProgressBar(
        1, length(fs), 0, title = paste("Progress for file suffix", suffix),
        style=3)
      for(i in seq_along(fs)){
        setTxtProgressBar(pb, i)
        f <- fs[i]
        f <- paste0("./", tmp_files$data_dir, "/", f)
        read.csv.sql(
          f, paste(
            "INSERT INTO drive_stats",
            # "SELECT date, serial_number, model, capacity_bytes, failure, smart_1_normalized, smart_1_raw, smart_2_normalized, smart_2_raw, smart_3_normalized, smart_3_raw, smart_4_normalized, smart_4_raw, smart_5_normalized, smart_5_raw, smart_7_normalized, smart_7_raw, smart_8_normalized, smart_8_raw, smart_9_normalized, smart_9_raw, smart_10_normalized, smart_10_raw, smart_11_normalized, smart_11_raw, smart_12_normalized, smart_12_raw, smart_13_normalized, smart_13_raw, smart_15_normalized, smart_15_raw, smart_183_normalized, smart_183_raw, smart_184_normalized, smart_184_raw, smart_187_normalized, smart_187_raw, smart_188_normalized, smart_188_raw, smart_189_normalized, smart_189_raw, smart_190_normalized, smart_190_raw, smart_191_normalized, smart_191_raw, smart_192_normalized, smart_192_raw, smart_193_normalized, smart_193_raw, smart_194_normalized, smart_194_raw, smart_195_normalized, smart_195_raw, smart_196_normalized, smart_196_raw, smart_197_normalized, smart_197_raw, smart_198_normalized, smart_198_raw, smart_199_normalized, smart_199_raw, smart_200_normalized, smart_200_raw, smart_201_normalized, smart_201_raw, smart_223_normalized, smart_223_raw, smart_225_normalized, smart_225_raw, smart_240_normalized, smart_240_raw, smart_241_normalized, smart_241_raw, smart_242_normalized, smart_242_raw, smart_250_normalized, smart_250_raw, smart_251_normalized, smart_251_raw, smart_252_normalized, smart_252_raw, smart_254_normalized, smart_254_raw, smart_255_normalized, smart_255_raw",
            "SELECT date, serial_number, model, capacity_bytes, failure, smart_1_raw, smart_2_raw, smart_3_raw, smart_4_raw, smart_5_raw, smart_7_raw, smart_8_raw, smart_9_raw, smart_10_raw, smart_11_raw, smart_12_raw, smart_13_raw, smart_15_raw, smart_183_raw, smart_184_raw, smart_187_raw, smart_188_raw, smart_189_raw, smart_190_raw, smart_191_raw, smart_192_raw, smart_193_raw, smart_194_raw, smart_195_raw, smart_196_raw, smart_197_raw, smart_198_raw, smart_199_raw, smart_200_raw, smart_201_raw, smart_223_raw, smart_225_raw, smart_240_raw, smart_241_raw, smart_242_raw, smart_250_raw, smart_251_raw, smart_252_raw, smart_254_raw, smart_255_raw",
            "FROM file"),
          dbname = db)
      }
      close(pb)
      unlink(tmp_files$data_dir, recursive = T, force = T)
    }

    #####
    # Make survival table
    cat("Making surival table in data base\n")

    queries <- sqlFromFile("create_surv_table.sql")
    dbSendQueries(conn, queries)

  } else
    cat("Db exists and will not be made again\n")

  #####
  # Make data frame for analysis and save as .RDS
  if(is.null(conn))
    eval(conn_exp)

  # surv_tbl <- readRDS("tmp.RDS")
  surv_tbl <- dbReadTable(conn, "drive_stats_survival")

  # Make stats functions and use it
  stats_func <- function(){
    cat("The total number of rows are", nrow(surv_tbl), "with ",
        length(unique(surv_tbl$serial_number)), "individuals",
        "\n")
    cat("Str on the table yields:\n")
    str(surv_tbl)
    cat("Summary stats are:\n")
    print(summary(surv_tbl))
  }
  stats_func()

  # Format data
  cat("Formatting data columns\n")
  for(n in c("date", "min_date", "max_date")){
    surv_tbl[[n]] <- as.Date(surv_tbl[[n]])
  }

  cat("There are a few different model where some are quite sparse:\n")
  sort(print(xtabs(~ model, surv_tbl, subset = !duplicated(surv_tbl$serial_number))))

  test_that("All disks only have one model", expect_true(
    all(tapply(surv_tbl$model, surv_tbl$serial_number, function(x) all(x[1] == x)))
  ))

  cat("Thus, we focus on the manufacturer:\n")
  surv_tbl$manufacturer <- str_extract(surv_tbl$model, "^([A-z]+(?=\ ))|(ST)")

  # "As noted in the BackBlaze article, models listed as HGST or Hitachi should be under the same label"
  # Source: http://blog.applied.ai/survival-analysis-part2/

  surv_tbl$manufacturer[surv_tbl$manufacturer == "Hitachi"] <- 'HGST'

  test_that("We found a manufacturer for all observations",
            expect_true(!any(is.na(surv_tbl$manufacturer))))
  tmp_tbl <- data.frame(ftable(xtabs(
    ~ manufacturer + model, surv_tbl, subset = !duplicated(surv_tbl$serial_number)),
    row.vars = c("manufacturer", "model")))
  tmp_tbl <- tmp_tbl[tmp_tbl$Freq != 0, ]
  print(tmp_tbl)
  print(xtabs(~ manufacturer, surv_tbl))

  cat("We remove SAMSUNG and TOSHIBA\n")
  surv_tbl <- surv_tbl[!surv_tbl$manufacturer %in% c("SAMSUNG", "TOSHIBA"), ]
  surv_tbl$model <- as.factor(surv_tbl$model)
  surv_tbl$manufacturer <- as.factor(surv_tbl$manufacturer)
  stats_func()

  cat("A few HDs have more than one size. These are:\n")

  all_sizes_equal_exp <- expression(
    tapply(surv_tbl$capacity_bytes_GB, surv_tbl$serial_number,
           function(x) all(x[1] == x)))
  has_more_sizes <- which(!eval(all_sizes_equal_exp))
  test_that("There are some HDs with more than onse size",
            expect_gt(length(has_more_sizes), 0))

  print(
    surv_tbl[surv_tbl$serial_number %in% names(has_more_sizes),
             c("serial_number", "capacity_bytes_GB")])

  cat("We handle these on a case by case basis\n")
  replacement_sizes <- c(
    "PL2331LAGPJ89J" = 4001,
    "PL2331LAGSUPTJ" = 4001,
    "PL2331LAGSU5RJ" = 4001,
    "W300935A" = 4001,
    "W300A7LV" = 4001,
    "W300A8MC" = 4001,
    "W300A9M2" = 4001,
    "W300BV4X" = 4001,
    "W300C6K2" = 4001,
    "W300CDJ1" = 4001,
    "W300CSMW" = 4001,
    "W300CSYX" = 4001,
    "W300F1LC" = 4001,
    "WD-WCAU45478169" = 1000,
    "WD-WCAU47760976" = 1000,
    "Z300JZGN" = 4001,
    "Z300K0NG" = 4001,
    "Z300K1HK" = 4001,
    "Z300K1L6" = 4001,
    "Z300K204" = 4001,
    "Z300K8MM" = 4001,
    "Z300K9LT" = 4001,
    "Z300KCGQ" = 4001,
    "Z300KCR8" = 4001,
    "Z300KCT4" = 4001,
    "Z300KMMG" = 4001,
    "Z300KN1Y" = 4001,
    "Z305GFM1" = 4001)
  for(i in seq_along(replacement_sizes)){
    indx <- which(surv_tbl$serial_number == names(replacement_sizes)[i])
    stopifnot(any(surv_tbl$capacity_bytes_GB[indx] == replacement_sizes[i]))
    surv_tbl$capacity_bytes_GB[indx] <- replacement_sizes[i]
  }

  test_that("There are not HDs with more than onse size",
            expect_true(all(eval(all_sizes_equal_exp))))

  cat("There are some very small HDS\n")
  xtabs(~surv_tbl$capacity_bytes_GB)

  cat("We remove the smallest\n")
  surv_tbl <- surv_tbl[surv_tbl$capacity_bytes_GB > 1e2, ]
  stats_func()

  cat("Some HDs fails more than once")
  xtabs(~ surv_tbl$n_fails, subset = !duplicated(surv_tbl$serial_number))

  cat("We remove these")
  surv_tbl <- surv_tbl[surv_tbl$n_fails < 2, ]
  stats_func()

  cat("We format the rest of the columns")
  surv_tbl$serial_number <- as.factor(surv_tbl$serial_number)
  for(n in c("smart_10_raw", "smart_12_raw",
             "smart_187_raw", "smart_188_raw", "smart_189_raw",
             "smart_196_raw", "smart_198_raw", "smart_201_raw")){
    is_missing <- surv_tbl[[n]] == ""
    surv_tbl[[n]] <- as.integer(surv_tbl[[n]])
    surv_tbl[[n]][is_missing] <- NA_integer_
  }

  test_that("We miss all values of 201", expect_true(
    all(is.na(surv_tbl$smart_201_raw))))
  surv_tbl$smart_201_raw <- NULL

  equal_min_n_max_h <- expression(
    sum(surv_tbl$min_hours == surv_tbl$max_hours))
  cat("There are", eval(equal_min_n_max_h), "HDs with equal min and max hours. We remove these\n")

  test_that("There are some cases with equal min and max hours",
            expect_gt(eval(equal_min_n_max_h), 0))
  surv_tbl <- surv_tbl[surv_tbl$min_hours < surv_tbl$max_hours, ]
  test_that("There are no longer cases with equal min and max hours",
            expect_equal(eval(equal_min_n_max_h), 0))

  cat("There are some HDs which seems to have run for more than 20 years. We remove these\n")

  test_that("Some have run for more than 20",
            expect_true(any(surv_tbl$max_hours > 365 * 24 * 20)))

  surv_tbl <- surv_tbl[surv_tbl$max_hours < 365 * 24 * 20, ]

  # par(mfcol = c(2, 1))
  # by(surv_tbl[!duplicated(surv_tbl$serial_number), ],
  #    surv_tbl$n_fails[!duplicated(surv_tbl$serial_number)],
  #    function(x) hist(x$max_hours))

  stats_func()

  # Check the smart_12
  surv_tbl <- surv_tbl[order(surv_tbl$serial_number, surv_tbl$smart_9_raw), ]
  test_that("smart_12 is sorted", expect_true(
    all(!unlist(tapply(surv_tbl$smart_12_raw, surv_tbl$smart_12_raw, is.unsorted)))))

  #####
  # Reproduce results from http://blog.applied.ai/survival-analysis-part-4/
  tmp_dat <- surv_tbl[!duplicated(surv_tbl$serial_number), ]
  tmp_dat$size <- round(tmp_dat$capacity_bytes_GB / 1e3, 1)
  tmp_dat <- tmp_dat[tmp_dat$size < 7 & tmp_dat$size > 1.4, ]

  tmp_dat$manufacturer <- relevel(
    tmp_dat$manufacturer, ref = "HGST")
  tmp_dat$size <- as.factor(tmp_dat$size)
  tmp_dat$size <- relevel(tmp_dat$size, ref = 1.5)

  tmp_fit <- coxph(
    Surv(max_hours, n_fails > 0) ~ manufacturer + size, tmp_dat)

  # tmp_fit$coefficients

  test_that("Roughly matches w/ 2015 example", {
    coefs <- tmp_fit$coefficients
    expect_true(2.5 < coefs["manufacturerST"] && coefs["manufacturerST"] < 3.5)
    expect_true(2 < coefs["manufacturerWDC"] && coefs["manufacturerWDC"] < 2.5)
    expect_true(1 < coefs["size2"] && coefs["size2"] < 1.3)
    expect_true(1.5 < coefs["size3"] && coefs["size3"] < 3)
  })

  # # Compare with size with http://blog.applied.ai/survival-analysis-part2/
  # # It seems like they have purchased some new HDs of size 4Tb and larger since
  # # this analysis
  # xtabs(~tmp_dat$size)

  #####
  # Make table for survival analysis

  # Need to sort first
  surv_tbl <- surv_tbl[order(surv_tbl$serial_number, surv_tbl$smart_9_raw), ]
  surv_tbl$size_tb <- round(surv_tbl$capacity_bytes_GB / 1e3, 1)

  # Make base frame
  tmp <- subset(
    surv_tbl, !duplicated(surv_tbl$serial_number),
    select = c(serial_number, model, manufacturer, n_fails, n_records, size_tb,
               min_date, max_date, min_hours, max_hours))
  final_tbl <- tmerge(
    tmp, tmp, id = serial_number, fails = event(max_hours, n_fails),
    tstart = min_hours, tstop = max_hours,
    options = list(idname = "serial_number"))

  # Add time-varying covariates
  covs <-
    c("smart_5", "smart_10", "smart_12",
      "smart_187", "smart_188", "smart_189",
      "smart_196", "smart_197", "smart_198")
  expres <- parse(
    text = paste(
      "tmerge(
      final_tbl, surv_tbl, id = serial_number,\n",
      paste0(covs, " = tdc(smart_9_raw, ", covs, "_raw)", collapse = ",\n"),
      ", options=list(na.rm=FALSE))"))

  final_tbl <- eval(expres)

  # summary(final_tbl)

  # TODO: Why are there extreme counts of smart_12? The disk cannot have had that many power cycles!
  # TODO: Should we be aware of extreme values for some of the others?
  # TODO: How to treat tha NAs / missing fields for smart_187 to smart_189
  # TODO: Are failures also when they change HDs due to some of the covariates? If so, which?
  # TODO: what is the scan errors in 'Failure Trends in a Large Disk Drive Population'
  # TODO: fill in the a zero dummy value for those w/ first observation that end within the first week of running
  # TODO: Some HDs seems to jump in time (like W300R8HJ). Smart_12 though seems consistent with smart_9 (and not the date)

  # Add binary features and cum sums
  for(n in covs[covs %in% c("smart_197")]){
    final_tbl[[paste0(n, "_bin")]] <- final_tbl[[n]] > 0

    final_tbl[[paste0(n, "_bin_cumsum")]] <- unlist(tapply(
      final_tbl[[paste0(n, "_bin")]], final_tbl$serial_number, cumsum))

    final_tbl[[paste0(n, "_cumsum")]] <- unlist(tapply(
      final_tbl[[n]], final_tbl$serial_number, cumsum))
  }

  saveRDS(final_tbl, final_tbl_name)

  #####
  # Check final tbl versus data from the orginal data base

  set.seed(467322)
  ids <- sample(unique(final_tbl$serial_number), 25)

  raw_data <- dbGetQuery(
    conn,
    paste0("Select * from drive_stats where serial_number IN (", paste0(
      "'", ids, "'", collapse = ", "
    ), ")"))

  raw_data <- raw_data[order(raw_data$serial_number, raw_data$smart_9_raw), ]

  n <- nrow(raw_data)
  keep <- c(
    T,
    raw_data$serial_number[-n] != raw_data$serial_number[-1] |
      raw_data$smart_5_raw[-n] != raw_data$smart_5_raw[-1] |
      raw_data$smart_10[-n] != raw_data$smart_10[-1] |
      raw_data$smart_12[-n] != raw_data$smart_12[-1] |
      raw_data$smart_187[-n] != raw_data$smart_187[-1] |
      raw_data$smart_188[-n] != raw_data$smart_188[-1] |
      raw_data$smart_189[-n] != raw_data$smart_189[-1] |
      raw_data$smart_196[-n] != raw_data$smart_196[-1] |
      raw_data$smart_197[-n] != raw_data$smart_197[-1] |
      raw_data$smart_198[-n] != raw_data$smart_198[-1])

  raw_data <- raw_data[keep, ]
  final_tbl_sub <- final_tbl[final_tbl$serial_number %in% ids, ]
  final_tbl_sub$serial_number <- as.character(final_tbl_sub$serial_number)
  final_tbl_sub <- final_tbl_sub[order(final_tbl_sub$serial_number), ]

  raw_data[, -1] <- apply(
    raw_data[, -1], 2, function(x){
      if(is.character(x)){
        x[x == ""] <- NA_character_
        x <- as.numeric(x)
        return(x)
      }
      x
    })

  row.names(raw_data) <- NULL
  row.names(final_tbl_sub) <- NULL

  test_that("Table in R match with final survival table", {
    expect_equal(nrow(raw_data), nrow(final_tbl_sub))
    expect_equal(raw_data$serial_number,
                 as.character(final_tbl_sub$serial_number))

    expect_equal(raw_data$smart_9_raw, final_tbl_sub$tstart)

    tmp_names <- colnames(final_tbl_sub)[grepl(
      "^smart_\\d+$", colnames(final_tbl_sub))]
    expect_equal(raw_data[, paste0(tmp_names, "_raw")],
                 final_tbl_sub[, tmp_names])
  })

}, finally = {
  if(!is.null(conn))
    dbDisconnect(conn)
  conn <- NULL

  for(l in tmp_files)
    unlink(l, recursive = T)
  tmp_files <- list()

  if(!is.null(pb)){
    close(pb)
    pb <- NULL
  }

  eval(reset_dir)

  sink()
  sink(file = NULL, type = "message")
  if(!is.null(log_file))
    close(log_file)

  options(warn = warn_old)

  rm(list = ls()[!ls() %in% cur_ls])
})


#####
# "VACUUM" the db to decrease the size
# See http://www.tutorialspoint.com/sqlite/sqlite_vacuum.htm
tryCatch({
  eval(set_dir)
  options(warn = 1)

  eval(conn_exp)

  dbSendQuery(conn, "VACUUM;")
}, finally = {
  if(!is.null(conn))
    dbDisconnect(conn)
})

# dbGetQuery(
#   conn,
#   "SELECT serial_number, date, smart_5_raw, smart_9_raw, smart_12_raw, smart_187_raw, smart_197_raw, smart_198_raw
# FROM drive_stats
#   where serial_number = 'WD-WCAU45409452'
#   ORDER BY smart_9_raw ASC")
