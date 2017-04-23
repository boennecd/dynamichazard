library(RSQLite)
library(sqldf)
library(tcltk)

cur_ls <- ls()
current_wd <- getwd()
conn <- NULL
tmp <- list()
new_folder <- NULL
pb <- NULL
do_remove <- TRUE

tryCatch({
  setwd(paste0(
    stringr::str_match(getwd(), ".+dynamichazard"), "/vignettes/jss/HDs/"))

  db <- "driver_stats.db"
  if(file.exists(db) && do_remove)
    unlink(db, force = T)

  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, db)

  create_query <- paste0(readLines("create_table.sql"), collapse = "\n")
  dbSendQuery(conn, create_query)

  suffixes <- c("2013", "2014", "2015", "Q1_2016", "Q2_2016", "Q3_2016")
  for(suffix in suffixes){
    # TODO: Clean up
    # doc_url <- "https://f001.backblaze.com/file/Backblaze-Hard-Drive-Data/docs_2013.zip"
    data_url <- paste0("https://f001.backblaze.com/file/Backblaze-Hard-Drive-Data/data_", suffix, ".zip")

    # # Dl import query file
    # tmp$query <- tempfile()
    # download.file(doc_url, tmp$query)
    # import_sql <- paste0(readLines(
    #   unz(tmp$query, "docs_2013/import_2013.sql")),
    #   collapse = "\n")
    # unlink(tmp$query)
    # tmp$query <- NULL

    # Dl data files and unzip
    cat("Downloading files for", suffix, "\n")
    tmp$data <- tempfile()
    download.file(data_url, tmp$data)
    tmp$data_dir <- suffix
    if(substr(suffix, 0, 1) == "Q")
      tmp$data_dir <- paste0("data_", tmp$data_dir)
    unzip(tmp$data)
    unlink(tmp$data)
    tmp$data <- NULL

    # Import
    cat("Adding data to database for", suffix, "\n")
    fs <- dir(tmp$data_dir)
    pb <- txtProgressBar(
      1, length(fs), 0, title = paste("Progress for file suffix", suffix),
      style=3)
    for(i in seq_along(fs)){
      setTxtProgressBar(pb, i)
      f <- fs[i]
      f <- paste0("./", tmp$data_dir, "/", f)
      read.csv.sql(
        f, paste(
          "INSERT INTO drive_stats",
          # "SELECT date, serial_number, model, capacity_bytes, failure, smart_1_normalized, smart_1_raw, smart_2_normalized, smart_2_raw, smart_3_normalized, smart_3_raw, smart_4_normalized, smart_4_raw, smart_5_normalized, smart_5_raw, smart_7_normalized, smart_7_raw, smart_8_normalized, smart_8_raw, smart_9_normalized, smart_9_raw, smart_10_normalized, smart_10_raw, smart_11_normalized, smart_11_raw, smart_12_normalized, smart_12_raw, smart_13_normalized, smart_13_raw, smart_15_normalized, smart_15_raw, smart_183_normalized, smart_183_raw, smart_184_normalized, smart_184_raw, smart_187_normalized, smart_187_raw, smart_188_normalized, smart_188_raw, smart_189_normalized, smart_189_raw, smart_190_normalized, smart_190_raw, smart_191_normalized, smart_191_raw, smart_192_normalized, smart_192_raw, smart_193_normalized, smart_193_raw, smart_194_normalized, smart_194_raw, smart_195_normalized, smart_195_raw, smart_196_normalized, smart_196_raw, smart_197_normalized, smart_197_raw, smart_198_normalized, smart_198_raw, smart_199_normalized, smart_199_raw, smart_200_normalized, smart_200_raw, smart_201_normalized, smart_201_raw, smart_223_normalized, smart_223_raw, smart_225_normalized, smart_225_raw, smart_240_normalized, smart_240_raw, smart_241_normalized, smart_241_raw, smart_242_normalized, smart_242_raw, smart_250_normalized, smart_250_raw, smart_251_normalized, smart_251_raw, smart_252_normalized, smart_252_raw, smart_254_normalized, smart_254_raw, smart_255_normalized, smart_255_raw",
          "SELECT date, serial_number, model, capacity_bytes, failure, smart_1_raw, smart_2_raw, smart_3_raw, smart_4_raw, smart_5_raw, smart_7_raw, smart_8_raw, smart_9_raw, smart_10_raw, smart_11_raw, smart_12_raw, smart_13_raw, smart_15_raw, smart_183_raw, smart_184_raw, smart_187_raw, smart_188_raw, smart_189_raw, smart_190_raw, smart_191_raw, smart_192_raw, smart_193_raw, smart_194_raw, smart_195_raw, smart_196_raw, smart_197_raw, smart_198_raw, smart_199_raw, smart_200_raw, smart_201_raw, smart_223_raw, smart_225_raw, smart_240_raw, smart_241_raw, smart_242_raw, smart_250_raw, smart_251_raw, smart_252_raw, smart_254_raw, smart_255_raw",
          "FROM file"),
        dbname = db)
    }
    close(pb)
    unlink(tmp$data_dir, recursive = T, force = T)
  }

  }, finally = {
    setwd(current_wd)
    if(!is.null(conn))
      dbDisconnect(conn)

    for(l in tmp)
      unlink(l, recursive = T)

    if(!is.null(pb))
      close(pb)

    rm(list = ls()[!ls() %in% cur_ls])
  })

# db <- "driver_stats.db"
# sqlite <- dbDriver("SQLite")
# conn <- dbConnect(sqlite, db)
# dbExistsTable(conn, "drive_stats")
#
# rs = dbSendQuery(conn, "select * from drive_stats  LIMIT 100")
# d = fetch(rs, n = -1)
#
# head(d, 100)
#
# dbDisconnect(conn)
