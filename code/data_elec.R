rm(list = ls())
library(readr)
library(forecast)
source("nicefigs.R")
source("config_paths.R")

filepath <- "../data/LD2011_2014.txt"
con = file(filepath, "r")

header <- readLines(con, n = 1)

i <- 0
line_read <- vector("list", 140256)
res <- vector("list", 140256)
dates_read <- numeric(140256)
while ( TRUE ) {
  if(i%%1000 == 0)
    print(i)
  x <- readLines(con, n = 1)
  if ( length(x) == 0 ) {
    break
  }
  i <- i + 1
  line_read[[i]] <- x
  z <- unlist(strsplit(x, ";"))
  res[[i]] <- as.numeric(sub(",", ".", z[-1], fixed = TRUE))
  dates_read[i]  <- z[1]
}
close(con)

mydates <- substring(dates_read, 4, 20)
library(lubridate)
mytime <- ymd_hms(mydates)


mat <- t(sapply(res, identity))

idselected <- which(year(mytime) == 2014)
y <- mat[idselected, ]
time_y <- mytime[idselected]
vec <- apply(y == 0, 2, sum)
obj_sort <- sort(vec, index = T)
id_series <- which(vec < 100)
Z <- y[, id_series]

time_y
Z
doy <- yday(time_y)
table(as.numeric(table(doy)))

newZ <- sapply(seq(ncol(Z)), function(j){
  as.numeric(tapply(Z[, j], factor(doy), sum))
})
id_remove <- c(5, 34, 37, 52, 67, 73, 100, 107, 114, 115, 125, 310, 313)
DT <- newZ[, -id_remove]

save(file = file.path(main.folder, "work/rdata", "elecDT.Rdata"), list = c("DT", "time_y"))
