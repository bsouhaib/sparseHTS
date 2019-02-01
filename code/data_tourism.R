rm(list = ls())
source("config_paths.R")
library(hts)

X <- read.csv("../data/Tourism data_v3.csv")
DT <- X[, -seq(2)]

vecnames <- colnames(DT)
vec_purpose <- substring(vecnames, 4, 6)
purposes <- unique(vec_purpose)
for(purpose in purposes){
  id <- which(vec_purpose == purpose)
  vec <- vecnames[id]
  #print(vec)
  series_name <- substr(vec, 1, 3)
  series_name <-sapply(seq(length(series_name)), function(j){
    paste(series_name[j], sprintf("%04d", j), sep = "")
  })
  Z <- DT[, id]
  colnames(Z) <- series_name
  #y <- hts(Z, characters = c(1, 1, 1, 4))
  nb_zeroes <- apply(Z == 0, 2, sum)
  print(purpose)
  print(sort(nb_zeroes))
  
  id_selected <- which(nb_zeroes < 100)
  Z <- Z[, id_selected]
  
  myfile <- file.path(rdata.folder, paste("tourism-", purpose, ".Rdata", sep = ""))
  save(file = myfile, list = c("Z"))
}