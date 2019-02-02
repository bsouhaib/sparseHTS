rm(list = ls())
source("config_paths.R")
library(readr)

load(file.path(main.folder, "work/rdata", "elecDT.Rdata"))
# DT
nseries <- ncol(DT)

id_available <- seq(nseries)
n <- 100
for(ihts in seq(50)){
  print(ihts)
  set.seed(1000 + ihts - 1)
  
  id <- sample(id_available, n)
  Z <- DT[, id]
  #cnames <- paste(sample(rep(c("AA", "AB", "AC", "BA", "BB", "BC"), each = 50)), sprintf("%04d", seq(n)), sep = "")
  cnames <- paste(sample(rep(c("AA", "AB", "BA", "BB"), each = 25)), sprintf("%04d", seq(n)), sep = "")
  
  colnames(Z) <- cnames
  #id_available <- setdiff(id_available, id)
  
  save(file = file.path(rdata.folder, paste("elec-", ihts,".Rdata", sep = "")), list = c("Z", "id"))
}


