rm(list = ls())
source("config_paths.R")
library(readr)

load(file.path(main.folder, "work/rdata", "elecDT.Rdata"))
# DT
nseries <- ncol(DT)

mycount <- numeric(nseries) + 1
n <- 100
for(ihts in seq(50)){
  print(ihts)
  set.seed(1000 + ihts - 1)
  
  idselected <- sample(nseries, size = n, prob = mycount/sum(mycount))
  mycount[idselected] <- mycount[idselected] + 1
  
  Z <- DT[, idselected]
  cnames <- paste(sample(rep(c("AA", "AB", "AC", "AD", "AE", "BA", "BB", "BC", "BD", "BE"), each = 10)), 
                  sprintf("%04d", seq(n)), sep = "")
  colnames(Z) <- cnames

  save(file = file.path(rdata.folder, paste("elec-", ihts,".Rdata", sep = "")), list = c("Z", "idselected"))
}


