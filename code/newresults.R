rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)

#experiment <- "tourism"
#experiment <- "wikipedia"

experiment <- "small-1"

lambda_selection <- "min"
add.bias <- FALSE
do.log <- ifelse(experiment == "tourism", FALSE, TRUE)
do.deseasonalization <- TRUE
do.scaling <- TRUE

info_file <- file.path(results.folder, paste("info_", experiment, ".Rdata", sep = ""))
load(info_file)


myfile <- file.path(results.folder, paste("resultsicml_", experiment, "_", lambda_selection, "_", 
                                          add.bias, "_", do.log, "_", 
                                          do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))
  if(!file.exists(myfile)){
    stop("FILE DOES NOT EXIST")
  }
  load(myfile) 
  
  FUTURE <- Y_test_h
  res <- sapply(results_h1, function(results_method){
    (results_method - FUTURE)^2
  }, simplify = "array")
  
  errors <- t(apply(res, c(2, 3), mean))
  
  allniveaus <- unique(niveaus)
  errors_agg <- matrix(NA, nrow = nrow(errors), ncol = length(allniveaus))
  for(i in  seq_along(allniveaus) ){
    niveau <- allniveaus[i]
    id <- which(niveau == niveaus)
    #print(id)
    mat <- apply(errors[, id, drop = F], 1, sum)
    errors_agg[, i] <- apply(errors[, id, drop = F], 1, sum)
  }
  row.names(errors_agg) <- row.names(errors)
  print(errors_agg[c("BASE", "BU", "REG", "REG-PBU", "MINTshr", "MINTols", "L1", "L1-PBU"), ])

  print(t(t(errors_agg) - errors_agg["BASE", ]))
  
  x <- t(t(errors_agg) - errors_agg["BASE", ])
  #if(experiment == "tourism")
  #  x <- x / 1000
  print(format(x[c("BASE", "BU", "REG", "REG-PBU", "MINTshr", "MINTols", "L1", "L1-PBU"), ], digits = 2))
  
  
