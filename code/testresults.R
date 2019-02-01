rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  experiment <- "tourism"
  #experiment <- "wikipedia"
  #experiment <- "small"; nbruns <- 100
  
  #algobf <- "ets"
  algobf <- "arima"
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  algobf     <- args[[2]]
  nbruns <- as.numeric(args[[3]])
  
}
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
library(stargazer)
library(Hmisc)

do.relative <- TRUE

info_file <- file.path(results.folder, 
                         paste("info_", experiment , ".Rdata", sep = ""))
load(info_file)

add.bias <- FALSE
lambda_selection <- "min"

if(grepl("small", experiment) || grepl("large", experiment)){
  do.log <- FALSE
  do.cleaning <- FALSE
  do.deseasonalization <- FALSE
  do.scaling <- FALSE
  
}else{
  do.log <- ifelse(experiment == "tourism", FALSE, TRUE)
  do.deseasonalization <- TRUE
  do.scaling <- FALSE
}

if(experiment %in% c("small", "large")){
  list_ids <- seq(nbruns)
  set_experiments <- paste(experiment, list_ids, sep = '-')
}else{
  set_experiments <- experiment
}

list_ferrors <- vector("list", length(set_experiments))
for(k in seq_along(set_experiments) ){
  exp <- set_experiments[k]
  myfile <- file.path(results.folder, paste("resultsicml_", exp, "_", lambda_selection, "_", 
                                            add.bias, "_", do.log, "_", 
                                            do.deseasonalization, "_", do.scaling, "_", algobf, ".Rdata", sep = ""))
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
      mat <- apply(errors[, id, drop = F], 1, sum)
      errors_agg[, i] <- apply(errors[, id, drop = F], 1, sum)
    }
    row.names(errors_agg) <- row.names(errors)
    list_ferrors[[k]] <- errors_agg
}

all_ferrors <- simplify2array(list_ferrors)

avg_ferrors <- apply(all_ferrors, c(1, 2), mean)
stderrors <- apply(all_ferrors, c(1, 2), sd)/sqrt(length(list_ferrors))
  
final_avg_ferrors <- avg_ferrors[c("BASE", "BU", "REG", "REG-PBU", 
                             "MINTshr", "MINTols"), ]
final_stderrors <- stderrors[c("BASE", "BU", "REG", "REG-PBU", 
                                 "MINTshr", "MINTols"), ]

colnames(final_avg_ferrors) <- paste("Level ", unique(niveaus), sep = "")
  

mymat <- final_avg_ferrors

if(do.relative){
  mymat <- 100 * t((t(mymat) - mymat["BASE", ])/mymat["BASE", ])
  mymat <- mymat[-which(row.names(mymat) == 'BASE'), ]
}


####################################
nbcol <- ncol(mymat)
nbrow <- nrow(mymat)
ids_row <- apply(mymat, 2, which.min)
#finalmat <- format(mymat, digits = 2, nsmall = 2)
finalmat <-mymat
for(j in seq(nbcol)){
  for(i in seq(nbrow)){
    finalmat[i, j] <- sprintf("%.2f", mymat[i, j])
    if(i == ids_row[j]){
      finalmat[i, j] <- paste("\\mathbf{", finalmat[i, j], "}", sep = "")
    }
    finalmat[i, j] <- paste("$", finalmat[i, j], "$", sep = "")
  }
}
####################################
print(finalmat)
#print(final_stderrors)

#tablefile <- paste("table", "_", experiment, "_", lambda_selection, "_", add.bias, "_", 
#             do.log, "_",  do.deseasonalization, "_", do.scaling, "_", algobf, ".tex", sep = "")
#stargazer(finalmat, title = paste(experiment, " - ", algobf, sep = ""), 
#          out = file.path(pdf.folder, tablefile),
#          float.env= "table*")

if(TRUE){
tablefile <- paste("_table", "_", experiment, "_", lambda_selection, "_", add.bias, "_", 
                   do.log, "_",  do.deseasonalization, "_", do.scaling, "_", algobf, ".tex", sep = "")
v <- latex(finalmat, star = TRUE,
      title = paste(experiment, " - ", algobf, sep = ""), 
      file = file.path(pdf.folder, tablefile))
}
stop("done")

  
print(errors_agg[c("BASE", "BU", "REG", "REG-PBU", "MINTshr", "MINTols", "L1", "L1-PBU"), ])
print(t(t(errors_agg) - errors_agg["BASE", ]))
x <- t(t(errors_agg) - errors_agg["BASE", ])
  #if(experiment == "tourism")
  #  x <- x / 1000
print(format(x[c("BASE", "BU", "REG", "REG-PBU", "MINTshr", "MINTols", "L1", "L1-PBU"), ], digits = 2))
  
  
