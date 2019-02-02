rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  #experiment <- "tourism-Hol"  # "tourism-Oth" "tourism-Vis"   # "tourism-Hol"  # "tourism-Bus" 
  #experiment <- "wikipedia26"
  #experiment <- "road_traffic-3"
  #experiment <- "elec1"
  
  
  #experiment <- "wikipedia" ; nbruns <- 100;
  #experiment <- "tourism" ; nbruns <- 3;
  experiment <- "small"; nbruns <- 100
  #experiment <- "road_traffic" ; nbruns <- 100;
  #experiment <- "elec" ; nbruns <- 45;
  
  algobf <- "all"

}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  algobf     <- args[[2]]
  
  if(length(args)==3)
  nbruns <- as.numeric(args[[3]])
  
}
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
library(stargazer)
library(Hmisc)

do.simplify <- TRUE
do.skipoutliers <- FALSE
do.diffbase <- F

info_file <- file.path(results.folder, 
                         paste("info_", experiment , ".Rdata", sep = ""))
if(!file.exists(info_file)){
  exp <- unlist(strsplit(experiment, "-"))[1]
  info_file <- file.path(results.folder, 
                         paste("info_", exp , ".Rdata", sep = ""))
}
load(info_file)



index_agg <- seq(nrow(A))
index_bottom <- seq(nrow(A) + 1, nrow(A) + ncol(A))

add.bias <- FALSE
lambda_selection <- "min"

if(experiment == "small" || experiment == "large"){
  do.log <- FALSE
  do.cleaning <- FALSE
  do.deseasonalization <- FALSE
  do.scaling <- FALSE
  
}else{
  #do.log <- ifelse(experiment == "tourism", FALSE, TRUE)
  do.log <- TRUE
  
  do.deseasonalization <- TRUE
  do.scaling <- FALSE
}


if(experiment %in% c("small", "large", "wikipedia", "tourism", "road_traffic", "elec")){
  list_ids <- seq(nbruns)
  set_experiments <- paste(experiment, list_ids, sep = '-')
}else{
  set_experiments <- experiment
}

setalgos <- algobf
if(algobf == "all")
  setalgos <- c("arima", "ets")

current_finalmat <- NULL
for(algo_bforecasting in setalgos){

set_experiments <- set_experiments[-c(40)] 
list_ferrors <- vector("list", length(set_experiments))

for(k in seq_along(set_experiments) ){
  print(k)
  exp <- set_experiments[k]
  myfile <- file.path(results.folder, paste("resultsicml_", exp, "_", lambda_selection, "_", 
                                            add.bias, "_", do.log, "_", 
                                            do.deseasonalization, "_", do.scaling, "_", algo_bforecasting, ".Rdata", sep = ""))
    if(!file.exists(myfile)){
      print(myfile)
      stop("FILE DOES NOT EXIST")
    }
    load(myfile) 
    
    if(experiment == "small" || experiment == "large"){
      mynames <- names(results_h1) 
      mynames[which(mynames == "REG")] <- "ERMreg"
      mynames[which(mynames == "REGBU")] <- "ERMregbu"
      names(results_h1) <- mynames
    }
    
    FUTURE <- Y_test_h
    
    if(do.skipoutliers){
      tag <- paste(do.log, "_", do.deseasonalization, "_", do.scaling, sep = "")
      datasave_file <- file.path(results.folder, paste("datasave_", exp, "_", tag ,".Rdata", sep = ""))
      load(datasave_file)
      
      idoutliers <- lapply(outindex, function(setoutliers){
        idx <- which(setoutliers > T_train + T_valid)
        ret <- NULL
        if(length(idx) > 0){
          ret <- setoutliers[idx] - (T_train + T_valid) 
        }
        ret
      })
      
      m <- lapply(seq(ncol(FUTURE)), function(j){ tsoutliers(ts(FUTURE[, j], freq = 7))$index }) 
    }
    #print(k)
    res <- sapply(results_h1, function(results_method){
      #print(dim(results_method))
    matres <- (results_method - FUTURE)^2
      if(do.skipoutliers){
        matres <- sapply(seq(ncol(matres)), function(j){ 
          x <- matres[, j]
          if(j > nrow(A)){
            idna <- idoutliers[[j - nrow(A)]] - 1 # because we start at two
            if(length(idna) > 0){
              x[idna] <- NA
            }
          }
          x
        })

      }
      matres
    }, simplify = "array")
    
    errors <- t(apply(res, c(2, 3), mean, na.rm = T))
    
    if(!do.simplify){
      allniveaus <- unique(niveaus)
      errors_agg <- matrix(NA, nrow = nrow(errors), ncol = length(allniveaus))
      for(i in  seq_along(allniveaus) ){
        niveau <- allniveaus[i]
        id <- which(niveau == niveaus)
        mat <- apply(errors[, id, drop = F], 1, sum)
        errors_agg[, i] <- apply(errors[, id, drop = F], 1, sum)
      }
      errors_agg <- cbind(apply(errors, 1, sum), errors_agg)
    }else{
      errors_agg <- apply(errors[, index_bottom], 1, sum)
      #errors_agg <- cbind(apply(errors[, index_agg], 1, sum), errors_agg)
      errors_agg <- cbind(apply(errors, 1, sum), errors_agg)
    }
    

    
    row.names(errors_agg) <- row.names(errors)
    # print("problem was here")
    # row.names(errors_agg) <-  c("ERM", "ERMreg", "ERMregbu", "BU", "MINTshr", "MINTsam", "MINTols")
    list_ferrors[[k]] <- errors_agg
}


# check errors wikipedia
if(FALSE){
  
  sort(sapply(list_ferrors, function(mat){mat["ERMregbu", 1]}) - 
    sapply(list_ferrors, function(mat){mat["MINTshr", 1]}) )
  
  tag <- paste(do.log, "_", do.deseasonalization, "_", do.scaling, sep = "")
  refit_step <- 14
  
  datasave_file <- file.path(results.folder, paste("datasave_", experiment, "_", tag ,".Rdata", sep = ""))
  load(datasave_file)
  
  obj <- sort(errors["MINTshr", index_bottom] - errors["REG-PBU", index_bottom], index = T)
  set <- head(obj$ix, 10)
  #par(mfrow = c(3, 1))
  for(j in set){
    #c("Z", "bts", "cleaned_bts")
    plot.ts(Z[, j]); abline(v = cumsum(c(86, 160, 120)), col = "red")
    plot.ts(cleaned_bts[, j]); abline(v = cumsum(c(86, 160, 120)), col = "red")
    plot.ts(bts[, j]); abline(v = cumsum(c(86, 160, 120)), col = "red")
    browser()
  }
  
  stop("done")
  
  
  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, "_", refit_step, "_", tag, "_", algo_bforecasting, ".Rdata", sep = ""))  
  load(file_bf)
  Yhat_valid_allh <- data_valid$Yhat
  Y_valid_allh     <- data_valid$Y
  sqd_errors <- t( (Yhat_valid_allh[1, , ] - Y_valid_allh[1, , ])^2 )
  obj <- sort(apply(sqd_errors, 2, mean), index = T, decreasing = T)
  for(j in obj$ix){
    if(j>111){
      plot.ts(Y_valid_allh[1, j, ], main = paste(j, j - 111))
      lines(Yhat_valid_allh[1, j, ], col = "red")
      browser()
    }
  }
  
  datasave_file <- file.path(results.folder, paste("datasave_", experiment, "_", tag ,".Rdata", sep = ""))
  load(datasave_file)
  obj <- sort(errors["MINTshr", index_bottom] - errors["REG-PBU", index_bottom], index = T)
  set <- head(obj$ix, 10)
  
  j <- set[2]

  plot.ts(FUTURE[, index_bottom[j]])
  lines(results_h1[["REG-PBU"]][, index_bottom[j]], col = "purple")
  lines(results_h1[["BASE"]][, index_bottom[j]], col = "blue")
  lines(results_h1[["MINTshr"]][, index_bottom[j]], col = "red")
   
  matplot(cbind(FUTURE[, index_bottom[j]], 
                 results_h1[["REG-PBU"]][, index_bottom[j]], 
                results_h1[["BASE"]][, index_bottom[j]]), type = 'l', col = c("black", "red", "blue"))
  
  plot.ts(Y_valid_allh[1, index_bottom[j], ])
  lines(Yhat_valid_allh[1, index_bottom[j], ], col = "red")
  
  stop("done")
  
  
  
  
}


if(do.diffbase){
  list_ferrors <- lapply(list_ferrors, function(mat){
    t(t(mat) - mat["BASE", ])
  })
}

all_ferrors <- simplify2array(list_ferrors)
avg_ferrors <- apply(all_ferrors, c(1, 2), mean)
stderrors <- apply(all_ferrors, c(1, 2), sd)/sqrt(length(list_ferrors))

if(experiment == "small" || experiment == "large"){
  order_methods <- c("BASE", "BU", "ERM", "ERMreg", "ERMregbu", "MINTsam", "MINTols", "MINTshr")
  nrgroup <- c(2, 3, 3)
}else{
  order_methods <- c("BASE", "BU", "ERMreg", "ERMregbu", "MINTols", "MINTshr")
  nrgroup <- c(2, 2, 2)
}

final_avg_ferrors <- avg_ferrors[order_methods, ]
final_stderrors <- stderrors[order_methods, ]
if(do.simplify){
  ###  colnames(final_avg_ferrors) <- c(paste("All levels ", "(", sum(table(niveaus)), ")", sep = ''), 
  ###                                 paste("Upper levels ", "(", length(index_agg), ")" , sep = ""), 
  ###                                 paste("Bottom level ", "(", length(index_bottom), ")" , sep = "") )
  
  colnames(final_avg_ferrors) <- c("All", "Bottom")
}else{
  colnames(final_avg_ferrors) <- c(paste("All levels ", "(", sum(table(niveaus)), ")", sep = ''), 
                                   paste("Level ", unique(niveaus), 
                                         " (", as.numeric(table(niveaus)), ")", sep = ""))
}

if(do.diffbase){
  final_avg_ferrors <- final_avg_ferrors[-which(row.names(final_avg_ferrors) == 'BASE'), ]
  final_stderrors   <- final_stderrors[-which(row.names(final_stderrors) == 'BASE'), ]
}

###if(do.relative){
###  mymat <-  t(t(mymat) - mymat["BASE", ])
###  mymat <- mymat[-which(row.names(mymat) == 'BASE'), ]
###}
####################################
nbcol <- ncol(final_avg_ferrors)
nbrow <- nrow(final_avg_ferrors)
ids_row <- apply(final_avg_ferrors, 2, which.min)
#finalmat <- format(mymat, digits = 2, nsmall = 2)
finalmat <- final_avg_ferrors
for(j in seq(nbcol)){
  for(i in seq(nbrow)){
    finalmat[i, j] <- sprintf("%.2f", final_avg_ferrors[i, j])
    if(i == ids_row[j]){
      finalmat[i, j] <- paste("\\mathbf{", finalmat[i, j], "}", sep = "")
    }
    
    stde_part <- NULL
    if(experiment == "small" || experiment == "large" || experiment == "wikipedia" 
       || experiment == "road_traffic" || experiment == "elec"){
      stde_part <- paste("~(", sprintf("%.2f", final_stderrors[i, j]), ")", sep = "")
    }
    finalmat[i, j] <- paste("$", finalmat[i, j], stde_part, "$", sep = "")
  }
}
####################################

current_finalmat <- cbind(current_finalmat, finalmat)

}

finalmat <- current_finalmat
print(finalmat)

tablefile <- paste("_table", "_", experiment, "_", lambda_selection, "_", add.bias, "_", 
                   do.log, "_",  do.deseasonalization, "_", do.scaling, "_", algobf, ".tex", sep = "")
v <- latex(finalmat, star = TRUE, 
           #title = paste(experiment, " - ", algo_bforecasting, sep = ""),
           cgroup = toupper(setalgos),
           n.cgroup = rep( ncol(final_avg_ferrors), length(setalgos)),
           rgroup=NULL, n.rgroup= nrgroup,
      title = "",
      file = file.path(pdf.folder, tablefile), booktabs = T, table.env = F, center = "none")

  

  
