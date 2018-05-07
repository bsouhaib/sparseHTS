rm(list = ls())
source("config_paths.R")
source("nicefigs.R")

#######
color_methods <- c("orange", "yellowgreen", "purple", "deeppink", "pink", "yellow",
                   "darkseagreen1", "darkblue", "darkgreen", "red", "cyan", "blue", "aquamarine", "darkseagreen1", "darkseagreen1", "darkseagreen1", "darkseagreen1", 
                   "slategray4", "slategray", "aquamarine")
nb_methods <- length(color_methods)
#######

methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "MINTsam",  "OLS", "L2-PBU", "L1-PBU")

horizon_of_interest <- 1
do.ratio <- FALSE
tag <- "nips"
experiment <- "tourism"
id_jobs <- 1 # seq(2000, 2060) #1986 #seq(400, 450) #420 #seq(200, 210) #420 #c(200, 210)
nb_simulations <- 1 #500
ids_simulations <- seq(nb_simulations) # 37
fmethod_agg <- "ETS"
fmethod_bot <- "ETS"
lambda_selection <- "1se" # "min"

info_file <- file.path(results.folder, paste("info_", experiment, "_", lambda_selection, ".Rdata", sep = ""))
load(info_file)

nbfiles <- 0
errors_simulations <- vector("list", length(id_jobs) * length(ids_simulations))

for(idjob in id_jobs){
  for(i in ids_simulations){
    if(i %% 50 == 0) print(i)
    myfile <- file.path(results.folder, paste("results_", experiment, "_", i, "_", lambda_selection, ".Rdata", sep = ""))
    if(file.exists(myfile)){
      nbfiles <- nbfiles + 1
      load(myfile) 
      
      H <- length(results_allh)
      errors <- lapply(seq(H), function(h){
        FUTURE <-t(Y_test_allh[h, , ])
        
        res <- sapply(results_allh[[h]], function(results_method){
          (results_method$predictions - FUTURE)^2
        }, simplify = "array")
        #browser()
        #print(dim(res))
        #print(length(names(results_allh[[h]])))
        #dimnames(res)[3] <- names(results_allh[[h]])
        res
        #PRED <- simplify2array(lapply(results_allh[[h]], "[[", "predictions"))
      })
      
      errors_simulations[[nbfiles]] <-  errors
      
    }
  }
}

errors_all <- errors_simulations[seq(nbfiles)] ## seq() !!!!!!!!!!!!!


h <- horizon_of_interest
errors_h <- lapply(errors_all, function(obj){
  obj[[h]]
})
errors_final <- lapply(errors_h, function(obj){
  apply(obj, c(2, 3), mean)
})

id.keep <- match(methods_toprint, dimnames(errors_final[[1]])[[2]])
id.base <- match("BASE2", dimnames(errors_final[[1]])[[2]])


myfile <- paste(tag, "_", ifelse(do.ratio, "ratio", "absolute"), "_", experiment, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.9)
par(mfrow = c(4, 2))
par(cex.axis=.4)


naggts <- nrow(A)
nbts <- ncol(A)
nts <- naggts + nbts
groupings <- list(all = rep(1, nts), agg = c(rep(1, naggts), rep(0, nbts)), bot = c(rep(0, naggts), rep(1, nbts)))

for(i in seq_along(groupings)){
  mygroup <- which(groupings[[i]] == 1)
  
  v <- sapply(errors_final, function(mat){
    err_all <- mat[mygroup, id.keep]
    err_ref <- mat[mygroup, id.base]
    if(do.ratio){
      err <- (apply(err_all, 2,  sum) - sum(err_ref))/sum(err_ref)
    }else{
      err <- apply(err_all, 2,  sum)
    }
    err
  }, simplify = "array")
  
  err_toplot <- t(v)
  
  boxplot(err_toplot, outline = F, 
          main = paste("Total - ", nbfiles), cex = .5, col = color_methods[id.keep])
}
dev.off()
