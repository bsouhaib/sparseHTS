rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
#######

#experiment <- #"tourism" #"meters" #

hierarchy_name <- "UKF13"
hierarchy_name <- "UKF15"
hierarchy_name <- "UKF16"
hierarchy_name <- "UKF21"
hierarchy_name <- "UKF23"

# hierarchy_names <- c("UKF23", "UKF16", "UKF21", "UKF13", "UKF15")

do.agg <- TRUE
hselected <- seq(48)

color_methods <- c("orange", "yellowgreen", "purple", "deeppink", "pink", "yellow",
                   "darkseagreen1", "darkblue", "darkgreen", "red", "cyan", "blue", "aquamarine", "darkseagreen1", "darkseagreen1", "darkseagreen1", "darkseagreen1", 
                   "slategray4", "slategray", "aquamarine")
nb_methods <- length(color_methods)
#######
methods_toprint <- c("BASE2", "BU", "MINTshr", "MINTols", "MINTsam",  "OLS", "L2-PBU", "L1-PBU", "L1", "L2")
methods_toprint <- c("BASE2", "BU", "MINTshr", "MINTols", "L2-PBU", "L1-PBU", "L1", "L2", "NEW")

#methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "L2-PBU", "L1-PBU")
#c("L1-PBU", "BU", "MINTols", "MINTshr", "BASE","BASE2", "L2-PBU")

horizon_of_interest <- 1
do.ratio <- FALSE
tag <- "NIPS"

lambda_selection <- "1se"
#lambda_selection <- "min"


info_file <- file.path(results.folder, paste("info_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
load(info_file)


    myfile <- file.path(results.folder, paste("results_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
    if(!file.exists(myfile)){
     stop("FILE DOES NOT EXIST")
    }
      load(myfile) 
     
      H <- length(results_allh)
      errors <- lapply(seq(H), function(h){
        FUTURE <-t(Y_test_allh[h, , ])
        
        res <- sapply(results_allh[[h]], function(results_method){
          #(results_method$predictions - FUTURE)^2	
	          (results_method - FUTURE)^2
        }, simplify = "array")
        #browser()
        #print(dim(res))
        #print(length(names(results_allh[[h]])))
        #dimnames(res)[3] <- names(results_allh[[h]])
        res
        #PRED <- simplify2array(lapply(results_allh[[h]], "[[", "predictions"))
      })
    
errors_all <- errors

#errors_all[[1]]
#print(apply(apply(errors_all[[1]], c(1, 3), sum), 2, mean))
#print(apply(apply(errors_all[[1]][, seq(17), ], c(1, 3), sum), 2, mean))
#print(apply(apply(errors_all[[1]][, seq(18, 55), ], c(1, 3), sum), 2, mean))


name_methods <- dimnames(errors_all[[1]])[[3]]
naggts <- nrow(A)
nbts <- ncol(A)
nts <- naggts + nbts

myfile <- paste(tag, "_", hierarchy_name, "_", do.agg, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.9, width = 21 * 0.9)
par(mfrow = c(4, 4))

for(h in hselected){
  
  if(do.agg){
    v <- errors_all[[h]][, seq(naggts), ]
  }else{
    v <- errors_all[[h]][, seq(naggts + 1, nts), ]
  }

u <- apply(v, c(1, 3), sum)
mse_mean_h <- apply(u, 2, mean)
mse_std_h <- apply(u, 2, std.error)

plotCI(mse_mean_h, uiw = mse_std_h, liw = mse_std_h, main = paste("Horizon h = ", h, sep = ""), ylab = "MSE", xaxt = "n", xlab = "")
axis(1, at = seq(length(mse_mean_h)) , labels = name_methods, col.axis="blue", las=2)
}
dev.off()

myfile <- paste(tag, "_", hierarchy_name, "_", "FINAL", sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.4, width = 21 * 0.9)
par(mfrow = c(1, 2))

#listhorizons <- list(night = seq(1, 12), day = seq(13, 48))
listhorizons <- list(night = seq(1, 12), day = seq(13, max(hselected) ))

#listgroups <- list(agg = seq(naggts) , bot = seq(naggts + 1, nts))
#listgroups <- list(all = seq(nts), agg = seq(naggts) , bot = seq(naggts + 1, nts))
listgroups <- list(all = seq(nts))

for(ihgroup in seq_along(listhorizons)){
  for(i in seq_along(listgroups)){
    mygroup <- listgroups[[i]]
    myhorizons <- listhorizons[[ihgroup]]
  
    m <- sapply(errors_all[hselected], function(v){ apply(v[, mygroup, ], c(1, 3), sum)   }, simplify = "array")
    res <- apply(m[, , myhorizons], c(1, 2), mean) # avg horizon
    mse_mean <- apply(res, 2, mean)
    mse_std <- apply(res, 2, std.error)
    plotCI(mse_mean, uiw = mse_std, liw = mse_std, ylab = "MSE", xaxt = "n", xlab = "")
  
    #mse_mean <- apply(m, 2, mean)
    #mse_std <- apply(m, 2, std.error)
    #plotCI(mse_mean, uiw = mse_std, liw = mse_std, ylab = "MSE", xaxt = "n", xlab = "")
  }
}
dev.off()

stop("done")

s <- apply(errors_all[[1]], c(1, 3), sum)
matplot(s[, c(1, 4)], type = 'l', col = c("red", "purple"))

if(FALSE){
  
  r <- t(sapply(errors_all, function(errors_h){
    apply(apply(errors_h, c(2, 3), mean), 2, sum)
  }))
  matplot(r, type = 'l', lty = 1)
}

####
h <- horizon_of_interest
errors_h <- errors_all[[h]]
 
errors_final <- apply(errors_h, c(2, 3), mean)

#errors_final <- lapply(errors_h, function(obj){
#  apply(obj, c(2, 3), mean)
#})

id.keep <- match(methods_toprint, dimnames(errors_final)[[2]])
id.base <- match("BASE2", dimnames(errors_final)[[2]])


myfile <- paste(tag, "_", ifelse(do.ratio, "ratio", "absolute"), "_", hierarchy_name, "_", lambda_selection, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.3)
par(mfrow = c(1, 3))
par(cex.axis=.4)


naggts <- nrow(A)
nbts <- ncol(A)
nts <- naggts + nbts
groupings <- list(all = rep(1, nts), agg = c(rep(1, naggts), rep(0, nbts)), bot = c(rep(0, naggts), rep(1, nbts)))

for(i in seq_along(groupings)){
  mygroup <- which(groupings[[i]] == 1)
  
  mat <- errors_final
  #v <- sapply(errors_final, function(mat){
    err_all <- mat[mygroup, id.keep]
    err_ref <- mat[mygroup, id.base]
    if(do.ratio){
      err <- (apply(err_all, 2,  sum) - sum(err_ref))/sum(err_ref)
    }else{
      err <- apply(err_all, 2,  sum)
    }
    err
  #}, simplify = "array")
  v <- err
  err_toplot <- t(v)
  
  boxplot(err_toplot, outline = F, 
          main = paste("Nb series: ", length(mygroup)), cex = .5, col = color_methods[id.keep])
}
dev.off()

stop("done")
#######
h <- 1
results <- results_allh[[1]]
#MAT <- results[[1]]$predictions
FUTURE <-t(Y_test_allh[1, , ])
v <- sapply(seq_along(results), function(imethod){results[[imethod]]$predictions}, simplify = "array")
for(j in seq(20, 100)){
  matplot(cbind(v[, j, c(1, 2, 3, 5)], FUTURE[, j]), col = c("darkblue", "red", "purple", "cyan", "black"), type = 'l', lty = 1)
  browser()
}


r <- sapply(seq(9), function(imethod){
  errors_h[, , imethod] - errors_h[, , 5]
}, simplify = "array")
dimnames(r)[[3]] <- dimnames(errors_h)[[3]]

for(j in seq(100)){
  boxplot(r[, j, ])
  abline(h = 0, col = "red")
  browser()
}


r <- sapply(seq(9), function(imethod){
  errors_h[, , imethod] - errors_h[, , 4]
})



u <- apply(errors_h, c(1, 3), sum)
r <- sapply(seq(9), function(imethod){
  u[, imethod] - u[, 5]
})
