rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
#######

#experiment <- #"tourism" #"meters" #

hierarchy_name <- "tourism"


do.agg <- TRUE
hselected <- seq(2)

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

#lambda_selection <- "1se"
lambda_selection <- "min"


info_file <- file.path(results.folder, paste("info_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
load(info_file)

if(FALSE){
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
}


###### INDEP
do.deseasonalization <- TRUE
add.bias <- TRUE
naggts <- nrow(A)
nbts <- ncol(A)
myfile <- file.path(results.folder, paste("resultsINDEP_", hierarchy_name, "_", 
                                          lambda_selection, "_", add.bias,  "_", 
                                          do.deseasonalization, ".Rdata", sep = ""))
load(myfile)

S <- rbind(A, diag(nbts))
ypredreg_test <- t(S %*% t(pred_reg))
ypredmintshrink_test <- t(S %*% t(pred_mintshrink))
ypredmintols_test <- t(S %*% t(pred_mintols))
ypredbu_test <- t(S %*% t(pred_bu))
ypredbase_test <- pred_base
ypredbasebiased_test <- pred_basebiased
ytrue_test <- t(Y_test_allh[horizon_of_interest, , ])
ypredjoint_test <- pred_joint

sqderrors_base <- (ytrue_test - ypredbase_test)^2
sqderrors_basebiased <- (ytrue_test - ypredbasebiased_test)^2
sqderrors_reg <- (ytrue_test - ypredreg_test)^2
sqderrors_mintshrink <- (ytrue_test - ypredmintshrink_test)^2
sqderrors_mintols <- (ytrue_test - ypredmintols_test)^2
sqderrors_bu <- (ytrue_test - ypredbu_test)^2
sqderrors_joint <- (ytrue_test - ypredjoint_test)^2

err_base <- apply(sqderrors_base, 2, mean)
err_basebiased <- apply(sqderrors_basebiased, 2, mean)
err_reg <- apply(sqderrors_reg, 2, mean)
err_mintshrink <- apply(sqderrors_mintshrink, 2, mean)
err_mintols <- apply(sqderrors_mintols, 2, mean)
err_bu <- apply(sqderrors_bu, 2, mean)
err_joint <- apply(sqderrors_joint, 2, mean)

# STD ERRORS
x <- apply(t(sqderrors_base), 2, sum)
print(paste(mean(x), " - ", std.error(x)))

x <- apply(t(sqderrors_reg), 2, sum)
print(paste(mean(x), " - ", std.error(x)))

x <- apply(t(sqderrors_mintshrink), 2, sum)
print(paste(mean(x), " - ", std.error(x)))

x <- apply(t(sqderrors_mintols), 2, sum)
print(paste(mean(x), " - ", std.error(x)))

x <- apply(t(sqderrors_bu), 2, sum)
print(paste(mean(x), " - ", std.error(x)))

x <- apply(t(sqderrors_joint), 2, sum)
print(paste(mean(x), " - ", std.error(x)))



#Yhat_test_h  <- t(Yhat_test_allh[horizon_of_interest, , ])
#errors_base <- (ytrue_test - Yhat_test_h)^2
#err_base <- apply(errors_base, 2, mean)

#err_base <- apply(errors[[1]][, , "BASE"], 2, mean)
#err_mintshr <- apply(errors[[1]][, , "MINTshr"], 2, mean)
#err_bu <- apply(errors[[1]][, , "BU"], 2, mean)
#err_l1 <- apply(errors[[1]][, , "L1"], 2, mean)

savepdf(file = "errors")
par(mfrow = c(2, 1))

index_bottom <- seq(naggts + 1, naggts + nbts)
index_agg    <- seq(1, naggts)
plot(err_reg[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "METHOD")
plot(err_reg[index_bottom] - err_base[index_bottom], type = 'h', , ylab = "err", main = "METHOD")
plot(err_reg[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "METHOD")
plot( (err_reg[index_bottom] - err_base[index_bottom])[-251], type = 'h', , ylab = "err", main = "METHOD")


#plot(err_l1[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "L1")
#plot(err_l1[index_bottom] - err_base[index_bottom], type = 'h', ylab = "err", main = "L1")

#plot(err_mintshr[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "MINTshr")
#plot(err_mintshr[index_bottom] - err_base[index_bottom], type = 'h', ylab = "err", main = "MINTshr")

plot(err_mintshrink[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "MINT")
plot(err_mintshrink[index_bottom] - err_base[index_bottom], type = 'h', ylab = "err", main = "MINT")

plot(err_bu[index_agg] - err_base[index_agg], type = 'h', ylab = "err", main = "BU")
plot(err_bu[index_bottom] - err_base[index_bottom], type = 'h', ylab = "err", main = "BU")


dev.off()

j <- 1
plot.ts(ytrue_test[, j], main = j)
lines(pred_base[, j], col = "blue")
lines(pred_basebiased[, j], col = "orange")
lines(ypredmintshrink_test[, j], col = "red")
lines(ypredmintols_test[, j], col = "darkgreen")
lines(ypredreg_test[, j], col = "purple")

head(err_mintshrink - err_base, 5)
head(err_mintols - err_base, 5)
head(err_reg - err_base, 5)

D <- cbind((ytrue_test[, j] - ypredreg_test[, j])^2,
           (ytrue_test[, j] - ypredmintshrink_test[, j])^2,
           (ytrue_test[, j] - ypredmintols_test[, j])^2,
           (ytrue_test[, j] - ypredbase_test[, j])^2)
matplot(x = seq(nrow(D)), D, pch = 21, col = c("purple", "red", "darkgreen", "black"), type = 'p', cex = 0.4)

for(k in seq(2)){
  if (k == 1){
    myindex <- index_agg
    print("Aggregated series")
  }else{
    myindex <- index_bottom
    print("Bottom level series")
  }
  print(paste("BU        : ", prettyNum(sum(err_bu[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("BASE      : ", prettyNum(sum(err_base[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("BASE bias : ", prettyNum(sum(err_basebiased[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("MINTSHRINK: ", prettyNum(sum(err_mintshrink[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("MINTOLS   : ", prettyNum(sum(err_mintols[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("REG       : ", prettyNum(sum(err_reg[myindex]), big.mark = ",", scientific = FALSE)) )
  print(paste("REG JOINT : ", prettyNum(sum(err_joint[myindex]), big.mark = ",", scientific = FALSE)) )
  
}


###### ###### ###### ###### ###### 
stop("done")

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
par(mfrow = c(1, 2))

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

stop("done")

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
