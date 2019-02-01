rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)

hierarchy_name <- "tourism"
#hierarchy_name <- "wikipedia"

do.log <- FALSE
do.deseasonalization <- TRUE
do.scaling <- FALSE
add.bias <- FALSE
#lambda_selection <- "1se"
lambda_selection <- "min"

horizon_of_interest <- 1

info_file <- file.path(results.folder, paste("info_", hierarchy_name, ".Rdata", sep = ""))
load(info_file)

naggts <- nrow(A)
nbts <- ncol(A)
myfile <- file.path(results.folder, paste("resultsicml_", hierarchy_name, "_", 
                                          lambda_selection, "_", add.bias,  "_", do.log, "_", 
                                          do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))
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

experiment <- hierarchy_name
obj_sort <- sort( (err_reg - err_mintshrink)[index_bottom], index = T)

datasave_file <- file.path(results.folder, paste("datasave_", experiment, "_", 
                                                 do.log, "_", do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))
load(datasave_file)
savepdf("test")
par(mfrow = c(3, 1))

ids <- rev(tail(obj_sort$ix, 10))
for(j in ids){
  plot.ts(Z[, j])
  plot.ts(cleaned_bts[, j])
  plot.ts(bts[, j])
  abline(v = c(86, 86 + 160, 86 + 160 + 120), col = "red")
}
dev.off()

########## P/C matrix

plotmat <- function(M){
  X <- as.matrix(M)
  X <- t(apply(X, 2, rev))
  Xmelted <- melt(X)
  ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= factor(value) )) + 
    geom_tile() + 
    scale_fill_manual(values = c('white','blue')) +
    theme_bw()
}

plotmat(P_BU)
C_REG <- P_REG - P_BU
plotmat( (C_REG != 0)+0 )
#C_MINTSHRINK <- P_MINTSHRINK - P_BU
#plotmat( (C_MINTSHRINK != 0)+0 )

###### Forecast bias
res <- apply((ytrue_test - ypredbase_test)/ytrue_test, 2, mean) 
#res <- apply((Y_valid_h - Yhat_valid_h)/Y_valid_h, 2, mean) 
myindex <- index_bottom
plot(res[myindex][-c(277, 280)], (err_reg[myindex] - err_mintshrink[myindex])[-c(277, 280)] )
abline(h = 0, v = 0)

#heatmap(C, Colv = NA, Rowv = NA, scale="column",  col = colorRampPalette(brewer.pal(8, "Blues"))(25))

