rm(list = ls())
source("nicefigs.R")

tag <- "biasVV"

experiment <- "small"
fmethod_agg <- "ARIMA"
fmethod_bot <- "ARIMA"

#fmethod_agg <- "AR1"
#lambda_selection <- "min"
lambda_selection <- "1se"


nb_methods <- 13
name_methods <- c("LAS", "LAS_PBU", "RID", "BU", "MINT", "LS", "NAIVE", "LSglmnet", "MINTOLS", "GGLAS", "GGLAS_PBU", "RID_PBU", "BASE")
color_methods <- c("red", "blue", "cyan", "orange", "purple", "green", "brown", "darkgreen", "grey", "pink", "darkseagreen1", "aquamarine", "black")

#nb_methods <- 10
#name_methods <- c("LASSO", "LASSO_TO_PBU", "RIDGE", "BU", "MINT", "LS", "NAIVE", "LSglmnet", "MINTOLS", "BASE")
#color_methods <- c("red", "blue", "cyan", "orange", "purple", "green", "brown", "darkgreen", "grey", "black")


id_jobs <- 1986 #seq(400, 450) #420 #seq(200, 210) #420 #c(200, 210)
nb_simulations <- 500
ids_simulations <- seq(nb_simulations)
#ids_simulations <- 37

print("-----")
print(paste("lambda_selection: ", lambda_selection, sep = ""))
print("-----")



nbfiles <- 0
results <- vector("list", length(id_jobs) * length(ids_simulations))

for(idjob in id_jobs){
  for(i in ids_simulations){
    if(i %% 50 == 0) print(i)
    myfile <- file.path(getwd(), "../work", paste("results_", experiment, "_", fmethod_agg, "_", fmethod_bot, "_", idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
    if(file.exists(myfile)){
      nbfiles <- nbfiles + 1
      load(myfile) 
      #print(dim(Ytilde_test_allh))
      res <- sapply(seq(nb_methods), function(imethod){
        apply((Ytilde_test_allh[, , , imethod] - Y_test_allh)^2, c(1, 2), mean)
      }, simplify = "array")
      results[[nbfiles]] <-  res
    }
  }
}

err <- sapply(seq(nbfiles), function(ifile){
  mean(results[[ifile]], na.rm = T)
})

#res_sort <- sort(err, decreasing = T, index = T)
#print(res_sort)
#myresults <- results[-res_sort$ix[seq(8)]]
res_sort <- sort(err, decreasing = T, index = T); print(res_sort)

myresults <- results[seq(nbfiles)] ## seq() !!!!!!!!!!!!!

err_all <- Reduce("+", myresults)/length(myresults)

err_sd <- sqrt(Reduce("+", lapply(myresults, function(mat){ (mat - err_all)^2}))/length(myresults))

err_avgnodes <- apply(err_all, c(1, 3), mean)
colnames(err_avgnodes) <- name_methods


v <- lapply(myresults, function(mat){ apply(mat, c(1, 3), mean) })
v_mean  <- Reduce("+", v)/length(v)
v_sd <- sqrt(Reduce("+", lapply(v, function(mat){ (mat - v_mean)^2}))/length(v))


#####
methods_toprint <- c("BU", "BASE", "MINT", "MINTOLS", "LAS_PBU", "LAS")
methods_toprint <- c("BU", "BASE", "MINT", "MINTOLS",  "LSglmnet", "RID", "LAS", "RID_PBU","LAS_PBU", "NAIVE")

id.keep <- match(methods_toprint, name_methods)


#####
htoplot <- 3
plot.splitted <- TRUE
if(FALSE){
if(TRUE){
  par(mfrow = c(3, 3))
  
  barplot(err_avgnodes[1, id.keep], cex.names = .5, main = "average over hts (h =1)")
  
  matplot(log(err_avgnodes[seq(htoplot), id.keep]), type = 'p', pch = 21, col = color_methods[id.keep], main = paste(lambda_selection, " - average over hts"), ylab = "MSE (log)")
  #matpoints(err_avgnodes[, id.keep], col = color_methods[id.keep], type = 'p')
  legend("bottomright", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
  
  for(j in seq(7)){
    
    matplot(log(err_all[seq(htoplot), j, id.keep]), type = 'p', pch = 21, col = color_methods[id.keep], main = j, ylab = "MSE (log)")
    if(j == 1){
      legend("bottomright", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
    }
    #barplot(err_all[1, j, id.keep])
  }
}else{
  #####
  matplot(err_avgnodes[seq(10), id.keep], type = 'p', pch = 21, col = color_methods[id.keep], main = paste(lambda_selection))
  #matpoints(err_avgnodes[, id.keep], col = color_methods[id.keep], type = 'p')
  legend("topleft", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
}
#u <- apply(results[[1]], c(1, 3), mean)
#matplot(u[seq(3), c(4, 5, 7, 8)])
# 
}

do.log <- FALSE
savepdf(paste("results_", experiment, "_", fmethod_agg, "_", fmethod_bot, "_", tag, sep = ""), height = 26 * 0.9)
par(mfrow = c(4, 2))
par(cex.axis=.4)

v <- sapply(myresults, function(mat){
  apply(mat[1, , id.keep], 2,  sum)
})
v <- t(v)
colnames(v) <- name_methods[id.keep]
final_v <- v
final_medians <- apply(v, 2, median)
if(do.log){
  final_v <- log(final_v)
  final_medians <- log(final_medians)
}
mymin <- min(final_medians)
boxplot(final_v, outline = F, 
        main = paste(fmethod_agg, "-", fmethod_bot, " - ", "Total - ", nbfiles), cex = .5, col = color_methods[id.keep])
abline(h = mymin, col = "red")

for(j in seq(7)){
  v <- sapply(myresults, function(mat){
    t(mat[1, j, id.keep])
  })
  v <- t(v)
  colnames(v) <- name_methods[id.keep]
  
  final_v <- v
  final_medians <- apply(v, 2, median)
  if(do.log){
    final_v <- log(final_v)
    final_medians <- log(final_medians)
  }
  mymin <- min(final_medians)
  
  boxplot(final_v, outline = F, main = j, col = color_methods[id.keep])
  abline(h = mymin, col = "red")
  #boxplot(v - v[, "BASE"], outline = F)
  #stop("done")
}
dev.off()

stop("done")

###############
res <- simplify2array(results)
res <- apply(res, c(1, 2, 3), mean)
par(mfrow = c(3, 3))
for(j in seq(nb_methods)){
  #matplot(res[, j, ], lty = 1, type = "l")
  matplot(res[, j, ])
}

res_overall <- apply(res, c(1, 3), mean)
matplot(res_overall)

# check aggregates and bottom
