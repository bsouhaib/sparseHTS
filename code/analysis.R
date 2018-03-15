rm(list = ls())

nb_methods <- 10
name_methods <- c("LASSO", "SGLROW", "SGLCOL", "BU", "MINT", "LS", "NAIVE", "LSglmnet", "MINTOLS", "BASE")
color_methods <- c("red", "blue", "blue", "orange", "purple", "green", "brown", "darkgreen", "grey", "black")

#nb_methods <- 9
#name_methods <- c("LASSO", "SGLROW", "SGLCOL", "BU", "MINT", "LS", "NAIVE", "LSglmnet", "BASE")
#color_methods <- c("red", "blue", "blue", "orange", "purple", "green", "brown", "darkgreen", "black")

towards_pbu <- FALSE
#towards_pbu <- FALSE
lambda_selection <- "min"
#lambda_selection <- "1se"

njobs <- 2 # 14
nb_simulations <- 130

print("-----")
print(paste("towards_pbu: ", towards_pbu, sep = ""))
print(paste("lambda_selection: ", lambda_selection, sep = ""))
print("-----")



nbfiles <- 0
results <- vector("list", njobs * nb_simulations)

for(idjob in seq(njobs)){
  for(i in seq(nb_simulations)){
    myfile <- file.path(getwd(), "../work", paste("results_", idjob, ".", i, "_", towards_pbu, "_", lambda_selection, ".Rdata", sep = ""))
    if(file.exists(myfile)){
        nbfiles <- nbfiles + 1
        load(myfile) 
        print(dim(Ytilde_test_allh))
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

#####
methods_toprint <- c("LASSO", "BU", "MINT", "LS", "NAIVE", "LSglmnet", "BASE")
methods_toprint <- c("LASSO", "BU", "MINT", "LS", "LSglmnet", "BASE")
methods_toprint <- c("LASSO", "BASE", "BU", "MINT", "MINTOLS", "LSglmnet")
methods_toprint <- c("LASSO", "BASE", "BU", "MINT", "LSglmnet")
#methods_toprint <- c("LASSO", "MINT")
id.keep <- match(methods_toprint, name_methods)


#####
plot.splitted <- TRUE
if(TRUE){
par(mfrow = c(3, 3))
for(j in seq(7)){
  
  matplot(err_all[seq(10), j, id.keep], type = 'p', pch = 21, col = color_methods[id.keep], main = paste(towards_pbu, " - ", lambda_selection))
  if(j == 1){
    legend("topleft", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
  }
  #barplot(err_all[1, j, id.keep])
}
}else{
#####
matplot(err_avgnodes[seq(10), id.keep], type = 'p', pch = 21, col = color_methods[id.keep], main = paste(towards_pbu, " - ", lambda_selection))
#matpoints(err_avgnodes[, id.keep], col = color_methods[id.keep], type = 'p')
legend("topleft", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
}
#u <- apply(results[[1]], c(1, 3), mean)
#matplot(u[seq(3), c(4, 5, 7, 8)])
# 
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