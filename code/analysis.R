rm(list = ls())

nb_methods <- 8

towards_pbu <- TRUE
#towards_pbu <- FALSE
lambda_selection <- "min"
#lambda_selection <- "1se"

njobs <- 5 # 14
nb_simulations <- 10

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
        res <- sapply(seq(nb_methods), function(imethod){
          apply((Ytilde_test_allh[, , , imethod] - Y_test_allh)^2, c(1, 2), mean)
         }, simplify = "array")
        results[[nbfiles]] <-  res
    }
  }
}

v <- Reduce("+", results[nbfiles])/nbfiles
v <- apply(v, c(1, 3), mean)
matplot(v[seq(3), c(4, 5, 7, 8)])

u <- apply(results[[3]], c(1, 3), mean)
matplot(u[seq(3), c(4, 5, 7, 8)])


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