rm(list = ls())



res <- sapply(seq(nb_methods), function(imethod){
  apply((Ytilde_test_allh[, , , imethod] - Y_test_allh)^2, c(1, 2), mean)
}, simplify = "array")


res <- simplify2array(results)
res <- apply(res, c(1, 2, 3), mean)
par(mfrow = c(3, 3))
for(j in seq(7)){
  #matplot(res[, j, ], lty = 1, type = "l")
  matplot(res[, j, ])
}

res_overall <- apply(res, c(1, 3), mean)
matplot(res_overall)

# check aggregates and bottom