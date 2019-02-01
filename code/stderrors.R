# myid <- sample(dim(res)[1], dim(res)[1], replace = T)
# errors <- t(apply(res[myid, , ], c(2, 3), mean))

myid <- sample(dim(res)[1], dim(res)[1], replace = T)
x <- sample(res_dimensions[1], res_dimensions[1], replace = T)
y <- seq(res_dimensions[2]) 
errors <- t(apply(res[x, y, ], c(2, 3), mean))
res_dimensions <- dim(res)



myarray <- sapply(seq(1000), function(k){
  allniveaus <- unique(niveaus)
  errors_agg <- matrix(NA, nrow = nrow(errors), ncol = length(allniveaus))
  for(i in  seq_along(allniveaus) ){
    niveau <- allniveaus[i]
    id <- which(niveau == niveaus)
    
    res_niveau <- res[, id, , drop = F]
    dimres <- dim(res_niveau)

    x <- sample(dimres[1], dimres[1], replace = T)
    y <- seq(dimres[2]) 
    errors <- t(apply(res_niveau[x, y, , drop = F], c(2, 3), mean))
    
    mat <- apply(errors, 1, sum)
    errors_agg[, i] <- mat
  }
  errors_agg
}, simplify = "array")

print(apply(myarray, c(1, 2), sd))


