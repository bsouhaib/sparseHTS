pbu <- function(objhts){
  
  #0_{m x (n-m)} I_m
  P_BU <- Matrix(0, nrow = objhts$nbts, ncol = objhts$nts, sparse = T)
  P_BU[cbind(seq(objhts$nbts), seq(objhts$nts - objhts$nbts + 1, objhts$nts)) ] <- 1
  P_BU
  #P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts))
}
  
fct_vec <- function(X){
  x <- X
  dim(x) <- c(dim(x)[1] * dim(x)[2], 1)
  x
}

fct_inv_vec <- function(x, nrows, ncolumns){
  matrix(x, nrow = nrows, ncol = ncolumns)
}

makeX <- function(Yhat, obj){
  Matrix::kronecker(obj$S, Yhat)
}
makey <- function(Y){
  y <- Y
  dim(y) <- c(dim(y)[1] * dim(y)[2], 1)
  y
}

makeY <- function(y, obj){
  if(!is.vector(y))
    stop("y must be a vector")
  matrix(y, ncol = obj$nts)
  # ANY PROBLEM HERE?? always byrom ???
}

getP <- function(beta, obj){
  Pprime <- matrix(beta, nrow = obj$nts, ncol = obj$nbts)
  t(Pprime)
}

#Pfrom_beta <- function(beta){
#  Pprime <- matrix(beta, nrow = n, ncol = m) # n x m
#  P <- t(Pprime)
#}


mint_betastar <- function(W, y_hat){
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  C2 <- MAT3 %*% y_hat
  adj <- C1 %*% C2
  -adj
}

pmint <- function(W){
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  P_BU - C1 %*% MAT3
}