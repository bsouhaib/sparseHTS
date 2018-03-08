

makeX <- function(Yhat, obj){
  Matrix(Matrix::kronecker(obj$S, Yhat), sparse = T)
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