

learnreg <- function(objreg, objhts, algo, config = NULL, selection = c("min", "1se"), towards_pbu = FALSE){
  
  X <- objreg$X
  if(towards_pbu){
    P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts))
    y <- objreg$y - X %*% makey(t(P_BU))
  }else{
    y <- objreg$y
  }
  
  
  lambda_final <- NULL
  s  <- NULL
  ############
  if(algo == "glmnet"){
    
    cvfit <- model <- do.call(cv.glmnet, c(list(x = X, y = y), config))
    s <- ifelse(selection == "min", "lambda.min", "lambda.1se")
    idmin_lambda <-  ifelse(selection == "min", match(model$lambda.min, model$lambda), match(model$lambda.1se, model$lambda))
    beta_hat <- as.numeric(coef(cvfit, s = s))
    if(!config$intercept)
      beta_hat <- beta_hat[-1]
  ############
  }else if(grepl("SGL", algo)){ 
    if(algo == "SGLrow"){
      grouping <- rep(seq(objhts$nts), objhts$nbts)
    }else if(algo == "SGLcol"){
      grouping <- rep(seq(objhts$nbts), each = objhts$nts)
    }
    
    data <- list(x = as.matrix(X) , y = y)
    param_list <- c(list(data = data, index = grouping, type = "linear"), config)
    cvfit <- do.call(cvSGL, param_list)
    
    idmin_lambda <- which.min(cvfit$lldiff)
    if(selection == "1se"){
      idmin_lambda <- min(which(cvfit$lldiff <= cvfit$lldiff[idmin_lambda] + cvfit$llSD))
    }
    
    param_list <- param_list[-which(names(param_list) == "nfold")]
    fit_SGL <- model <- do.call(SGL, param_list)
    beta_hat <- fit_SGL$beta[, idmin_lambda]
    lambda_final <- fit_SGL$lambdas[idmin_lambda]
  }
  
  obj_return <- list(algo = algo, model = model, idmin_lambda = idmin_lambda, s = s, towards_pbu = towards_pbu)
  if(towards_pbu){
    c_hat <- beta_hat
    C <- getP(c_hat, objhts)
    obj_return <- c(list(C = C), obj_return)
  }else{
    P <- getP(beta_hat, objhts)
    obj_return <- c(list(P = P), obj_return)
  }
  
  return(obj_return) 
}


predtest <- function(objreg, objhts, objlearn = NULL){
  X_test <- as.matrix(objreg$X)
  algo <- objlearn$algo
  
  if(algo == "glmnet" || grepl("SGL", algo)){  
    if(algo == "glmnet"){
      vecpred <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = objlearn$s))
    }else if(grepl("SGL", algo)){
      idmin_lambda <- objlearn$idmin_lambda
      vecpred <- as.numeric(predictSGL(objlearn$model, X_test, lam = idmin_lambda))
    }
    #else{
    #  vecpred <- X_test %*% t(objlearn$P)
    #}
    matpred <- makeY(vecpred, objhts)
    
    if(objlearn$towards_pbu){
      matpred <- matpred + predtest(objreg, objhts, bu(objhts))
    }
  }else{
    matpred <-  objreg$Yhat %*% t(objlearn$P) %*% t(objhts$S)
  }
  as.matrix(matpred)
}


bu <- function(objhts){
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) # use sparse matrices here
  list(algo = "BU", P = P_BU)
}

mls <- function(objreg, objhts){
  S <- objhts$S
  P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(objreg$Y) %*% objreg$Yhat %*% solve(t(objreg$Yhat) %*% objreg$Yhat)
  list(algo = "LS", P = P_LS)
}

mint <- function(objhts, method = NULL, residuals = NULL, h = NULL){
  J <- Matrix(cbind(matrix(0, nrow = objhts$nbts, ncol = objhts$nts - objhts$nbts), diag(objhts$nbts)), sparse = TRUE)
  U <- Matrix(rbind(diag(objhts$nts - objhts$nbts), -t(objhts$A)), sparse = TRUE)
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) 

  if(is.null(h)){
    h <- 1  
  }
  R1 <- t(residuals[h, , ])
  
  if(is.null(method))
    method <- "diagonal"
  
  if(method == "diagonal"){
    # Diagonal matrix
    W <- Diagonal(x = vec_w(R1))
  }else if(method == "shrink"){
    # Shrunk matrix
    target_diagonal <- lowerD(R1)
    shrink_results <- shrink.estim(R1, target_diagonal)
    W <- shrink_results$shrink.cov
  }
 
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  P_MINT <- P_BU - C1 %*% MAT3
  
  
  list(algo = "MINT", P = P_MINT)
}
