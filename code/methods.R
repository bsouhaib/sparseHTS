new_learnreg <- function(objreg, objhts, objmethod){
  
  Y <- objreg$Y
  Yhat <- objreg$Yhat
  
  # Penalize towards a certain P matrix?
  if(!is.null(objmethod$Ptowards)){
    Y <- as.matrix(Y - Yhat %*% t(objmethod$Ptowards) %*% t(objhts$S))
  }
  

  #Y_centered <- scale(Y, center = TRUE, scale = FALSE)
  #center_Y <- attr(Y_centered, "scaled:center")
  
  #Yhat_scaled <- scale(Yhat, center = TRUE, scale = TRUE)
  #center_Yhat <- attr(Yhat_scaled, "scaled:center")
  #scale_Yhat  <- attr(Yhat_scaled, "scaled:scale")
  
  #
  mu_Y <- colMeans(Y)
  Y_centered <- t(t(Y) - mu_Y)
  
  mu_Yhat <- colMeans(Yhat)
  Yhat_centered <- t(t(Yhat) - mu_Yhat)
  sd_Yhat <- apply(Yhat_centered, 2, sd)
  
  # Dealing with constant variables
  id_constvar <- which(sd_Yhat == 0)
  if(length(id_constvar) > 0)
    sd_Yhat[id_constvar] <- 1
  
  Yhat_scaled <- t(t(Yhat_centered)/sd_Yhat)
  
  var_X_exclude <- NULL
  if(length(id_constvar) > 0){
    var_yhat_exclude <- rep(0, ncol(Yhat))
    var_yhat_exclude[id_constvar] <- 1
    var_X_exclude <- which(kronecker(rep(1, objhts$nbts), var_yhat_exclude) == 1)
    
    config$glmnet <- c(config$glmnet, list(exclude = var_X_exclude))
    config$cvglmnet <- c(config$cvglmnet, list(exclude = var_X_exclude))
  }
  
  scaling_info <- list(mu_Y = mu_Y, mu_Yhat = mu_Yhat, sd_Yhat = sd_Yhat)
    
  y <- fct_vec(Y_centered)
  X <- makeX(Yhat_scaled, objhts)

  algo <- objmethod$algo
  if(algo == "LASSO" || algo == "OLS"){
    config <- objmethod$config
    
    model <- do.call(glmnet, c(list(x = X, y = y), config$glmnet)) 
    if(algo == "OLS"){
      s <- 0
    }else if(algo == "LASSO"){
      mylambdas <- c(model$lambda, seq(tail(model$lambda, 1), 0,length.out = 25))
      model <- do.call(cv.glmnet, c(list(x = X, y = y), config$cvglmnet , list(lambda = mylambdas)))
      s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
      
      # SHOULD WE use the weight argument of glmnet?????
    }
  }else if(algo == "GGLASSO"){
    group1 <- rep(seq(objhts$nts), objhts$nbts)
    #browser()
    model <- cv.gglasso(x = as.matrix(X), y = y, group=group1, pred.loss="L2", foldid = config$cvglmnet$foldid, intercept = FALSE)
    s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
  }
  
  obj_return <- list(model = model, s = s, scaling_info = scaling_info, objmethod = objmethod)
} 

new_predtest <- function(objreg, objhts, objlearn = NULL){
  
  objmethod <- objlearn$objmethod
  algo <- objmethod$algo
  
  if(algo == "LASSO" || algo == "OLS" || algo == "GGLASSO"){  
    Yhat <- objreg$Yhat
    
    scaling_info <- objlearn$scaling_info
    Yhat_scaled <- t((t(Yhat) - scaling_info$mu_Yhat)/scaling_info$sd_Yhat)
    
    X_test <- makeX(Yhat_scaled, my_bights)
    X_test <- as.matrix(X_test) # why not keep it sparese?
    
    ##
    if(algo == "GGLASSO"){
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, s = objlearn$s))
    }else{
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = objlearn$s))
    }
    
    #Ytilde_centered_test <- makeY(vec_ytile_test, objhts)  # USE fct_inv_vec
    Ytilde_centered_test <- fct_inv_vec(vec_ytile_test, nrows = nrow(Yhat_scaled), ncolumns = ncol(Yhat_scaled)) 
    
    Ytilde_uncentered_test <- t((t(Ytilde_centered_test) + scaling_info$mu_Y))
    
    if(!is.null(objmethod$Ptowards)){
      Ytilde_test <- Ytilde_uncentered_test + 
        new_predtest(objreg, objhts, list(objmethod = list(algo = "BU"), P = objmethod$Ptowards))
    }else{
      Ytilde_test <- Ytilde_uncentered_test
    }
  }else if(algo %in% c("BU", "MINT") ){
    Ytilde_test <-  objreg$Yhat %*% t(objlearn$P) %*% t(objhts$S)
  }else{
    stop("ERROR IN METHOD'S NAME")
  }
  as.matrix(Ytilde_test)
}





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
  if(algo == "glmnet" || algo == "glmnetOLS"){
    
    if(algo == "glmnetOLS"){
      cvfit <- model <- do.call(glmnet, c(list(x = X, y = y), config))
      s <- 0
      idmin_lambda <- NULL
    }else{
      model <- do.call(glmnet, c(list(x = X, y = y), config[!names(config) %in% c("foldid", "nfolds") ]))
      mylambdas <- c(model$lambda, seq(tail(model$lambda, 1), 0,length.out = 25))
      cvfit <- model <- do.call(cv.glmnet, c(list(x = X, y = y), config, list(lambda = mylambdas)))
      s <- ifelse(selection == "min", "lambda.min", "lambda.1se")
      idmin_lambda <-  ifelse(selection == "min", match(model$lambda.min, model$lambda), match(model$lambda.1se, model$lambda))
    }
    
    browser()
    #slm.model <- slm.fit(x = X, y = as.numeric(y) )
    #X2 <- new("matrix.csr", X2)
    
    beta_hat <- as.numeric(coef(cvfit, s = s))
    if(!config$intercept)
      beta_hat <- beta_hat[-1]
    number_zeroes <- length(which(abs(beta_hat) < 10^-16))
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
    #beta_hat <- fit_SGL$beta[, idmin_lambda]
    lambda_final <- fit_SGL$lambdas[idmin_lambda]
  }
  
  obj_return <- list(algo = algo, model = model, idmin_lambda = idmin_lambda, s = s, towards_pbu = towards_pbu,
                     number_zeroes = number_zeroes)
  
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
  
  if(algo == "glmnet" || algo == "glmnetOLS" || grepl("SGL", algo)){  
    if(algo == "glmnet" || algo == "glmnetOLS"){
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
  }else if(objlearn$algo %in% c("BU", "MINT") ){
    matpred <-  objreg$Yhat %*% t(objlearn$P) %*% t(objhts$S)
  }else if(objlearn$algo == "LS"){
    # INTERCEPT IS INCLUDED HERE
    matpred <-  objreg$Yhat_intercept %*% t(objlearn$P) %*% t(objhts$S)
  }else{
    stop("ERROR IN METHOD'S NAME")
  }
  as.matrix(matpred)
}


bu <- function(objhts){
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) # use sparse matrices here
  list(objmethod = list(algo = "BU"), P = P_BU)
}

mls <- function(objreg, objhts){
  S <- objhts$S
  #P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(objreg$Y) %*% objreg$Yhat %*% solve(t(objreg$Yhat) %*% objreg$Yhat)
  P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(objreg$Y) %*% objreg$Yhat_intercept %*% solve(t(objreg$Yhat_intercept) %*% objreg$Yhat_intercept)
  list(objmethod = list(algo = "LS"), P = P_LS)
}

mint <- function(objhts, method = NULL, residuals = NULL, h = NULL){
  J <- Matrix(cbind(matrix(0, nrow = objhts$nbts, ncol = objhts$nts - objhts$nbts), diag(objhts$nbts)), sparse = TRUE)
  U <- Matrix(rbind(diag(objhts$nts - objhts$nbts), -t(objhts$A)), sparse = TRUE)
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) 

  if(is.null(h)){
    h <- 1  
  }
  
  if(!is.null(residuals))
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
  }else if(method == "ols"){
    W <- diag(objhts$nts)
  }
 
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  P_MINT <- P_BU - C1 %*% MAT3
  
  
  list(objmethod = list(algo = "MINT"), P = P_MINT)
}
