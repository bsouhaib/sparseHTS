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
  config <- objmethod$config
  model <- NULL
  s <- NULL
  if(algo == "LASSO" || algo == "gOLS"){
    model <- do.call(glmnet, c(list(x = X, y = y), config$glmnet)) 
    if(algo == "gOLS"){
      s <- 0
    }else if(algo == "LASSO"){
      #mylambdas <- c(model$lambda, seq(tail(model$lambda, 1), 0,length.out = 10))
      mylambdas <- model$lambda
      
      do.parallel <- FALSE
      if(nb.cores.cv > 1){
        do.parallel <- TRUE
        registerDoMC(cores = nb.cores.cv)
      }
      model <- do.call(cv.glmnet, 
                       c(list(x = X, y = y), 
                         config$cvglmnet , 
                         list(lambda = mylambdas), 
                         list(parallel = do.parallel)))
      s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
      
      #browser()
    }
  }else if(algo == "GGLASSO"){
    group1 <- rep(seq(objhts$nts), objhts$nbts)
    #model <- cv.gglasso(x = as.matrix(X), y = y, group=group1, pred.loss="L2", foldid = config$cvglmnet$foldid, intercept = FALSE)
    model <- do.call(cv.gglasso, c(list(x = as.matrix(X), y = y, group = group1, pred.loss = "L2"), config))
    s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
  }
  
  obj_return <- list(model = model, s = s, scaling_info = scaling_info, objmethod = objmethod)
} 

new_predtest <- function(objreg, objhts, objlearn = NULL){
  
  objmethod <- objlearn$objmethod
  algo <- objmethod$algo
  
  if(algo == "LASSO" || algo == "gOLS" || algo == "GGLASSO"){  
    Yhat <- objreg$Yhat
    
    scaling_info <- objlearn$scaling_info
    Yhat_scaled <- t((t(Yhat) - scaling_info$mu_Yhat)/scaling_info$sd_Yhat)
    
    X_test <- makeX(Yhat_scaled, my_bights)
    #X_test <- as.matrix(X_test) # why not keep it sparese?
    
    ##
    if(algo == "GGLASSO"){
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, s = objlearn$s))
    }else{
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = objlearn$s))
    }
    
    #Ytilde_centered_test <- makeY(vec_ytile_test, objhts)  # USE fct_inv_vec
    Ytilde_centered_test <- fct_inv_vec(vec_ytile_test, 
                                        nrows = nrow(Yhat_scaled), ncolumns = ncol(Yhat_scaled)) 

    Ytilde_uncentered_test <- t((t(Ytilde_centered_test) + scaling_info$mu_Y))
    if(!is.null(objmethod$Ptowards)){
      Ytilde_test <- Ytilde_uncentered_test + 
        new_predtest(objreg, objhts, list(objmethod = list(algo = "BU"), P = objmethod$Ptowards))
    }else{
      Ytilde_test <- Ytilde_uncentered_test
    }

  }else if(algo %in% c("BU", "MINT", "OLS") ){
    Ytilde_test <-  objreg$Yhat %*% t(objlearn$P) %*% t(objhts$S)
  }else if(algo %in% c("OLSS")){
      Yhat <- objreg$Yhat
      scaling_info <- objlearn$scaling_info
      Yhat_scaled <- t((t(Yhat) - scaling_info$mu_Yhat)/scaling_info$sd_Yhat)
      Ytilde_centered_test <-  Yhat_scaled %*% t(objlearn$P) %*% t(objhts$S)
      Ytilde_uncentered_test <- t((t(Ytilde_centered_test) + scaling_info$mu_Y))
      Ytilde_test <- Ytilde_uncentered_test
  }else{
    stop("ERROR IN METHOD'S NAME")
  }
  
  as.matrix(Ytilde_test)
}


bu <- function(objhts){
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) # use sparse matrices here
  list(objmethod = list(algo = "BU"), P = P_BU)
}


ols <- function(objreg, objhts){
  S <- objhts$S
  P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(objreg$Y) %*% objreg$Yhat %*% solve(t(objreg$Yhat) %*% objreg$Yhat)
  list(objmethod = list(algo = "OLS"), P = P_LS)
}

olss <- function(objreg, objhts){
  
  Y <- objreg$Y
  mu_Y <- colMeans(Y)
  Y_centered <- t(t(Y) - mu_Y)
  
  Yhat <- objreg$Yhat
  mu_Yhat <- colMeans(Yhat)
  Yhat_centered <- t(t(Yhat) - mu_Yhat)
  sd_Yhat <- apply(Yhat_centered, 2, sd)
  Yhat_scaled <- t(t(Yhat_centered)/sd_Yhat)
  
  S <- objhts$S
  P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(Y_centered) %*% Yhat_scaled %*% solve(t(Yhat_scaled) %*% Yhat_scaled)
  list(objmethod = list(algo = "OLSS"), P = P_LS)
}

mint <- function(objhts, method = NULL, e_residuals = NULL, h = NULL){
  J <- Matrix(cbind(matrix(0, nrow = objhts$nbts, ncol = objhts$nts - objhts$nbts), diag(objhts$nbts)), sparse = TRUE)
  U <- Matrix(rbind(diag(objhts$nts - objhts$nbts), -t(objhts$A)), sparse = TRUE)
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) 

  if(is.null(h)){
    h <- 1  
  }
  
  if(!is.null(e_residuals))
  R1 <- e_residuals
  
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
  }else if(method == "sample"){
    n <- nrow(R1)
    W <- crossprod(R1) / n
    if(is.positive.definite(W)==FALSE)
    {
      stop("MinT needs covariance matrix to be positive definite.", call. = FALSE)
    }
  }
 
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  P_MINT <- P_BU - C1 %*% MAT3
  
  
  list(objmethod = list(algo = "MINT"), P = P_MINT)
}
