

new_learnreg <- function(objreg, objhts, objmethod, standardizeX = NULL, centerY = NULL){
  Y <- objreg$Y
  Yhat <- objreg$Yhat
  algo <- objmethod$algo
  config <- objmethod$config
  
  if(algo != "ERMreg" && algo != "MRCE"){
    # Penalize towards a certain P matrix?
    if(!is.null(objmethod$Ptowards)){
      Y <- as.matrix(Y - Yhat %*% t(objmethod$Ptowards) %*% t(objhts$S))
    }
    
    scaling_info <- NULL
    if(centerY){
      #
      mu_Y <- colMeans(Y)
      Y_centered <- t(t(Y) - mu_Y)
      y <- fct_vec(Y_centered)
      scaling_info$mu_Y <- mu_Y
    }else{
      y <- fct_vec(Y)
    }
    
    if(standardizeX){
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
      
      scaling_info$mu_Yhat <- mu_Yhat
      scaling_info$sd_Yhat <- sd_Yhat
      
      X <- makeX(Yhat_scaled, objhts)
    }else{
      X <- makeX(Yhat, objhts)
    }
    
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
                           list(parallel = do.parallel) ))
                           #,
                           #list(penalty.factor = weights_shrinkage)  ))
        s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
        #browser()
        
      }
    }else if(algo == "GGLASSO"){
      group1 <- rep(seq(objhts$nts), objhts$nbts)
      #model <- cv.gglasso(x = as.matrix(X), y = y, group=group1, pred.loss="L2", foldid = config$cvglmnet$foldid, intercept = FALSE)
      model <- do.call(cv.gglasso, c(list(x = as.matrix(X), y = y, group = group1, pred.loss = "L2"), config))
      s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
    }
    obj_return <- list(model = model, s = s, scaling_info = scaling_info, objmethod = objmethod, 
                       standardizeX = standardizeX, centerY = centerY)
  }else if(algo == "MRCE"){
    
    S <- objhts$S
    B <- Y[, -seq(objhts$naggts)]
    omegaS <- t(S) %*% S
    lam2.vec <- rev(10^seq(from=-2, to=0, by=0.5))
    
    browser()
    
    fit=mrce(X = Yhat, Y = B, lam1=10^(-1.5), lam2=10^(-0.5), method="single")
    fit
    lam2.mat=1000*(fit$Bhat==0)
    refit=mrce(X = Yhat, Y = B, lam2=lam2.mat, method="fixed.omega", omega=omegaS, tol.in=1e-12) 
    refit
    
    # ce code renvoit un NA/NaN/Inf in foreign function call 
    res <- mrce(X = Yhat, Y = B, method = "fixed.omega", omega = omegaS, silent  = F, lam2.vec=lam2.vec)
    
    
  }else{
    s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
    models <- variables <- vector("list", objhts$nbts)
    P_REF <- objmethod$Ptowards
    
    for(j in seq(objhts$nbts)){
      if(j%%10 == 0)
        print(j)
      
      #if(j == 261)
      #  browser()
      
      k <- objhts$naggts + j
      Bk <- Y[, k]
      X <- Yhat

      B_towards_valid <- X %*% t(P_REF)
      y <- Bk - B_towards_valid[, j]
      
      model <- do.call(glmnet, c(list(x = X, y = y), config$glmnet)) 
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
                         list(parallel = do.parallel) ))
      models[[j]] <- model
      variables[[j]] <- coef(model, s = s)
    }
    C_REG_withbias <- t(sapply(variables, function(p){ as.vector(p) }))
    Pmatrix <- P_REF + C_REG_withbias[, -1] 
    obj_return <- list(models = models, s = s, objmethod = objmethod, variables = variables, P = Pmatrix)
  }
  
  obj_return
} 

new_predtest <- function(objreg, objhts, objlearn = NULL, svalue = NULL){
  
  objmethod <- objlearn$objmethod
  algo <- objmethod$algo
  
  
  if(algo == "LASSO" || algo == "gOLS" || algo == "GGLASSO"){  
    Yhat <- objreg$Yhat
    standardizeX <- objlearn$standardizeX
    centerY      <- objlearn$centerY
    
    if(standardizeX){
      scaling_info <- objlearn$scaling_info
      Yhat_scaled <- t((t(Yhat) - scaling_info$mu_Yhat)/scaling_info$sd_Yhat)
      
      X_test <- makeX(Yhat_scaled, my_bights)
      #X_test <- as.matrix(X_test) # why not keep it sparse?
    }else{
      X_test <- makeX(Yhat, my_bights)
    }
    
    ##
    if(algo == "GGLASSO"){
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, s = objlearn$s))
    }else{
      if(is.null(svalue)){
        vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = objlearn$s))
      }else{
        vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = svalue))
      }
    }
    
    Ytilde_test <- fct_inv_vec(vec_ytile_test, 
                nrows = nrow(Yhat), ncolumns = ncol(Yhat)) 
    if(centerY){
      Ytilde_test <- t((t(Ytilde_test) + scaling_info$mu_Y))
    }
    
    if(!is.null(objmethod$Ptowards)){
      Ytilde_test <- Ytilde_test + 
        new_predtest(objreg, objhts, list(objmethod = list(algo = "BU"), P = objmethod$Ptowards))
    }
    
    #Ytilde_centered_test <- fct_inv_vec(vec_ytile_test, 
    #                                    nrows = nrow(Yhat_scaled), ncolumns = ncol(Yhat_scaled)) 
    #Ytilde_uncentered_test <- t((t(Ytilde_centered_test) + scaling_info$mu_Y))
    #if(!is.null(objmethod$Ptowards)){
    #  Ytilde_test <- Ytilde_uncentered_test + 
    #    new_predtest(objreg, objhts, list(objmethod = list(algo = "BU"), P = objmethod$Ptowards))
    #}else{
    #  Ytilde_test <- Ytilde_uncentered_test
    #}

  }else if(algo %in% c("BU", "MINT", "ERM") ){
    if(algo == "ERM"){
      X <- objreg$Yhat
      if(!is.null(objlearn$scale_info)){
        X_center <- objlearn$scale_info$X_center
        X_scale <- objlearn$scale_info$X_scale
        B_center <- objlearn$scale_info$B_center
        X <- t((t(X) - X_center)/X_scale)
      }
      
      Btilde_test <- X %*% t(objlearn$P)
      if(!is.null(objlearn$scale_info)){
        Btilde_test <- t(t(Btilde_test) +  B_center)
      }
      Ytilde_test <-  Btilde_test %*% t(objhts$S)
    }else{
      Ytilde_test <-  objreg$Yhat %*% t(objlearn$P) %*% t(objhts$S)
    }

  }else if(algo %in% c("OLSS")){
      Yhat <- objreg$Yhat
      scaling_info <- objlearn$scaling_info
      Yhat_scaled <- t((t(Yhat) - scaling_info$mu_Yhat)/scaling_info$sd_Yhat)
      Ytilde_centered_test <-  Yhat_scaled %*% t(objlearn$P) %*% t(objhts$S)
      Ytilde_uncentered_test <- t((t(Ytilde_centered_test) + scaling_info$mu_Y))
      Ytilde_test <- Ytilde_uncentered_test
  }else if(algo == "ERMreg"){
    P_REF <- objmethod$Ptowards
    X_test <- objreg$Yhat
    B_towards_test <- X_test %*% t(P_REF)
    Btilde_test <- matrix(NA, nrow = nrow(X_test), ncol = objhts$nbts)
    for(j in seq(objhts$nbts)){
      Btilde_test[, j] <- predict(objlearn$models[[j]], X_test, s = objlearn$s) + B_towards_test[, j]
    }
    
    Ytilde_test <- Btilde_test %*% t(objhts$S)
    
  }else{
    stop("ERROR IN METHOD'S NAME")
  }
  
  as.matrix(Ytilde_test)
}


bu <- function(objhts){
  P_BU <- cbind(matrix(0, objhts$nbts, objhts$nts - objhts$nbts), diag(objhts$nbts)) # use sparse matrices here
  list(objmethod = list(algo = "BU"), P = P_BU)
}


erm <- function(objreg, objhts){
  #S <- objhts$S
  #P_LS <- solve(t(S) %*% S) %*% t(S) %*% t(objreg$Y) %*% objreg$Yhat %*% solve(t(objreg$Yhat) %*% objreg$Yhat)
  #list(objmethod = list(algo = "OLS"), P = P_LS)
  

  S <- objhts$S
  B <- objreg$Y[, seq(objhts$naggts + 1, objhts$nts)]
  X <- objreg$Yhat
  
  scale_info <- NULL
  do.scaling <- TRUE
  if(do.scaling){
    Bscaled <- scale(B, center = T, scale = F)
    Xscaled <- scale(X, center = T, scale = T)
    B_center <- attr(Bscaled, "scaled:center")
    X_center <- attr(Xscaled, "scaled:center")
    X_scale <- attr(Xscaled, "scaled:scale")
    if(any(X_scale == 0) ){
      id <- which(X_scale == 0)
      X_scale[id] <- 1
      Xscaled[, id] <- scale(X[, id, drop = F], center = T, scale = F)
      X_center[id] <- apply(X[, id, drop = F], 2, mean)
    }
    scale_info <- list(B_center = B_center, X_center = X_center, X_scale = X_scale)
    X <- Xscaled
    B <- Bscaled
  }
  #browser()
  P_ERM <- t(B) %*% X %*% ginv(t(X) %*% X)
  #P_ERM <- t(B) %*% X %*% ginv(t(X) %*% X + 0.1 * diag(ncol(X)))
  list(objmethod = list(algo = "ERM"), P = P_ERM, scale_info = scale_info)
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
    #if(is.positive.definite(W)==FALSE)
    #{
    #  stop("MinT needs covariance matrix to be positive definite.", call. = FALSE)
    #}
  }
 
  MAT1 <- W %*% U
  MAT2 <- crossprod(U,MAT1)
  MAT3 <- tcrossprod(ginv(as.matrix(MAT2)), U)
  # MAT3 <- tcrossprod(solve(MAT2), U)
  C1 <- J %*% MAT1
  P_MINT <- P_BU - C1 %*% MAT3
  
  
  list(objmethod = list(algo = "MINT"), P = P_MINT)
}
