new_learnreg <- function(objreg, objhts, objmethod, standardizeX = TRUE, centerY = TRUE){
  Y <- objreg$Y
  Yhat <- objreg$Yhat
  
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
  
  algo <- objmethod$algo
  config <- objmethod$config
  model <- NULL
  s <- NULL
  if(algo == "LASSO" || algo == "NEW"){
    #model <- do.call(glmnet, c(list(x = X, y = y), config$glmnet)) 
    if(algo == "LASSO" || algo == "NEW"){
      #mylambdas <- c(model$lambda, seq(tail(model$lambda, 1), 0,length.out = 10))
      #mylambdas <- model$lambda
      
      if(algo == "NEW"){
        n_before <- nrow(X)
      
        set_lambda1 <- objmethod$set_lambda1
        #err_lambda1 <- vector("list", length(set_lambda1))
        min_err <- Inf
        for(ilambda1 in seq_along(set_lambda1) ){  
        
            lambda1 <- set_lambda1[[ilambda1]]
          
        Xbis <- sqrt(lambda1) * kronecker(objhts$A, diag(objhts$nts))
        X_augm <- rbind(X, Xbis)
        n_after <- nrow(X_augm)
        
        k <- objhts$nts - objhts$nbts
        P_constraint <- Matrix(0, nrow = k, ncol = objhts$nts, sparse = T)
        P_constraint[cbind(seq(k), seq(k)) ] <- 1
        
        if(!is.null(objmethod$Ptowards)){
          P_constraint <- P_constraint - objhts$A %*% objmethod$Ptowards
        }
        P_constraint <- t(P_constraint)
        
        ybis <-  sqrt(lambda1) * fct_vec(P_constraint)
        y_augm <- rbind(y, ybis)
        
        model <- do.call(glmnet, c(list(x = X_augm, y = y_augm), config$glmnet)) 
        foldid <- config$cvglmnet$foldid
        folds <- unique(foldid)
        foldbinary <- sapply(seq_along(folds), function(ifold){
          fold <- folds[ifold]
          id <- which(foldid == fold)
          vec <- numeric(length(y))
          vec[id] <- 1
          vec[seq(n_before + 1, n_after)] <- 1
          vec
        })
        nfolds <- ncol(foldbinary)
        lambdas <- model$lambda
        #err_folds <- array(NA, c(n_before, length(lambdas)) )
        err_folds <- vector("list", nfolds)
        
        for(j in seq(nfolds)){
          idtrain <- which(foldbinary[, j] == 1)
          idvalid <- which(foldbinary[, j] == 0)
          model_cv <- do.call(glmnet, c(list(x = X_augm[idtrain, ], y = y_augm[idtrain, ]), config$glmnet, list(lambda = lambdas))) 
          pred <- predict(model_cv, X_augm[idvalid, ])
          #err_folds[idvalid, ]  <-  (pred - y[idvalid, ])^2
          err_folds[[j]] <- (pred - y_augm[idvalid, ])^2
        }
        #err_lambda1[[ilambda1]] <- err_folds
        
        res  <- t(sapply(err_folds, function(mat){
          apply(mat, 2, mean)
        }))
        err_lambdas2 <- apply(res, 2, mean)
        local_err <- min(err_lambdas2)
        
        print(local_err)
        
        if(local_err < min_err){
          min_err <- local_err
          model_save <- model
          lambda1_min <- lambda1
          err_folds_save <- err_folds
        }
          
      }
 
        err_folds <- err_folds_save
        #err_cv <- apply(err_folds, 2, mean)
        res_cv <- t(sapply(seq(length(lambdas)), function(ilambda){ 
          cv_mean_eachfold <- sapply(err_folds, function(mat){
            mean(mat[, ilambda])
          })  
          cv_sd <- sd(cv_mean_eachfold)/sqrt(nfolds)
          cv_mean <- mean(cv_mean_eachfold)
          list(cv_mean = cv_mean, cv_sdplus = cv_mean + cv_sd, cv_sdminus = cv_mean - cv_sd)
        }))
        # matplot(res_cv)
        
        cv_mean <- unlist(res_cv[, "cv_mean"])
        if(objmethod$selection == "1se"){
          cv_sdplus <- res_cv[, "cv_sdplus"]
          i <- which.min(cv_mean)
          lambda_min <- lambdas[which.min(cv_mean)]
          set_lambdas <- lambdas[which(cv_mean <= cv_sdplus[i] & lambdas > lambda_min)]
          stopifnot(length(set_lambdas)>0)
          s <- max(set_lambdas)
        }else if(objmethod$selection == "min"){
          s <- lambdas[which.min(cv_mean)]
        }else{
          stop("error in lambda selection !")
        }
        
        #browser()
        
      }else{
      
        do.parallel <- FALSE
        if(nb.cores.cv > 1){
          do.parallel <- TRUE
          registerDoMC(cores = nb.cores.cv)
        }
        model <- do.call(cv.glmnet, 
                         c(list(x = X, y = y), 
                           config$cvglmnet ,
                           list(parallel = do.parallel)))
        s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
        #browser()
        #browser()
      }
    }
  }else if(algo == "GGLASSO"){
    group1 <- rep(seq(objhts$nts), objhts$nbts)
    #model <- cv.gglasso(x = as.matrix(X), y = y, group=group1, pred.loss="L2", foldid = config$cvglmnet$foldid, intercept = FALSE)
    model <- do.call(cv.gglasso, c(list(x = as.matrix(X), y = y, group = group1, pred.loss = "L2"), config))
    s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
  }
  
  obj_return <- list(model = model, s = s, scaling_info = scaling_info, objmethod = objmethod, 
                     standardizeX = standardizeX, centerY = centerY)
} 

new_predtest <- function(objreg, objhts, objlearn = NULL){
  
  objmethod <- objlearn$objmethod
  algo <- objmethod$algo
  
  
  if(algo == "LASSO" || algo == "NEW" || algo == "GGLASSO"){  
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
      vec_ytile_test <- as.numeric(predict(objlearn$model, newx = X_test, type = "response", s = objlearn$s))
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
