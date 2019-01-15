rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
source("config_paths.R")
source("packages.R")
source("bights.R")
source("methods.R")
source("simulate.R")
source("simulate_large.R")
source("utils.R")
source("hts.R")
source("code.R")
library(methods)
library(doMC)
library(gglasso)

experiment <- "tourism"
#experiment <- "wikipedia"

do.joint <- FALSE
#do.joint <-  experiment != "wikipedia"

if(experiment == "tourism"){
  do.log <- FALSE
}else{
  do.log <- TRUE
}
do.scaling <- TRUE
do.cleaning <- TRUE
do.deseasonalization <- TRUE
do.loo <- FALSE
use.intercept <- FALSE
add.bias <- FALSE
nobiasforwho <- "agg"  # "bottom" "agg"
regularization <- "lasso" # ridge

print(do.deseasonalization)
print(do.cleaning)
print(add.bias)
print(regularization)

H <- 2


#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
config_main <- list(intercept = use.intercept, standardize = FALSE, alpha = ifelse(regularization == "lasso", 1, 0), 
                    thresh = 10^-6, nlambda = 50)
config_cvglmnet <- c(config_main, list(nfolds = 3))

#lambda_selection <- "1se"
lambda_selection <- "min"


mc.cores.basef <- 6 #30 # 20
mc.cores.methods <- 1 #4 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 6 #3
do.save <- TRUE


config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
                        param_fit_fct = list(seasonal = TRUE, ic = "aic",  max.p = 1, max.q = 0, max.P = 1, max.Q = 0,
                                             approximation = TRUE, stationary = FALSE, num.cores = 2), 
                        param_forecast_fct = list(use.initial.values = TRUE))

config_forecast_agg <- config_forecast_bot <- config_forecast

if(experiment == "tourism"){
  X <- read.csv("../data/Tourism data_v3.csv")
  Z <- X[, -seq(2)]
  obj_sort <- sort(apply(Z == 0, 2, sum), index = T, decreasing = T)
  Z <- Z[, -head(obj_sort$ix, 34)]
  if(TRUE){
    series_name <- substr(colnames(Z), 1, 3)
    series_name <-sapply(seq(length(series_name)), function(j){
      paste(series_name[j], sprintf("%04d", j), sep = "")
    })
    colnames(Z) <- series_name
    y <- hts(Z, characters = c(1, 1, 1, 4))
  }else{
    y <- gts(Z, characters = list(c(1, 1, 1), c(3)))
  }
  refit_step <- 12
}else if(experiment == "wikipedia"){
  source("wikipedia.R")
  
  set.outliers <- c(185L, 277L, 4L, 182L, 65L, 77L, 138L, 128L, 30L, 274L)
  Z <- Z[, -set.outliers]
  
  y <- hts(Z, characters = c(2, 3, 3, 4) )
  refit_step <- 14
}

vec <- sapply(y$nodes, length)
niveaus <- unlist(sapply(seq(length(vec)), function(i){
  rep(i, vec[i])
}))
niveaus <- c(niveaus, rep(length(vec) + 1, ncol(y$bts)))

S <- smatrix(y)
A <- head(S, nrow(S) - ncol(Z))


POLS <- solve(t(S) %*% S) %*% t(S)

if(do.log){
  bts <- log(1 + Z)
}else{
  bts <- Z
}

if(do.cleaning){
  cleaned_bts <- sapply(seq(ncol(bts)), function(j){
    tsclean(bts[, j])
  })
  bts <- cleaned_bts
}

if(do.deseasonalization){
  yts <- t(S %*% t(bts))
  
  if(experiment == "tourism"){
    yts <- ts(yts, c(1998, 1), freq = 12)
  }else if(experiment == "wikipedia"){
    yts <- ts(yts, freq = 7)
  }
  
  zts <- sapply(seq(ncol(yts)), function(j){
    obj <- stl(yts[, j], s.window = "periodic", robust = TRUE)
    obj$time.series[, c("remainder")]
    #yts[, j] - obj$time.series[, c("seasonal")]
  })
  
  nbts <- ncol(bts)
  nyts <- ncol(yts)
  nats <- nyts - nbts
  
  zts_bottom <- zts[, seq(nats +1, nyts)]
  
  #P_OLS <- solve(t(S) %*% S) %*% t(S)
  #zts_bottom_tilde <- t(P_OLS %*% t(zts))
  
  zts_bottom_tilde <- zts_bottom
  bts <- zts_bottom_tilde
}

if(do.scaling){
  bts <- apply(bts, 2, scale, center = T, scale = T)
}

datasave_file <- file.path(results.folder, paste("datasave_", experiment, "_", 
                                                 do.log, "_", do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))
save(file = datasave_file, list = c("Z", "bts", "cleaned_bts"))


if(experiment == "tourism"){
  bts <- ts(bts, c(1998, 1), freq = 12)
}else if(experiment == "wikipedia"){
  bts <- ts(bts, freq = 7)
}


if(experiment == "tourism"){
  # PREVIOUSLY
  #T_train <- 9 * 12 #108 # 96 #7
  #T_valid <- 5 * 12 #60 # 5
  
  T_train <- 5 * 12
  T_valid <- 9 * 12
  T_test <- 60
}else if(experiment == "wikipedia"){
  T_train <- 86
  T_valid <- 160
  T_test <- 120
}

T_learn <- T_train + T_valid
T_all <- T_learn + T_test
stopifnot(T_all == nrow(bts))

my_bights <- bights(bts, A)
info_file <- file.path(results.folder, paste("info_", experiment, ".Rdata", sep = ""))
save(file = info_file, list = c("A", "niveaus"))
 
  
  print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))
  print(paste("Tvalid = ", T_valid, sep = ""))
  print(paste("N = ", my_bights$nts * T_valid, " - p = ", my_bights$nbts * my_bights$nts, sep = ""))
  

  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, "_", refit_step, "_", do.log, "_", do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))  
 
   
  if(do.save && file.exists(file_bf)){
    load(file_bf)
    print("loading")
  }else{
    print(base::date())
    print("base forecasting ...")
    ### valid
    list_subsets_valid <- lapply(seq(T_train, T_learn - H), function(i){c(i - T_train + 1, i)})
    data_valid <- makeMatrices(my_bights, list_subsets_valid, H = H, 
                               config_forecast_agg = config_forecast_agg, config_forecast_bot = config_forecast_bot, 
                               refit_step = refit_step, mc.cores = mc.cores.basef)
    ### test
    list_subsets_test <- lapply(seq(T_learn, T_all - H), function(i){c(i - T_learn + 1, i)}) # I CHANGED it from T_TRAIN TO T_LEARN !!
    data_test <- makeMatrices(my_bights, list_subsets_test, H = H, 
                              config_forecast_agg = config_forecast_agg, config_forecast_bot = config_forecast_bot, 
                              refit_step = refit_step, mc.cores = mc.cores.basef)
    if(do.save){
      save(file = file_bf, list = c("data_valid", "data_test", "list_subsets_valid", "list_subsets_test"))
    }
    print(base::date())
  }
  
  ### valid
  Yhat_valid_allh <- data_valid$Yhat
  Y_valid_allh     <- data_valid$Y
  
  ### test
  Yhat_test_allh  <- data_test$Yhat
  Y_test_allh    <- data_test$Y
  
  Eresiduals <- data_test$Eresiduals

  save_Yhat_test_allh <- Yhat_test_allh
  
  
  if(add.bias){
    nvalid <- dim(Yhat_valid_allh)[3]
    ntest  <- dim(Yhat_test_allh)[3]
    ninsample <- nrow(data_test$IN_y)
    
    n <- my_bights$nts
    #u <- runif(2 * n, 0.8, 0.95)
    u <- runif(2 * n, 0.05, 0.1)
    
    vec <- seq(1, 2*n, 2)
    
    ubias <- lapply(vec, function(j){
      myinterval <- sort(u[j:(j+1)])
      delta_valid <- runif(nvalid, myinterval[1], myinterval[2])
      delta_test <- runif(ntest, myinterval[1], myinterval[2])
      #delta_insample <- runif(ninsample, myinterval[1], myinterval[2])
      list(delta_valid = delta_valid, delta_test = delta_test)
    })

    delta_valid <- simplify2array(lapply(ubias, "[[", "delta_valid"))
    delta_test <- simplify2array(lapply(ubias, "[[", "delta_test"))
    #delta_insample <- simplify2array(lapply(ubias, "[[", "delta_insample"))
    
    if(nobiasforwho == "bottom"){
      idb <- seq(my_bights$naggts + 1, my_bights$naggts + my_bights$nbts)
      
      delta_valid[, idb] <- 1
      delta_test[, idb]  <- 1
    }else if(nobiasforwho == "agg"){
      idb <- seq(my_bights$naggts)  
      
      delta_valid[, idb] <- 1
      delta_test[, idb]  <- 1
    }
    

    
    for(h in seq(H)){
      Yhat_valid_allh[h, , ] <- Yhat_valid_allh[h, , ] * t(delta_valid)
      Yhat_test_allh[h, , ]  <- Yhat_test_allh[h, , ]  * t(delta_test)
    }
    
    #Eresiduals <- data_test$IN_y - (data_test$IN_yhat * delta_insample)
    
  }
  
  
  config <- list(glmnet = config_main, cvglmnet = config_cvglmnet)
  objmethod <- list(selection = lambda_selection)
  
  
  
  P_BU <- Matrix(0, nrow = my_bights$nbts, ncol = my_bights$nts, sparse = T)
  P_BU[cbind(seq(my_bights$nbts), seq(my_bights$nts - my_bights$nbts + 1, my_bights$nts)) ] <- 1
  P_OLS <- solve(t(S) %*% S) %*% t(S)
  P_0 <- Matrix(0, nrow = my_bights$nbts, ncol = my_bights$nts, sparse = T)
  
  
  wmethod <- substr("MINTshr", 5, 7)
  cov_method <-  switch(wmethod, ols = "ols", shr = "shrink", sam = "sample", NA)
  e_residuals <- switch(wmethod, ols = NULL, shr = Eresiduals, sam = Eresiduals, NA)
  obj_learn <- mint(my_bights, method = cov_method, e_residuals = e_residuals , h = 1)
  P_MINTSHRINK <- obj_learn$P
  P_MINTOLS <- P_OLS
  

  #P_REF <- P_0
  P_REF <- P_BU
  #P_REF <- P_OLS
  #P_REF  <- P_MINT
  
  Ytilde_test_allh <- array(NA, dim(Yhat_test_allh))
  results_allh <- vector("list", H)
  h <- 1
  Yhat_valid_h <- t(Yhat_valid_allh[h, , ])
  Y_valid_h    <- t(Y_valid_allh[h, , ])
  Yhat_test_h  <- t(Yhat_test_allh[h, , ])
  save_Yhat_test_h <- t(save_Yhat_test_allh[h, , ])

  
  nValid <- nrow(Yhat_valid_h)
  if(do.loo){
    config$cvglmnet$nfolds <- nValid
    config$cvglmnet$foldid <- NULL
  }else{
    nfolds <- config$cvglmnet$nfolds
    nbperfold <- floor(nValid/nfolds)
    foldid    <- rep(seq(nfolds), each = nbperfold)
    remainder <- nValid %% nfolds
    foldid    <- c(foldid, rep(nfolds, remainder))
    config$cvglmnet$foldid <- foldid
  }
  
  pred_reg <- matrix(NA, nrow = nrow(Yhat_test_h), ncol = my_bights$nbts)
  variables <- vector("list", my_bights$nbts)
  for(j in seq(my_bights$nbts)){
    if(j%%10 == 0)
    print(j)
    
    
    k <- my_bights$naggts + j
    X <- Yhat_valid_h
    
    B_towards_valid <- X %*% t(P_REF)
    y <- Y_valid_h[, k] - B_towards_valid[, j]

    
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
    s <- ifelse(objmethod$selection == "min", "lambda.min", "lambda.1se")
    
    variables[[j]] <- coef(model, s = s)
    #browser()
    #ypred <- predict(model, X, s = s)
    #matplot(cbind(ypred, y), lty = 1, type = 'l', col = c("purple", "black")); lines(Yhat_valid_h[, k], col = "blue")
    
    do.varsel <- FALSE
    if(do.varsel){
      tmp_coeffs <- coef(model, s = "lambda.min")[-1]
      idvar <- which(as.numeric(tmp_coeffs) != 0)
      B_towards_test <- Yhat_test_h %*% t(P_REF)
      
      if(length(idvar) == 0){
        pred_reg[, j] <- 0 + B_towards_test[, j] # 0 INSTEAD OF THE INTERCEPT
      }else{
        newXvalid <- X[, idvar, drop = F]
        colnm <- paste("var", seq(ncol(newXvalid)), sep = "")
        colnames(newXvalid) <- colnm
        
        model <- lm("y ~ .", data.frame(y , newXvalid))

        X_test <- data.frame(Yhat_test_h[, idvar, drop = F])
        colnames(X_test) <- colnm
        pred_reg[, j] <- predict(model, X_test) + B_towards_test[, j]
      }
     
    }else{
      X_test <- Yhat_test_h
      B_towards_test <- X_test %*% t(P_REF)
      pred_reg[, j] <- predict(model, X_test, s = s) + B_towards_test[, j]
      
      # NO INTERCEPT IF ZERO MODEL
      #tmp_coeffs <- coef(model, s = s)[-1]
      #nb_nonzero <- length(which(as.numeric(tmp_coeffs) != 0))
      #if(nb_nonzero == 0){
      #  pred_reg[, j] <- 0 + B_towards_test[, j]
      #}
      
    }
  }
  
  C_REG_withbias <- t(sapply(variables, function(p){ as.vector(p) }))
  P_REG <- P_REF + C_REG_withbias[, -1] 

  pred_mintshrink <- t(P_MINTSHRINK %*% t(Yhat_test_h))
  pred_mintols <- t(P_MINTOLS %*% t(Yhat_test_h))
  pred_bu   <- t(P_BU %*% t(Yhat_test_h)) 
  pred_base <-  save_Yhat_test_h
  pred_basebiased <- Yhat_test_h
  
  
  #### JOINT METHOD
  if(do.joint){
    print("doing joint ...")
    print(base::date())
    objreg_valid <- list(y = makey(Y_valid_h), X = makeX(Yhat_valid_h, my_bights), Y = Y_valid_h, Yhat = Yhat_valid_h)
    objreg_test <- list(X = makeX(Yhat_test_h, my_bights), Yhat = Yhat_test_h)
    nValid <- nrow(Yhat_valid_h)
    foldid <- make_foldid(nValid, my_bights$nts, config$cvglmnet$nfolds)
    config$cvglmnet$foldid <- foldid
    if(regularization == "ridge"){
      config_ridge <- config
      config_ridge$glmnet$alpha <- 0
      config_ridge$cvglmnet$alpha <- 0
      config <- config_ridge
    }
    obj_method <- list(algo = "LASSO", Ptowards = P_REF, config = config, selection = lambda_selection)
    obj_learn <- new_learnreg(objreg_valid, my_bights, obj_method, TRUE, TRUE)
    pred_joint <- new_predtest(objreg_test, my_bights, obj_learn)
    print("joint done...")
    print(base::date())
  }else{
    pred_joint <- pred_base
  }
  

  myfile <- file.path(results.folder, paste("resultsicml_", experiment, "_", lambda_selection, "_", 
                                            add.bias, "_", do.log, "_", 
                                            do.deseasonalization, "_", do.scaling, ".Rdata", sep = ""))
  save(file = myfile, list = c("P_REF", "P_REG", "C_REG_withbias", "P_BU", "P_MINTSHRINK", "P_MINTOLS", 
                               "variables", "pred_reg", "pred_mintshrink", "pred_mintols", "pred_bu", 
                               "pred_base", "pred_basebiased", "pred_joint", "Y_test_allh"))

