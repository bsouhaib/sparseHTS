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

do.deseasonalization <- TRUE
do.cleaning <- TRUE
add.bias <- TRUE
nobiasforwho <- "agg"  # "bottom" "agg"
regularization <- "lasso" # ridge
#regularization <- "ridge"

print(do.deseasonalization)
print(do.cleaning)
print(add.bias)
print(regularization)

do.conditioning <- FALSE
#Tyear <- rep(seq(12), 19)
#Tyear <- rep( c(1, rep(0, 11)), 19)
Tyear <- seq(12*19)

H <- 2





# INCLUDE INTERCEPT ?????
# P_REF ?????#
#stop("PUT INTERCEPT AT TRUE")
#stop("PUT standardize AT TRUE")


#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
config_main <- list(intercept = TRUE, standardize = FALSE, alpha = ifelse(regularization == "lasso", 1, 0), 
                    thresh = 10^-6, nlambda = 50)
config_cvglmnet <- c(config_main, list(nfolds = 6))

#lambda_selection <- "1se"
lambda_selection <- "min"
experiment <- "tourism"

mc.cores.basef <- 6 #30 # 20
mc.cores.methods <- 1 #4 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 6 #3
do.save <- TRUE
refit_step <- 12

config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
                        param_fit_fct = list(seasonal = TRUE, ic = "aic",  max.p = 1, max.q = 0, max.P = 1, max.Q = 0,
                                             approximation = TRUE, stationary = FALSE, num.cores = 2), 
                        param_forecast_fct = list(use.initial.values = TRUE))

config_forecast_agg <- config_forecast_bot <- config_forecast

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
S <- smatrix(y)
A <- head(S, nrow(S) - ncol(Z))


POLS <- solve(t(S) %*% S) %*% t(S)

bts <- Z



if(do.deseasonalization){
  yts <- t(S %*% t(bts))
  yts <- ts(yts, c(1998, 1), freq = 12)
  
  zts <- sapply(seq(ncol(yts)), function(j){
    obj <- stl(yts[, j], s.window = "periodic", robust = TRUE)
    obj$time.series[, c("remainder")]
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

if(do.cleaning){
  cleaned_bts <- sapply(seq(ncol(bts)), function(j){
    tsclean(bts[, j])
  })
  bts <- cleaned_bts
}

bts <- ts(bts, c(1998, 1), freq = 12)

# PREVIOUSLY
#T_train <- 9 * 12 #108 # 96 #7
#T_valid <- 5 * 12 #60 # 5

T_train <- 5 * 12
T_valid <- 9 * 12

T_learn <- T_train + T_valid
T_test <- 60
T_all <- T_learn + T_test
stopifnot(T_all == nrow(bts))

my_bights <- bights(bts, A)
info_file <- file.path(results.folder, paste("info_", experiment, "_", lambda_selection, ".Rdata", sep = ""))
save(file = info_file, list = c("A"))
 
  
  print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))
  print(paste("Tvalid = ", T_valid, sep = ""))
  print(paste("N = ", my_bights$nts * T_valid, " - p = ", my_bights$nbts * my_bights$nts, sep = ""))
  

  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, "_", refit_step, "_", do.deseasonalization, ".Rdata", sep = ""))  
 
   
  if(do.save && file.exists(file_bf)){
    load(file_bf)
    print("loading")
  }else{
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
    u <- runif(2 * n, 0.4, 0.6)
    
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
  
  if(do.conditioning){
    P_BU <- cbind(P_BU, 0)
    P_OLS <- cbind(P_OLS, 0)
    P_0 <- cbind(P_0, 0)
  }
  
  wmethod <- substr("MINTshr", 5, 7)
  cov_method <-  switch(wmethod, ols = "ols", shr = "shrink", sam = "sample", NA)
  e_residuals <- switch(wmethod, ols = NULL, shr = Eresiduals, sam = Eresiduals, NA)
  obj_learn <- mint(my_bights, method = cov_method, e_residuals = e_residuals , h = 1)
  P_MINTSHRINK <- obj_learn$P
  P_MINTOLS <- P_OLS
  
  if(do.conditioning){
    P_MINTSHRINK <- cbind(P_MINTSHRINK, 0)
    P_MINTOLS    <- cbind(P_MINTOLS, 0)
  }
  
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
  
  if(do.conditioning){
    idvalid <- seq(T_train + 1, T_learn - H + 1)
    idtest <- seq(T_learn + 1, T_learn + T_test - H + 1)

    Yhat_valid_h <- cbind(Yhat_valid_h, Tyear[idvalid])
    Yhat_test_h  <- cbind(Yhat_test_h, Tyear[idtest])
    
    penalty.factor	<- rep(1, ncol(Yhat_valid_h))
    penalty.factor[ncol(Yhat_valid_h)] <- 0

    config$glmnet <- c(config$glmnet, list(penalty.factor = penalty.factor))
    config$cvglmnet <- c(config$cvglmnet, list(penalty.factor = penalty.factor))
    
  }
  
  pred_reg <- matrix(NA, nrow = nrow(Yhat_test_h), ncol = my_bights$nbts)
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

  pred_mintshrink <- t(P_MINTSHRINK %*% t(Yhat_test_h))
  pred_mintols <- t(P_MINTOLS %*% t(Yhat_test_h))
  pred_bu   <- t(P_BU %*% t(Yhat_test_h)) 
  pred_base <-  save_Yhat_test_h
  pred_basebiased <- Yhat_test_h
  
  
  #### JOINT METHOD
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
  
  

  myfile <- file.path(results.folder, paste("resultsINDEP_", experiment, "_", 
                                            lambda_selection, "_", add.bias, "_", do.deseasonalization, ".Rdata", sep = ""))
  save(file = myfile, list = c("pred_reg", "pred_mintshrink", "pred_mintols", "pred_bu", 
                               "pred_base", "pred_basebiased", "pred_joint", "Y_test_allh"))


  
  
  
  
