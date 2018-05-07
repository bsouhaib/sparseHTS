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


#######
#name_methods <- c("L1-PBU", "L1", "BU", "MINTols", "MINTshr", "BASE","BASE2", 
#                  "L2", "L2-PBU")
name_methods <- c("L1-PBU", "BU", "MINTols", "MINTshr", "BASE","BASE2", "L2-PBU")
#, "G-L1-PBU"
nb_methods <- length(name_methods)
#######

lambda_selection <- "1se"
#lambda_selection <- "min"
experiment <- "tourism"
H <- 1

X <- read.csv("../data/Tourism data_v3.csv")
Z <- X[, -seq(2)]
vec <- colnames(Z)
y <- gts(Z, characters = list(c(1, 1, 1), c(3)))
S <- smatrix(y)
ally <- allts(y)
bts <- Z
bts <- ts(bts, c(1998, 1), freq = 12)
A <- head(S, nrow(S) - ncol(Z))

# CONSIDER ETS ????
config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
                        param_fit_fct = list(seasonal = TRUE, ic = "aic",  max.p = 2, max.q = 2,
                                             approximation = TRUE, stationary = FALSE, num.cores = 2), 
                        param_forecast_fct = list(use.initial.values = TRUE))
config_forecast_agg <- config_forecast_bot <- config_forecast

T_train <- 96 # 7
T_valid <- 60 # 5
T_learn <- T_train + T_valid
T_test <- 72
T_all <- T_learn + T_test
stopifnot(T_all == nrow(bts))


mc.cores.basef <- 30 # 20
mc.cores.methods <- 7 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 3

do.save <- TRUE
refit_step <- 1
sameP_allhorizons <- TRUE
nb_simulations <- 1

for(i in seq(nb_simulations)){
  
  my_bights <- bights(bts, A)
  
  if(i == 1){
    info_file <- file.path(results.folder, paste("info_", experiment, "_", lambda_selection, ".Rdata", sep = ""))
    #save(file = info_file, list = c("name_methods", "A"))
    save(file = info_file, list = c("A"))
  }
  
  print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))
  print(paste("Tvalid = ", T_valid, sep = ""))
  print(paste("N = ", my_bights$nts * T_valid, " - p = ", my_bights$nbts * my_bights$nts, sep = ""))
  
  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, "_", i, "_", refit_step, ".Rdata", sep = ""))  
  
  if(do.save && file.exists(file_bf)){
    load(file_bf)
  }else{
    ### valid
    list_subsets_valid <- lapply(seq(T_train, T_learn - H), function(i){c(i - T_train + 1, i)})
    data_valid <- makeMatrices(my_bights, list_subsets_valid, H = H, 
                               config_forecast_agg = config_forecast_agg, config_forecast_bot = config_forecast_bot, 
                               refit_step = refit_step, mc.cores = mc.cores)
    ### test
    list_subsets_test <- lapply(seq(T_learn, T_all - H), function(i){c(i - T_learn + 1, i)}) # I CHANGED it from T_TRAIN TO T_LEARN !!
    data_test <- makeMatrices(my_bights, list_subsets_test, H = H, 
                              config_forecast_agg = config_forecast_agg, config_forecast_bot = config_forecast_bot, 
                              refit_step = refit_step, mc.cores = mc.cores)
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
  
  
  config_basic <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
  config <- list(glmnet = config_basic, cvglmnet = c(config_basic, list(nfolds = 3)))
  config_gglasso <- list(eps = 10^-6, nlambda = 50, intercept = FALSE, maxit = 1e+7)
  
  #sgl_config    <- list(nfold = 3, standardize = FALSE, alpha = .8)
  #print(glmnet_config)
  #print(sgl_config)
  
  # naive predictions
  #id <- sapply(list_subsets_test, function(vec){ vec[2] })
  #predictions_naive <- my_bights$yts[id, ]
  
  #predictions_avg <- t(sapply(list_subsets_test, function(vec){ 
  #  apply(my_bights$yts[vec[1]:vec[2], ], 2, mean)
  #}))
  
  Ytilde_test_allh <- array(NA, c(dim(Yhat_test_allh), nb_methods) )
  
  results_allh <- vector("list", H)
  for(h in seq(H)){
    print(paste("h = ", h, sep = ""))
    
    Yhat_valid_h <- t(Yhat_valid_allh[h, , ])
    Y_valid_h    <- t(Y_valid_allh[h, , ])
    Yhat_test_h  <- t(Yhat_test_allh[h, , ])
    
    objreg_valid <- list(y = makey(Y_valid_h), X = makeX(Yhat_valid_h, my_bights), Y = Y_valid_h, Yhat = Yhat_valid_h, 
                         Yhat_intercept = cbind(Yhat_valid_h, rep(1, nrow(Yhat_valid_h)) ) )
    objreg_test <- list(X = makeX(Yhat_test_h, my_bights), Yhat = Yhat_test_h,
                        Yhat_intercept = cbind(Yhat_test_h, rep(1, nrow(Yhat_test_h)) ) )
    
    nValid <- nrow(Yhat_valid_h)
    foldid <- rep(sample(seq(config$cvglmnet$nfolds), nValid, replace = T), my_bights$nts)
    #glmnet_config_local <- c(glmnet_config, list(foldid = foldid))
    #glmnet_config_local$nfolds <- NA
    config$cvglmnet$foldid <- foldid
    config$cvglmnet$nfolds <- NA
    config_gglasso$foldid <- foldid
    
    print("Starting methods")
    
    results <- mclapply(name_methods, function(current_method){
      #print(current_method)
      obj_method <- NULL
      if(current_method == "BU"){
        obj_learn <- bu(my_bights)
      }else if(grepl("MINT", current_method)){
        wmethod <- substr(current_method, 5, 7)
        cov_method <-  switch(wmethod, ols = "ols", shr = "shrink", sam = "sample", NA)
        e_residuals <- switch(wmethod, ols = NULL, shr = Eresiduals, sam = Eresiduals, NA)
        obj_learn <- mint(my_bights, method = cov_method, e_residuals = e_residuals , h = h) # J U and W
      }else if(current_method == "OLS"){
        obj_learn <- ols(objreg_valid, my_bights) 
      }else if(current_method == "L1"){
        obj_method <- list(algo = "LASSO", Ptowards = NULL, config = config, selection = lambda_selection)
      }else if(current_method == "L1-PBU"){
        obj_method <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
      }else if(current_method == "L2" || current_method == "L2-PBU"){
        config_ridge <- config
        config_ridge$glmnet$alpha <- 0
        config_ridge$cvglmnet$alpha <- 0 # !!!!!!!
        if(current_method == "L2"){
          obj_method <- list(algo = "LASSO", Ptowards = NULL, config = config_ridge, selection = lambda_selection)
        }else{
          obj_method <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config_ridge, selection = lambda_selection)
        }
      }else if(current_method == "G-L1-PBU"){
        obj_method <- list(algo = "GGLASSO", Ptowards = pbu(my_bights), config = config_gglasso, selection = lambda_selection)
      }else{
        obj_learn <- NULL
      }
      
      if(!is.null(obj_method)){
        if(!sameP_allhorizons || (sameP_allhorizons && h == 1)){
          obj_learn <- new_learnreg(objreg_valid, my_bights, obj_method)
        }else{
          obj_learn <- results[[current_method]]$obj_learn
        }
      }
      predictions <- NULL
      if(!is.null(obj_learn)){
        predictions <- new_predtest(objreg_test, my_bights, obj_learn)
      }
      
      print(current_method)
      
      list(obj_learn = obj_learn, predictions = predictions)
    }, mc.cores = mc.cores.methods)
    names(results) <- name_methods
    
    results[["BASE"]]   <- list(obj_learn = NULL, predictions = t(Yhat_test_allh[h, , ]))
    results[["BASE2"]]   <- list(obj_learn = NULL, predictions = t(save_Yhat_test_allh[h, , ]))
    
    results_allh[[h]] <- results
    
  } # HORIZON
  
  
  myfile <- file.path(results.folder, paste("results_", experiment, "_", i, "_", lambda_selection, ".Rdata", sep = ""))
  save(file = myfile, list = c("results_allh", "Y_test_allh"))
} # simulation


