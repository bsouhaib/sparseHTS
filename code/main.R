rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  
  experiment <- "small-unbiased"
  idjob <- 1
  nb_simulations <- 10
  lambda_selection <- "1se"
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  idjob <- as.numeric(args[[2]])
  nb_simulations <- as.numeric(args[[3]])
  lambda_selection <- args[[4]]
}

#set.seed(6120)
set.seed(idjob)

stopifnot(lambda_selection %in% c("min", "1se"))

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

res_experiment <- unlist(strsplit(experiment, "-"))
stopifnot(res_experiment[2] %in% c("biased", "unbiased"))
DGP <- res_experiment[1]
add.bias <- (res_experiment[2] == "biased")

#######
name_methods <- c("BU", "MINTols", "MINTshr", "BASE","BASE2", 
                  "L1", "L2",  "L1-PBU", "L2-PBU", "G-L1-PBU")
if(DGP != "large"){
  name_methods <- c(name_methods, "MINTsam", "OLS")
}
nb_methods <- length(name_methods)
# "NAIVE", "AVG", 
#######

config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
                        param_fit_fct = list(seasonal = FALSE, ic = "aic", max.p = 2, max.q = 2,  
                                     approximation = TRUE, stationary = FALSE), 
                        param_forecast_fct = list(use.initial.values = TRUE))
config_forecast_agg <- config_forecast_bot <- config_forecast

mc.cores.basef <- 10
mc.cores.methods <- 10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 3

H <- 3
do.save <- TRUE
refit_step <- 40
sameP_allhorizons <- TRUE
T_learn <- 600
T_train <- floor((2 * T_learn)/3)
T_valid <- T_learn - T_train
T_test <- 200
T_all <- T_train + T_valid + T_test
n_warm <- 300
n_simul <- n_warm + T_all

for(i in seq(nb_simulations)){
  
  print(paste(i, " - Start ALL -", base::date(), sep = ""))

  if(DGP == "small"){
    A <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1))
    obj_simul <- simulate_hts(n_simul)
    bts <- obj_simul$bts
  }else  if(DGP == "large"){
    res <- simulate_large_hts(n_simul)
    A <- res$A
    btsWithWarm <- res$bts
    bts <- tail(btsWithWarm, -n_warm)
  }else  if(DGP == "medium"){
    
    ngroups <- 20
    allbts <- vector("list", ngroups)
    for(igroup in seq(ngroups)){
      # save parameters in a list
      allbts[[igroup]] <- simulate_hts(n_simul)$bts
    }
    
    bts <- do.call(cbind, allbts)
    nbts <- ncol(bts)
    
    my_aggregation <- c(1, 4, 20)
    Alist <- lapply(my_aggregation, function(igroup){
      Stemp <- diag(igroup)
      res <- t(sapply(seq(nrow(Stemp)), function(i){
        rep(Stemp[i, ], each = nbts/igroup)
      }))
      res
    })
    A <- do.call(rbind, Alist)
    #A <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1))
  }else{
    stop("error in DGP !")
  }
  
  my_bights <- bights(bts, A)
 
  if(i == 1){
    info_file <- file.path(results.folder, paste("info_", experiment, "_",
                                                 idjob, "_", lambda_selection, ".Rdata", sep = ""))
    #save(file = info_file, list = c("name_methods", "A"))
    save(file = info_file, list = c("A"))
  }
  
  print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))
  print(paste("Tvalid = ", T_valid, sep = ""))
  print(paste("N = ", my_bights$nts * T_valid, " - p = ", my_bights$nbts * my_bights$nts, sep = ""))

  file_bf <- file.path(bf.folder, 
              paste("bf_", experiment, "_", idjob, ".", i, "_", refit_step, ".Rdata", sep = ""))  
  
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
  
  stop("done")
  
if(add.bias){
  #hat_all <- Y_valid_allh[1, , ]
  #PBU <- pbu(my_bights)
  #hat_bottom <- PBU %*% hat_all
  #mu_hat_bottom <- apply(hat_bottom, 1, mean) 
  #sd_hat_bottom <- apply(hat_bottom, 1, sd)
  
  #bias_mu_bottom <- rep(2, nbts); bias_sd_bottom <- rep(1, nbts)
  #bias_mu_bottom <- mu_hat_bottom/2; bias_sd_bottom <- sd_hat_bottom/2
  nbts <- my_bights$nbts
  if(DGP == "small"){
    bias_mu_bottom <- rep(1, nbts); bias_sd_bottom <- rep(1, nbts)
  }else if(DGP == "large"){
    # bias_mu_bottom <- rep(1, nbts); bias_sd_bottom <- rep(1, nbts) # error at 100
    # bias_mu_bottom <- rep(0.1, nbts); bias_sd_bottom <- rep(0.25, nbts) # error at 6
    # bias_mu_bottom <- rep(0.1, nbts); bias_sd_bottom <- rep(0.25, nbts) # error at 1.5
    # bias_mu_bottom <- rep(0.1, nbts); bias_sd_bottom <- rep(0.05, nbts) # error at 1
    bias_mu_bottom <- rep(0.05, nbts); bias_sd_bottom <- rep(0.05, nbts) # error at ?
  }
  
  bias_mu_agg <- A %*% bias_mu_bottom
  bias_sd_agg <- A %*% bias_sd_bottom
  bias_mu <- c(bias_mu_agg, bias_mu_bottom)
  bias_sd <- c(bias_sd_agg, bias_sd_bottom)

  nvalid <- ncol(Yhat_valid_allh[1, , ])
  MYBIAS_valid <- t(sapply(seq(my_bights$nts), function(j){ 
    rnorm(nvalid, mean = bias_mu[j], sd = bias_sd[j])
  }))
  Yhat_valid_allh[1, , ]  <- Yhat_valid_allh[1, , ]  + MYBIAS_valid
  
  ntest <- ncol(Yhat_test_allh[1, , ])
  MYBIAS_test <- t(sapply(seq(my_bights$nts), function(j){ 
    rnorm(ntest, mean = bias_mu[j], sd = bias_sd[j])
  }))
  Yhat_test_allh[1, , ] <- Yhat_test_allh[1, , ] + MYBIAS_test  
}

config_basic <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
config <- list(glmnet = config_basic, cvglmnet = c(config_basic, list(nfolds = 3)))
config_gglasso <- list(eps = 10^-6, nlambda = 50, intercept = FALSE, maxit = 1e+7)

#sgl_config    <- list(nfold = 3, standardize = FALSE, alpha = .8)
#print(glmnet_config)
#print(sgl_config)

# naive predictions
id <- sapply(list_subsets_test, function(vec){ vec[2] })
predictions_naive <- my_bights$yts[id, ]

predictions_avg <- t(sapply(list_subsets_test, function(vec){ 
         apply(my_bights$yts[vec[1]:vec[2], ], 2, mean)
}))

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
      
      list(obj_learn = obj_learn, predictions = predictions)
    }, mc.cores = mc.cores.methods)
    names(results) <- name_methods
    
    results[["BASE"]]   <- list(obj_learn = NULL, predictions = t(Yhat_test_allh[h, , ]))
    results[["BASE2"]]   <- list(obj_learn = NULL, predictions = t(save_Yhat_test_allh[h, , ]))

    results_allh[[h]] <- results
  
  } # HORIZON


  myfile <- file.path(results.folder, paste("results_", experiment, "_",
                                            idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
  save(file = myfile, list = c("results_allh", "Y_test_allh", "config_forecast_agg", "config_forecast_bot"))
  
} # simulation


