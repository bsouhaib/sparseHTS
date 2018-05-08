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

source("../../PROJ/code/config_general.R")

#######
#name_methods <- c("L1-PBU", "L1", "BU", "MINTols", "MINTshr", "BASE","BASE2", 
#                  "L2", "L2-PBU")
#name_methods <- c("L1-PBU", "BU", "MINTols", "MINTshr", "BASE","BASE2", "L2-PBU")
#, "G-L1-PBU"

#name_methods <- c("L1-PBU", "BU", "MINTshr", "MINTols", "BASE2")
#name_methods <- c("L1-PBU", "BU", "MINTshr", "MINTols", "BASE2")
name_methods <- c("L1-PBU", "BU", "MINTols", "MINTshr", "BASE","BASE2", "L2-PBU", "L1", "L2")
nb_methods <- length(name_methods)
#######

#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
config_main <- list(intercept = FALSE, standardize = FALSE, alpha = 1, thresh = 10^-6, nlambda = 50)
config_cvglmnet <- c(config_main, list(nfolds = 3))

lambda_selection <- "1se"
#lambda_selection <- "min"
experiment <- "meters"
H <- 5

mc.cores.methods <- 5 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 3

sameP_allhorizons <- FALSE


#load(file.path(file.path("/home/rstudio/PROJ", "work"), "myinfo.Rdata"))
sparsehts.work.folder <- file.path("/home/rstudio/sparseHTS", "work")
load(file.path(sparsehts.work.folder, "myinfo.Rdata"))

source(file.path("/home/rstudio/PROJ/code", "config_splitting.R"))


A <- Sagg
bts <- matrix(NA, nrow = length(learn$id) + length(test$id), ncol = n_bottom)
my_bights <- bights(bts, A)

info_file <- file.path(results.folder, paste("info_", experiment, "_", lambda_selection, ".Rdata", sep = ""))
#save(file = info_file, list = c("name_methods", "A"))
save(file = info_file, list = c("A"))
  
  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, ".Rdata", sep = ""))  
  load(file_bf)
  
  ### valid
  Yhat_valid_allh <- data_valid$Yhat
  Y_valid_allh     <- data_valid$Y
  
  ### test
  Yhat_test_allh  <- data_test$Yhat
  Y_test_allh    <- data_test$Y
  

  #Eresiduals <- data_test$Eresiduals
  ids <- which(calendar$periodOfDay[tail(learn$id, - n_past_obs_kd)] == H)
  Eresiduals <- data_test$Eresiduals[ids, ]
  
  save_Yhat_test_allh <- Yhat_test_allh
  
  # 
  #stop("HERE")
  #data_valid$IN_y
  #data_valid$IN_yhat
  
  config_basic <- config_main
  config <- list(glmnet = config_basic, cvglmnet = config_cvglmnet)
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
    
    objreg_valid <- list(y = makey(Y_valid_h), X = makeX(Yhat_valid_h, my_bights), Y = Y_valid_h, Yhat = Yhat_valid_h)
    objreg_test <- list(X = makeX(Yhat_test_h, my_bights), Yhat = Yhat_test_h)
    
    #browser()
    #foldid <- rep(sample(seq(config$cvglmnet$nfolds), nValid, replace = T), my_bights$nts)
    #glmnet_config_local <- c(glmnet_config, list(foldid = foldid))
    #glmnet_config_local$nfolds <- NA
    nValid <- nrow(Yhat_valid_h)
    foldid <- make_foldid(nValid, my_bights$nts, config$cvglmnet$nfolds)
    
    config$cvglmnet$foldid <- foldid
    #config$cvglmnet$nfolds <- NA
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
  
  
  myfile <- file.path(results.folder, paste("results_", experiment, "_", lambda_selection, ".Rdata", sep = ""))
  save(file = myfile, list = c("results_allh", "Y_test_allh"))

