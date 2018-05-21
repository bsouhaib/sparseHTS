rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  # hierarchy_names <- c("UKF23", "UKF16", "UKF21", "UKF13", "UKF15")
  hierarchy_name <- "UKF13"
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  hierarchy_name <- args[[1]]
}
source("config_paths.R")
source("packages.R")
source("bights.R")
source("methods_meters.R") ##############
source("simulate.R")
source("simulate_large.R")
source("utils.R")
source("hts.R")
source("code.R")
library(methods)
library(doMC)
library(gglasso)
source("meters_utils.R")


do.diff <- TRUE
if(do.diff){
  print("DIFF IS APPLIED !!!!!")
}

set.seed(1986)

source("../../PROJ/code/config_general.R")
source(file.path("/home/rstudio/PROJ/code", "config_splitting.R"))

#load(file.path(file.path("/home/rstudio/PROJ", "work"), "myinfo.Rdata"))
sparsehts.work.folder <- file.path("/home/rstudio/sparseHTS", "work")
load(file.path(sparsehts.work.folder, paste("myinfo_", hierarchy_name,".Rdata", sep = "")))


add.bias <- FALSE
if(add.bias){
  print("BIAS IS ADDED !!!")
}

lambda_selection <- "1se"
#lambda_selection <- "min"
H <- 48
mc.cores.methods <- 1 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <- 3
sameP_allhorizons <- FALSE

file_bf <- file.path(bf.folder,  paste("bf_", hierarchy_name, ".Rdata", sep = ""))  
load(file_bf)

#######
name_methods <- c("BU", "MINTshr", "MINTols",  "L1-PBU", "L1")
nb_methods <- length(name_methods)
#######

#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = 1, thresh = 10^-6, nlambda = 50, weights = weights_glmnet)
config_main <- list(intercept = FALSE, standardize = FALSE, alpha = 1, thresh = 10^-6, nlambda = 100)
config_cvglmnet <- c(config_main, list(nfolds = 3))




A <- Sagg
bts <- matrix(NA, nrow = length(learn$id) + length(test$id), ncol = n_bottom)
my_bights <- bights(bts, A)

print(my_bights$nbts)
print(my_bights$nts)

info_file <- file.path(results.folder, paste("info_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
#save(file = info_file, list = c("name_methods", "A"))
save(file = info_file, list = c("A"))
  
#file_bf <- file.path(bf.folder,  paste("bf_", experiment, ".Rdata", sep = ""))  
#load(file_bf)
  
  ### valid
  Yhat_valid_allh <- data_valid$Yhat
  Y_valid_allh     <- data_valid$Y
  
  ### test
  Yhat_test_allh  <- data_test$Yhat
  Y_test_allh    <- data_test$Y
  
  Eresid <- data_test$Eresiduals
  
  save_Yhat_test_allh <- Yhat_test_allh
  
  if(do.diff){
    makediff <- function(myarray){
      res <- sapply(seq(H), function(h){
        apply(myarray[h, , ], 1, diff)
      }, simplify = "array")
      res <- aperm(res, c(3, 2, 1))
    }
    
    Yhat_valid_allh  <- makediff(Yhat_valid_allh)
    Y_valid_allh     <- makediff(Y_valid_allh)
    Yhat_test_allh   <- makediff(Yhat_test_allh)
    Y_test_allh      <- makediff(Y_test_allh)
    save_Yhat_test_allh <- makediff(save_Yhat_test_allh)
    
    Eresid <- apply(Eresid, 2, diff, lag = 48)
    pdays <- tail(calendar$periodOfDay[tail(learn$id, -n_past_obs_kd)], -48)
    
  }
  
  if(add.bias){
    nbts <- my_bights$nbts
    nts <- my_bights$nts
    A <- my_bights$A
    bias_mu_bottom <- rep(0.3, nbts); bias_sd_bottom <- rep(0.05, nbts) # error at ?
    bias_mu_agg <- as.numeric(A %*% bias_mu_bottom)
    bias_sd_agg <- as.numeric(A %*% bias_sd_bottom)
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
  
  if(do.diff){
    alldim <- dim(Yhat_test_allh)
    alldim[3] <- alldim[3] - 1
    Ytilde_test_allh <- array(NA, c(alldim, nb_methods) )
  }else{
    Ytilde_test_allh <- array(NA, c(dim(Yhat_test_allh), nb_methods) )
  }
  results_allh <- vector("list", H)
  for(h in seq(H)){
    print(paste("h = ", h, sep = ""))
    
    if(do.diff){
      Eresiduals <- Eresid[which(pdays == h), ]
    }
    
    Yhat_valid_h <- t(Yhat_valid_allh[h, , ])
    Y_valid_h    <- t(Y_valid_allh[h, , ])
    Yhat_test_h  <- t(Yhat_test_allh[h, , ])
    
    
    #stop("done")
    #ids <- which(calendar$periodOfDay[tail(learn$id, - n_past_obs_kd)] == pday_touse(h))
    #Eresiduals <- data_test$Eresiduals[ids, ]
    #stop("done")
    
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
      }else if(current_method == "NEW"){
        #obj_method <- list(algo = "NEW", lambda1 = 0.2, Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
        #obj_method <- list(algo = "NEW", lambda1 = 0, Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
        #obj_method <- list(algo = "NEW", set_lambda1 = c(0, 0.0001, 0.001, 0.01, 0.02, 0.05, 0.09, 0.1), Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
        #obj_method <- list(algo = "NEW", set_lambda1 = c(0.1, 0.2, 0.3, 0.4, 1, 2, 3), Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
        obj_method <- list(algo = "NEW", set_lambda1 = c(1, 2, 3, 5, 6, 10), Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
        
      }else if(current_method == "L1"){
        obj_method <- list(algo = "LASSO", Ptowards = NULL, config = config, selection = lambda_selection)
      }else if(current_method == "L1-PBU"){
        obj_method <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
      }else if(current_method == "L1-POLS"){
        pols <- solve(t(my_bights$S) %*% my_bights$S) %*% t(my_bights$S)
        obj_method <- list(algo = "LASSO", Ptowards = pols, config = config, selection = lambda_selection)
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
          #obj_learn <- results[[current_method]]$obj_learn
          obj_learn <- results_h1[[current_method]]
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
    
    if(sameP_allhorizons && h == 1){
      results_h1 <-  lapply(results , "[[", "obj_learn")
    }
    
    results[["BASE"]]   <- list(obj_learn = NULL, predictions = t(Yhat_test_allh[h, , ]))
    results[["BASE2"]]   <- list(obj_learn = NULL, predictions = t(save_Yhat_test_allh[h, , ]))
    
    # results_allh[[h]] <-  results ONLY NEED PREDICTIONS + SAVE SPACE
    results_allh[[h]] <- lapply(results , "[[", "predictions")
    
    if(h%%5 == 0){
      myfile <- file.path(results.folder, paste("results_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
      save(file = myfile, list = c("results_allh", "Y_test_allh"))
    }
    
  } # HORIZON
  
  myfile <- file.path(results.folder, paste("results_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
  save(file = myfile, list = c("results_allh", "Y_test_allh"))


