rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  
  experiment <- "large-biased"
  idjob <- 1
  nb_simulations <- 100
  fmethod_agg <- "ARIMA"
  fmethod_bot <- "ARIMA"
  lambda_selection <- "1se"
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  idjob <- as.numeric(args[[2]])
  nb_simulations <- as.numeric(args[[3]])
  fmethod_agg <- args[[4]]
  fmethod_bot <- args[[5]]
  lambda_selection <- args[[6]]
}

#set.seed(6120)
set.seed(idjob)

#print(experiment)
#print(idjob)
#print(nb_simulations)
#print(fmethod_agg)
#print(fmethod_bot)
#print(lambda_selection)

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

name_methods <- c("BU", "MINTOLS", "MINTSAM", "MINTSHR", "NAIVE", "AVG", 
                  "OLS", "OLSS", "gOLS",  "L1", "L2",  "L1-PBU", "L2-PBU", "SGLR", "SGLC", "GGLAS", "GGLAS_PBU", 
                  "BASE","BASE2", "OPT")
nb_methods <- length(name_methods)

res_experiment <- unlist(strsplit(experiment, "-"))
stopifnot(res_experiment[2] %in% c("biased", "unbiased"))
DGP <- res_experiment[1]
add.bias <- (res_experiment[2] == "biased")

mc.cores <- 1
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
 

#M <- 8 * 15
#results <- mclapply(seq(M), function(iteration){
#results <- sapply(seq(M), function(iteration){

results <- vector("list", nb_simulations)
nbzeroes <- vector("list", nb_simulations)
for(i in seq_along(results)){
  
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
  H <- 1
  print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))
  print(paste("Tvalid = ", T_valid, sep = ""))
  print(paste("N = ", my_bights$nts * T_valid, " - p = ", my_bights$nbts * my_bights$nts, sep = ""))

  file_bf <- file.path(rdata.folder, 
              paste("bf_", experiment, "_", fmethod_agg, "_", 
              fmethod_bot, "_", idjob, ".", i, "_", refit_step, ".Rdata", sep = ""))  
  
  if(do.save && file.exists(file_bf)){
    load(file_bf)
  }else{
    ### valid
    list_subsets_valid <- lapply(seq(T_train, T_learn - H), function(i){c(i - T_train + 1, i)})
    data_valid <- makeMatrices(my_bights, list_subsets_valid, H = H, 
                               fmethod_agg = fmethod_agg, fmethod_bot = fmethod_bot, refit_step = refit_step, mc.cores = mc.cores)
    Yhat_valid_allh <- data_valid$Yhat
    Y_valid_allh     <- data_valid$Y
    
    ### test
    list_subsets_test <- lapply(seq(T_learn, T_all - H), function(i){c(i - T_learn + 1, i)}) # I CHANGED it from T_TRAIN TO T_LEARN !!
    data_test <- makeMatrices(my_bights, list_subsets_test, H = H, 
                              fmethod_agg = fmethod_agg, fmethod_bot = fmethod_bot, refit_step = refit_step, mc.cores = mc.cores)
    Yhat_test_allh  <- data_test$Yhat
    Y_test_allh    <- data_test$Y
    
    if(do.save){
      save(file = file_bf, list = c("data_valid", "Yhat_valid_allh", 
                                  "Y_valid_allh", "data_test", "Yhat_test_allh", "Y_test_allh"))
    }
  }

save_Yhat_test_allh <- Yhat_test_allh

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

objlearn_glmnet <- NULL
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

  ##
  if(!sameP_allhorizons || (sameP_allhorizons && h == 1) ){
    
    #objmethod <- list(algo = "GGLASSO", Ptowards = NULL, config = config_gglasso, selection = lambda_selection)
    #objlearn_GGLASSO <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #objmethod <- list(algo = "GGLASSO", Ptowards = pbu(my_bights), config = config_gglasso, selection = lambda_selection)
    #objlearn_GGLASSO_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #print(paste("Start LASSO -", base::date(), sep = ""))
    # LASSO
    objmethod <- list(algo = "LASSO", Ptowards = NULL, config = config, selection = lambda_selection)
    objlearn_LASSO <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #print(paste("Start LASSO-PBU -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
    objlearn_LASSO_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #print(paste("Start LS -", base::date(), sep = ""))
    # OLS
    #objmethod <- list(algo = "gOLS", Ptowards = NULL, config = config, selection = lambda_selection)
    #objlearn_gOLS <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    config_ridge <- config
    config_ridge$glmnet$alpha <- 0
    config_ridge$cvglmnet$alpha <- 0 # !!!!!!!
    
    #print(paste("Start RIDGE -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = NULL, config = config_ridge, selection = lambda_selection)
    objlearn_RIDGE <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #print(paste("Start RIDGE-PBU -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config_ridge, selection = lambda_selection)
    objlearn_RIDGE_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
  }
  
  predictions_LASSO <- new_predtest(objreg_test, my_bights, objlearn_LASSO)
  predictions_LASSO_towardspbu <- new_predtest(objreg_test, my_bights, objlearn_LASSO_towardspbu)
  #predictions_gOLS <- new_predtest(objreg_test, my_bights, objlearn_gOLS)
  predictions_RIDGE <- new_predtest(objreg_test, my_bights, objlearn_RIDGE)
  predictions_RIDGE_towardspbu <- new_predtest(objreg_test, my_bights, objlearn_RIDGE_towardspbu)
  
  #print(paste("END PREDICTIONS -", base::date(), sep = ""))
  
  #predictions_GGLASSO <- 0
  #predictions_GGLASSO_towardspbu <- 0
  #if(DGP == "small"){
  #  predictions_GGLASSO <- new_predtest(objreg_test, my_bights, objlearn_GGLASSO)
  #  predictions_GGLASSO_towardspbu <- new_predtest(objreg_test, my_bights, objlearn_GGLASSO_towardspbu)
  #}
  
  #print(date())
  
  ##
  if(FALSE){
    objlearn_sglrow <- learnreg(objreg_valid, my_bights, "SGLrow", 
                         config = sgl_config, 
                         selection = lambda_selection, towards_pbu = towards_pbu)
    predictions_sglrow <- predtest(objreg_test, my_bights, objlearn_sglrow)
    
    ##
    objlearn_sglcol <- learnreg(objreg_valid, my_bights, "SGLcol", 
                                config = sgl_config, 
                                selection = lambda_selection, towards_pbu = towards_pbu)
    predictions_sglcol <- predtest(objreg_test, my_bights, objlearn_sglcol)
  }
  
  # BU
  obj_bu <- bu(my_bights)
  predictions_bu <- new_predtest(objreg_test, my_bights, obj_bu)
  
  #print(date())
  
  if(DGP == "large"){
    predictions_mintsample <- NA
  }else{
    # MINTsample
    obj_mintsample <- mint(my_bights, method = "sample", residuals = residuals, h = h) # J U and W
    predictions_mintsample <- new_predtest(objreg_test, my_bights, obj_mintsample)
    #print(date())
  }
  
  # MINTshrink
  residuals <- Y_valid_allh - Yhat_valid_allh
  obj_mintshrink <- mint(my_bights, method = "shrink", residuals = residuals, h = h) # J U and W
  predictions_mintshrink <- new_predtest(objreg_test, my_bights, obj_mintshrink)
  
  # MINTols
  obj_mintols <- mint(my_bights, method = "ols", h = h) # J U and W
  predictions_mintols <- new_predtest(objreg_test, my_bights, obj_mintols)
  
  if(DGP == "large"){
    predictions_ols <- NA
  }else{
    # OLS
    obj_ols <- ols(objreg_valid, my_bights) 
    predictions_ols <- new_predtest(objreg_test, my_bights, obj_ols)
  }
  
  # OLSS
  #objmethod <- list(algo = "OLSS", Ptowards = NULL, config = NULL, selection = NULL)
  #objlearn_olss <- new_learnreg(objreg_valid, my_bights, objmethod)
  #objlearn_olss_bis <- ols(objreg_valid, my_bights) 
  #objlearn_olss <- c(objlearn_olss, objlearn_olss_bis)
  #predictions_olss<- new_predtest(objreg_test, my_bights, objlearn_olss)

  Ytilde_test_allh[h, , , 1] <- t(predictions_bu) 
  
  Ytilde_test_allh[h, , , 2] <- t(predictions_mintols)
  Ytilde_test_allh[h, , , 3] <- t(predictions_mintsample)
  Ytilde_test_allh[h, , , 4] <- t(predictions_mintshrink)
  Ytilde_test_allh[h, , , 5] <- t(predictions_naive)
  Ytilde_test_allh[h, , , 6] <- t(predictions_avg)
  
  Ytilde_test_allh[h, , , 7] <- t(predictions_ols)
  #Ytilde_test_allh[h, , , 8] <- t(predictions_olss)
  #Ytilde_test_allh[h, , , 9]  <- t(predictions_gOLS)
  
  Ytilde_test_allh[h, , , 10] <- t(predictions_LASSO)
  Ytilde_test_allh[h, , , 11] <- t(predictions_RIDGE)
  
  Ytilde_test_allh[h, , , 12]  <- t(predictions_LASSO_towardspbu)
  Ytilde_test_allh[h, , , 13]  <- t(predictions_RIDGE_towardspbu)
  
  #Ytilde_test_allh[h, , , 14]   <- t(predictions_sglrow)
  #Ytilde_test_allh[h, , , 15]   <- t(predictions_sglcol)
  #Ytilde_test_allh[h, , , 16]   <- t(predictions_GGLASSO)
  #Ytilde_test_allh[h, , , 17]   <- t(predictions_GGLASSO_towardspbu)
}
Ytilde_test_allh[, , , 18] <- Yhat_test_allh
Ytilde_test_allh[, , , 19] <- save_Yhat_test_allh
#Ytilde_test_allh[, , , 20] <- t(predictions_optimal)
dimnames(Ytilde_test_allh)[[4]] <- name_methods

  #results[[i]] <- list(Ytilde_test_allh = Ytilde_test_allh, Y_test_allh = Y_test_allh)
  myfile <- file.path(rdata.folder, paste("results_", experiment, "_", 
                                          fmethod_agg, "_", fmethod_bot, "_", idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
  save(file = myfile, list = c("Ytilde_test_allh", "Y_test_allh"))
}

#print(as.numeric(nbzeroes))
#}, simplify = "array")
#}, mc.cores = 8)
#stop("done")
#myfile <- file.path(getwd(), "../work", paste("results_", idjob, "_", towards_pbu, "_", lambda_selection, ".Rdata", sep = ""))
#save(file = myfile, list = c("results"))


