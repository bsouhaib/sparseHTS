rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0){
  
  experiment <- "small"
  idjob <- 220
  nb_simulations <- 500
  fmethod_agg <- "AR1"
  fmethod_bot <- "AR1"
  lambda_selection <- "1se"
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  idjob <- as.integer(args[[2]])
  nb_simulations <- as.integer(args[[3]])
  fmethod_agg <- args[[4]]
  fmethod_bot <- args[[5]]
  lambda_selection <- args[[6]]
}

#set.seed(6120)
set.seed(idjob)

#print(idjob)
#print(lambda_selection)
#print(towards_pbu)
#print(nb_simulations)

stopifnot(lambda_selection %in% c("min", "1se"))


source("packages.R")
source("bights.R")
source("methods.R")
source("simulate.R")
source("simulate_large.R")
source("utils.R")
source("hts.R")
source("code.R")

nb_methods <- 13

sameP_allhorizons <- TRUE
refit_step <- 40

T_train <- 300  
T_valid <- 200
#T_valid <- 300
T_test <- 200
T_learn <- T_train + T_valid
T_all <- T_train + T_valid + T_test
n_warm <- 300
n_simul <- n_warm + T_all
 

#M <- 8 * 15
#results <- mclapply(seq(M), function(iteration){

#results <- sapply(seq(M), function(iteration){

results <- vector("list", nb_simulations)
nbzeroes <- vector("list", nb_simulations)
for(i in seq_along(results)){
  
  print(i)

  if(experiment == "small"){
    A <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1))
    bts <- simulte_hts(n_simul)
  }else  if(experiment == "large"){
    res <- simulte_large_hts(n_simul)
    A <- res$A
    btsWithWarm <- res$bts
    bts <- tail(btsWithWarm, -n_warm)
  }
  
  my_bights <- bights(bts, A)
  H <- 2

# OLD: data_valid <- make.data(my_bights, list_subsets_valid, H = H)
# OLD: data_test <- make.data(my_bights, list_subsets_test, H = H)

print(paste("Start forecasting in validation -", base::date(), sep = ""))
### valid
list_subsets_valid <- lapply(seq(T_train, T_learn - H), function(i){c(i - T_train + 1, i)})
data_valid <- makeMatrices(my_bights, list_subsets_valid, H = H, fmethod_agg = fmethod_agg, fmethod_bot = fmethod_bot, refit_step = refit_step)
Yhat_valid_allh <- data_valid$Yhat
Y_valid_allh     <- data_valid$Y

print(paste("Start forecasting in testing -", base::date(), sep = ""))
### test
list_subsets_test <- lapply(seq(T_learn, T_all - H), function(i){c(i - T_learn + 1, i)}) # I CHANGED it from T_TRAIN TO T_LEARN !!
data_test <- makeMatrices(my_bights, list_subsets_test, H = H, fmethod_agg = fmethod_agg, fmethod_bot = fmethod_bot, refit_step = refit_step)
Yhat_test_allh  <- data_test$Yhat
Y_test_allh    <- data_test$Y

# save this in rdata? + tag ?

#glmnet_config <- list(intercept = TRUE, standardize = TRUE, alpha = .98, nfolds = 3, thresh = 10^-5)
#config <- list(intercept = FALSE, standardize = FALSE)
#glmnet_config <- c(config, list(alpha = .98, nfolds = 3))
#glmnet_configOLS <- config

config_basic <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6)
config <- list(glmnet = config_basic, cvglmnet = c(config_basic, list(nfolds = 3)))

#sgl_config    <- list(nfold = 3, standardize = FALSE, alpha = .8)
#print(glmnet_config)
#print(sgl_config)


# naive predictions
id <- sapply(list_subsets_test, function(vec){vec[2]})
predictions_naive <- my_bights$yts[id, ]

# add mean predictions


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
  
  #print(date())
  
  nValid <- nrow(Yhat_valid_h)
  foldid <- rep(sample(seq(config$cvglmnet$nfolds), nValid, replace = T), my_bights$nts)
  #glmnet_config_local <- c(glmnet_config, list(foldid = foldid))
  #glmnet_config_local$nfolds <- NA
  config$cvglmnet$foldid <- foldid
  config$cvglmnet$nfolds <- NA
  
 # print(date())
  #stop("done")
  ##
  if(!sameP_allhorizons || (sameP_allhorizons && h == 1) ){
    
    #objmethod <- list(algo = "GGLASSO", Ptowards = NULL, config = config, selection = lambda_selection)
    #objlearn_GGLASSO <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    #objmethod <- list(algo = "GGLASSO", Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
    #objlearn_GGLASSO_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    print(paste("Start LASSO -", base::date(), sep = ""))
    # LASSO
    objmethod <- list(algo = "LASSO", Ptowards = NULL, config = config, selection = lambda_selection)
    objlearn_LASSO <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    print(paste("Start LASSO-PBU -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config, selection = lambda_selection)
    objlearn_LASSO_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    print(paste("Start LS -", base::date(), sep = ""))
    # OLS
    objmethod <- list(algo = "OLS", Ptowards = NULL, config = config, selection = lambda_selection)
    objlearn_OLS <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    config_ridge <- config
    config_ridge$glmnet$alpha <- 0
    config_ridge$cvglmnet$alpha <- 0 # !!!!!!!
    
    print(paste("Start RIDGE -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = NULL, config = config_ridge, selection = lambda_selection)
    objlearn_RIDGE <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    print(paste("Start RIDGE-PBU -", base::date(), sep = ""))
    objmethod <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config_ridge, selection = lambda_selection)
    objlearn_RIDGE_towardspbu <- new_learnreg(objreg_valid, my_bights, objmethod)
    
    
  
  }
  
  predictions_LASSO <- new_predtest(objreg_test, my_bights, objlearn_LASSO)
  predictions_LASSO_towardspbu <- new_predtest(objreg_test, my_bights, objlearn_LASSO_towardspbu)
  predictions_OLS <- new_predtest(objreg_test, my_bights, objlearn_OLS)
  predictions_RIDGE <- new_predtest(objreg_test, my_bights, objlearn_RIDGE)
  predictions_RIDGE_towardspbu <- new_predtest(objreg_test, my_bights, objlearn_RIDGE_towardspbu)
  
  print(paste("END PREDICTIONS -", base::date(), sep = ""))
  
  #predictions_GGLASSO <- 0
  #predictions_GGLASSO_towardspbu <- 0
  #if(experiment == "small"){
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
  
  # MINT
  residuals <- Y_valid_allh - Yhat_valid_allh
  obj_mint <- mint(my_bights, method = "shrink", residuals = residuals, h = h) # J U and W
  predictions_mint <- new_predtest(objreg_test, my_bights, obj_mint)
  
  #print(date())
  
  obj_mintols <- mint(my_bights, method = "ols", h = h) # J U and W
  predictions_mintols <- new_predtest(objreg_test, my_bights, obj_mintols)
  
  #print(date())
  
  #if(i == 37)
  #  stop("STOP HERE")
  
  # LS
  #obj_ls <- mls(objreg_valid, my_bights) 
  #predictions_ls <- predtest(objreg_test, my_bights, obj_ls)

  Ytilde_test_allh[h, , , 1] <- t(predictions_LASSO)
  Ytilde_test_allh[h, , , 2] <- t(predictions_LASSO_towardspbu)
  Ytilde_test_allh[h, , , 3] <- t(predictions_RIDGE)
  
  #Ytilde_test_allh[h, , , 2] <- t(predictions_sglrow)
  #Ytilde_test_allh[h, , , 3] <- t(predictions_sglcol)
  Ytilde_test_allh[h, , , 4] <- t(predictions_bu)
  Ytilde_test_allh[h, , , 5] <- t(predictions_mint)
  #Ytilde_test_allh[h, , , 6] <- t(predictions_ls)
  Ytilde_test_allh[h, , , 7] <- t(predictions_naive)
  Ytilde_test_allh[h, , , 8] <- t(predictions_OLS)
  Ytilde_test_allh[h, , , 9] <- t(predictions_mintols)
  
  #Ytilde_test_allh[h, , , 10] <- t(predictions_GGLASSO)
  #Ytilde_test_allh[h, , , 11] <- t(predictions_GGLASSO_towardspbu)
  Ytilde_test_allh[h, , , 12] <- t(predictions_RIDGE_towardspbu)
  #print(date())
  
  print(paste("FINISH ALL -", base::date(), sep = ""))
}
Ytilde_test_allh[, , , 13] <- Yhat_test_allh



#results[[i]] <- list(Ytilde_test_allh = Ytilde_test_allh, Y_test_allh = Y_test_allh)
myfile <- file.path(getwd(), "../work", paste("results_", experiment, "_", fmethod_agg, "_", fmethod_bot, "_", idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
save(file = myfile, list = c("Ytilde_test_allh", "Y_test_allh"))
}

#print(as.numeric(nbzeroes))

#}, simplify = "array")
#}, mc.cores = 8)
#stop("done")

#myfile <- file.path(getwd(), "../work", paste("results_", idjob, "_", towards_pbu, "_", lambda_selection, ".Rdata", sep = ""))
#save(file = myfile, list = c("results"))


