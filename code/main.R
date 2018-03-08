rm(list = ls())
args = (commandArgs(TRUE))
if(length(args) == 0){
  
  idjob <- 1
  towards_pbu <- TRUE
  lambda_selection <- "min"
  nb_simulations <- 2
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  idjob <- args[[1]]
  lambda_selection <- args[[2]]
  towards_pbu <- args[[3]]
  nb_simulations <- args[[4]]
}

stopifnot(lambda_selection %in% c("min", "1se"))


source("packages.R")
source("bights.R")
source("methods.R")
source("simulate.R")
source("simulate_large.R")
source("utils.R")
source("hts.R")
source("code.R")

#set.seed(6120)
set.seed(idjob)

experiment <- "2"

fmethod <- "ARIMA"
fmethod <- "ETS"

refit_step <- 40

T_train <- 300  
T_valid <- 200
T_test <- 200
T_learn <- T_train + T_valid
T_all <- T_train + T_valid + T_test
n_warm <- 300
n_simul <- n_warm + T_all
 

#M <- 8 * 15
#results <- mclapply(seq(M), function(iteration){

#results <- sapply(seq(M), function(iteration){

results <- vector("list", nb_simulations)
for(i in seq_along(results)){
  
  print(i)

  if(experiment == "1"){
    A <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1))
    bts <- simulte_hts(n_simul)
  }else  if(experiment == "2"){
    res <- simulte_large_hts(n_simul)
    A <- res$A
    btsWithWarm <- res$bts
    bts <- tail(btsWithWarm, -n_warm)
  }
  
  my_bights <- bights(bts, A)
  H <- 10

# OLD: data_valid <- make.data(my_bights, list_subsets_valid, H = H)
# OLD: data_test <- make.data(my_bights, list_subsets_test, H = H)

### valid
list_subsets_valid <- lapply(seq(T_train, T_learn - H), function(i){c(i - T_train + 1, i)})
data_valid <- makeMatrices(my_bights, list_subsets_valid, H = H, fmethod = fmethod, refit_step = refit_step)
Yhat_valid_allh <- data_valid$Yhat
Y_valid_allh     <- data_valid$Y

### test
list_subsets_test <- lapply(seq(T_learn, T_all - H), function(i){c(i - T_learn + 1, i)}) # I CHANGED it from T_TRAIN TO T_LEARN !!
data_test <- makeMatrices(my_bights, list_subsets_test, H = H, fmethod = fmethod, refit_step = refit_step)
Yhat_test_allh  <- data_test$Yhat
Y_test_allh    <- data_test$Y

# save this in rdata? + tag ?

glmnet_config <- list(intercept = FALSE, standardize = FALSE, alpha = 1, nfolds = 3, thresh = 10^-5)
sgl_config    <- list(nfold = 3, standardize = FALSE, alpha = .8)



nb_methods <- 7

Ytilde_test_allh <- array(NA, c(dim(Yhat_test_allh), nb_methods) )

for(h in seq(H)){
  print(paste("h = ", h, sep = ""))

  Yhat_valid_h <- t(Yhat_valid_allh[h, , ])
  Y_valid_h    <- t(Y_valid_allh[h, , ])
  Yhat_test_h  <- t(Yhat_test_allh[h, , ])
  
  objreg_valid <- list(y = makey(Y_valid_h), X = makeX(Yhat_valid_h, my_bights), Y = Y_valid_h, Yhat = Yhat_valid_h)
  objreg_test <- list(X = makeX(Yhat_test_h, my_bights), Yhat = Yhat_test_h)
  
  
  nValid <- nrow(Yhat_valid_h)
  foldid <- rep(sample(seq(glmnet_config$nfolds), nValid, replace = T), my_bights$nts)
  glmnet_config_local <- c(glmnet_config, list(foldid = foldid))
  glmnet_config_local$nfolds <- NA
  
  ##
  objlearn_glmnet <- learnreg(objreg_valid, my_bights, "glmnet", 
                       config = glmnet_config_local, 
                       selection = lambda_selection, towards_pbu = towards_pbu)
  predictions_glmnet <- predtest(objreg_test, my_bights, objlearn_glmnet)
  
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
  predictions_bu <- predtest(objreg_test, my_bights, obj_bu)
  
  # MINT
  residuals <- Y_valid_allh - Yhat_valid_allh
  obj_mint <- mint(my_bights, method = "shrink", residuals = residuals, h = h) # J U and W
  predictions_mint <- predtest(objreg_test, my_bights, obj_mint)
  
  # LS
  obj_ls <- mls(objreg_valid, my_bights) 
  predictions_ls <- predtest(objreg_test, my_bights, obj_ls)

  Ytilde_test_allh[h, , , 1] <- t(predictions_glmnet)
  #Ytilde_test_allh[h, , , 2] <- t(predictions_sglrow)
  #Ytilde_test_allh[h, , , 3] <- t(predictions_sglcol)
  Ytilde_test_allh[h, , , 4] <- t(predictions_bu)
  Ytilde_test_allh[h, , , 5] <- t(predictions_mint)
  Ytilde_test_allh[h, , , 6] <- t(predictions_ls)
  
}
Ytilde_test_allh[, , , 7] <- Yhat_test_allh

# compare with naive forecasts
# v <- matrix(rep(apply(my_bights$yts[seq(T_learn), ], 2, mean), length(list_subsets_valid)), ncol = 7, byrow = T)

  results[[i]] <- list(Ytilde_test_allh = Ytilde_test_allh, Y_test_allh = Y_test_allh)
}

#}, simplify = "array")
#}, mc.cores = 8)
#stop("done")

myfile <- file.path(getwd(), "../work", paste("results_", idjob, "_", towards_pbu, "_", lambda_selection, ".Rdata", sep = ""))
save(file = myfile, list = c("results"))

