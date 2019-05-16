rm(list = ls())
assign("last.warning", NULL, envir = baseenv())
args = (commandArgs(TRUE))
if(length(args) == 0 || length(args) == 2){
  
  if(TRUE){
    if(length(args) == 2){
      experiment <- args[[1]]
      algobf <- args[[2]]
    }else{
      #experiment <- "tourism-Vis" #"tourism-Bus" #"tourism-Oth"  #"tourism-Hol"  
      experiment <- "wikipedia-1"
      #experiment <- "road_traffic-1"
      #experiment <- "elec-1"
      algobf <- "arima"
    }
    add.bias <- FALSE
    
    #do.log <- ifelse(experiment == "tourism", FALSE, TRUE)
    do.log <- TRUE
    
    do.cleaning <- TRUE
    do.deseasonalization <- TRUE
    do.scaling <- FALSE
  }else{
    experiment <- "small-1"
    add.bias <- FALSE
    algobf <- "arima"
    
    do.log <- FALSE
    do.cleaning <- FALSE
    do.deseasonalization <- FALSE
    do.scaling <- FALSE
  }
 
}else{
  
  for(i in 1:length(args)){
    print(args[[i]])
  }
  
  experiment <- args[[1]]
  add.bias   <- as.logical(args[[2]])
  algobf     <- args[[3]]
  
  do.log <- FALSE
  do.cleaning <- FALSE
  do.deseasonalization <- FALSE
  do.scaling <- FALSE

}

print(base::date())
print(experiment)

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
library(igraph)
library(MASS)

library(MRCE)

print(experiment)
print(add.bias)
print(algobf)

if(grepl("small", experiment) || grepl("large", experiment)){
  if(grepl("small", experiment)){
    idjob <- unlist(strsplit(experiment, "small-"))[2]
  }else{
    idjob <- unlist(strsplit(experiment, "large-"))[2]
  }
  print(idjob)
  set.seed(idjob)
}

do.usetrainingdata <- TRUE
redo.basef <- FALSE

use.intercept <- FALSE
do.standardize <- ifelse(experiment %in% c("small", "large"), TRUE, FALSE)

mc.cores.basef <-   1 #6          #30 # 20
mc.cores.methods <- 1 #1          #4 #10 # mc.cores.basef  AND mc.cores.methods
nb.cores.cv <-      1 #6          #3

tag <- paste(do.log, "_", do.deseasonalization, "_", do.scaling, sep = "")


H <- 2

#name_methods <- c("REG", "REGBU", "BU", "MINTshr")
#name_methods <- c("L1", "REG", "REG-PBU", "L1-PBU", "BU", "MINTshr", "MINTols")
# KDD -> # name_methods <- c("ERM", "ERMreg", "ERMregbu", "BU", "MINTshr", "MINTsam", "MINTols")
name_methods <- c("MRCE", "ERM")


lambda_selection <- "min"  #"1se"
do.loo <- FALSE
nobiasforwho <- NULL # "agg"  # "bottom" "agg"


#config_main <- list(intercept = FALSE, standardize = FALSE, alpha = .98, thresh = 10^-6, nlambda = 50)
config_main <- list(intercept = use.intercept, standardize = do.standardize, alpha = 1, thresh = 10^-6, nlambda = 50)
config_cvglmnet <- c(config_main, list(nfolds = 3))



if(algobf == "arima"){
#  config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
#                        param_fit_fct = list(seasonal = TRUE, ic = "aic",  max.p = 1, max.q = 0, max.P = 1, max.Q = 0,
#                                           approximation = TRUE, stationary = FALSE, num.cores = 2), 
#                        param_forecast_fct = list(use.initial.values = TRUE))
  
  config_forecast <- list(fit_fct = auto.arima, forecast_fct = Arima, 
                          param_fit_fct = list(seasonal = FALSE, ic = "aic", max.p = 2, max.q = 2,
                                               approximation = TRUE, stationary = FALSE, num.cores = 2), 
                          param_forecast_fct = list(use.initial.values = TRUE))
}else if (algobf == "ets"){
  config_forecast <- list(fit_fct = ets, forecast_fct = ets, 
                        param_fit_fct = list(), 
                        param_forecast_fct = list(use.initial.values = FALSE))
}
config_forecast_agg <- config_forecast_bot <- config_forecast

#######################################################################################
if(grepl("tourism", experiment)){
  ##X <- read.csv("../data/Tourism data_v3.csv")
  ##Z <- X[, -seq(2)]
  ##obj_sort <- sort(apply(Z == 0, 2, sum), index = T, decreasing = T)
  ##Z <- Z[, -head(obj_sort$ix, 34)]

  myfile <- file.path(rdata.folder, paste(experiment, ".Rdata", sep = ""))
  load(myfile)

  Z <- ts(Z, freq = 12) ## IMPORTANT FOR CLEANING
  y <- hts(Z, characters = c(1, 1, 1, 4))
  refit_step <- 12
  
  ##if(TRUE){
    ##series_name <- substr(colnames(Z), 1, 3)
    ##series_name <-sapply(seq(length(series_name)), function(j){
    ##  paste(series_name[j], sprintf("%04d", j), sep = "")
    ##})
    ##colnames(Z) <- series_name
  ##  y <- hts(Z, characters = c(1, 1, 1, 4))
  ##}else{
  ##  y <- gts(Z, characters = list(c(1, 1, 1), c(3)))
  ##}

}else if(grepl("wikipedia", experiment)){

  myfile <- file.path(rdata.folder, paste(experiment, ".Rdata", sep = ""))
  load(myfile)
  
  Z <- ts(Z, freq = 7)  ## IMPORTANT FOR CLEANING
  y <- hts(Z, characters = c(2, 3, 3, 4) )
  refit_step <- 14
}else if(grepl("road_traffic", experiment)){
  myfile <- file.path(rdata.folder, paste(experiment, ".Rdata", sep = ""))
  load(myfile)
  
  Z <- ts(Z, freq = 7)  ## IMPORTANT FOR CLEANING
  y <- hts(Z, characters = c(1, 1, 4) )
  refit_step <- 14
}else if(grepl("elec", experiment)){
  myfile <- file.path(rdata.folder, paste(experiment, ".Rdata", sep = ""))
  load(myfile)
  
  Z <- ts(Z, freq = 7)  ## IMPORTANT FOR CLEANING
  y <- hts(Z, characters = c(1, 1, 4) )
  refit_step <- 14
}  


if(grepl("tourism", experiment)){
  # PREVIOUSLY
  #T_train <- 9 * 12 #108 # 96 #7
  #T_valid <- 5 * 12 #60 # 5
  
  T_train <- 5 * 12
  T_valid <- 9 * 12
  T_test <- 60
}else if(grepl("wikipedia", experiment)){
  T_train <- 86
  T_valid <- 160
  T_test <- 120
}else if(grepl("road_traffic", experiment)){
  T_train <- 7*17 + 1
  T_valid <- 7*17 + 1
  T_test <- 7*18
}else if(grepl("elec", experiment)){
  T_train <- 7*17 
  T_valid <- 7*17 + 1
  T_test <- 7*18
}

if(grepl("tourism", experiment) || grepl("wikipedia", experiment) 
   || grepl("road_traffic", experiment) || grepl("elec", experiment)){
  vec <- sapply(y$nodes, length)
  niveaus <- unlist(sapply(seq(length(vec)), function(i){
    rep(i, vec[i])
  }))
  niveaus <- c(niveaus, rep(length(vec) + 1, ncol(y$bts)))
  S <- smatrix(y)
  A <- head(S, nrow(S) - ncol(Z))
  #######################################################################################
  if(do.log){
    bts <- log(1 + Z)
  }else{
    bts <- Z
  }
}


if(do.cleaning){
  istart <- 1
  iend <- T_train + T_valid + T_test
  
  #iend <- T_train + T_valid
  #if(grepl("road_traffic", experiment) || grepl("elec", experiment)){
  #iend <- T_train + T_valid + T_test
  #}
  #idx <- seq(istart, iend)
  obj_clean <- lapply(seq(ncol(bts)), function(j){
    #x[idx] <- tsclean(subset(x, start = istart, end = iend))
    #x
    x <- bts[, j]
    outliers <- tsoutliers(x)
    x[outliers$index] <- outliers$replacements

    x <- na.interp(x)
    list(outindex = outliers$index, x = x)
  })
  cleaned_bts <- sapply(obj_clean, "[[", "x")
  outindex <- lapply(obj_clean, "[[", "outindex")
  bts <- cleaned_bts

}

if(do.deseasonalization){
  yts <- t(S %*% t(bts))
  
  if(grepl("tourism", experiment)){
    yts <- ts(yts, c(1998, 1), freq = 12)
  }else if(grepl("wikipedia", experiment) || grepl("road_traffic", experiment) || grepl("elec", experiment)){
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

if(grepl("tourism", experiment)){
  bts <- ts(bts, c(1998, 1), freq = 12)
}else if(grepl("wikipedia", experiment) || grepl("road_traffic", experiment) || grepl("elec", experiment)){
  bts <- ts(bts, freq = 7)
}


if(grepl("tourism", experiment) || grepl("wikipedia", experiment) 
   || grepl("road_traffic", experiment) || grepl("elec", experiment)){
  datasave_file <- file.path(results.folder, paste("datasave_", experiment, "_", tag ,".Rdata", sep = ""))
  save(file = datasave_file, list = c("Z", "bts", "cleaned_bts", "outindex", "T_train", "T_valid", "T_test"))
}



if(grepl("small", experiment) || grepl("large", experiment)){
  refit_step <- 40
  T_learn <- 600
  T_train <- floor((2 * T_learn)/3)
  T_valid <- T_learn - T_train
  T_test <- 200
  T_all <- T_train + T_valid + T_test
  n_warm <- 300
  n_simul <- T_all
  
  if(grepl("small", experiment)){
    A <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1))
    obj_simul <- simulate_hts(n_simul)
    bts <- obj_simul$bts
    niveaus <- c(1, 2, 2, 3, 3, 3, 3)
  }else if(grepl("large", experiment)){
    obj_simul <- list()
    obj_simul$param <- NULL
    obj_simul$param$ar_param <- obj_simul$param$ma_param <- NULL
    
    ngroups <- 25
    list_params <- vector("list", ngroups)
    
    res <- NULL
    for(igroup in seq(ngroups)){
      NM <- cbind(paste("T", igroup, sep = "") , 
                  c( paste("A", igroup, sep = "") , paste("B", igroup, sep = "") , 
                     paste("C", igroup, sep = "") , paste("D", igroup, sep = "")   ))
      res <- rbind(res, NM)
    }
    res <- cbind("T", res)
    
    tags <- cbind(res[, 1], sapply(seq(2, ncol(res)), function(j)(
      apply(res[, seq(j)], 1, paste, collapse = "")
    )))
    
    myinfo <- makeINFO2(tags)
    A <- myinfo$A
    
    bts <- matrix(NA, nrow = n_simul, ncol = ngroups * 4)
    for(igroup in seq(ngroups)){
      res_sim <- simulate_hts(n_simul)
      obj_simul$param$ar_param <- c(obj_simul$param$ar_param, res_sim$param$ar_param)
      obj_simul$param$ma_param <- c(obj_simul$param$ma_param, res_sim$param$ma_param)
      id <- seq((igroup - 1) * 4 + 1, (igroup - 1) * 4 + 4 )
      bts[, id] <- res_sim$bts
    }
    niveaus <- c(1, rep(2, 25), rep(3, 100))
  }
  else if(grepl("largeOLD", experiment)){
    res <- simulate_large_hts(n_simul)
    A <- res$A
    btsWithWarm <- res$bts
    bts <- tail(btsWithWarm, -n_warm)
  }
  S <- rbind(A, diag(ncol(bts)))
}

if(experiment == "small-1"){
  myfiledata <- file.path(rdata.folder, paste(experiment, ".Rdata", sep = ""))
  Z <- bts
  save(file = myfiledata, list = c("Z"))
}

if(grepl("tourism", experiment) || grepl("wikipedia", experiment) 
   || grepl("road_traffic", experiment) || grepl("elec", experiment) ){
  T_learn <- T_train + T_valid
  T_all <- T_learn + T_test
}
stopifnot(T_all == nrow(bts))

my_bights <- bights(bts, A)


exp <- unlist(strsplit(experiment, "-"))[1]
info_file <- file.path(results.folder, paste("info_", exp, ".Rdata", sep = ""))
save(file = info_file, list = c("A", "niveaus"))

print(paste(" m = ", my_bights$nbts, " - n = ", my_bights$nts, sep = ""))

file_bf <- file.path(bf.folder, 
                     paste("bf_", experiment, "_", refit_step, "_", tag, "_", algobf, ".Rdata", sep = ""))  


if(file.exists(file_bf) && !redo.basef){
  load(file_bf)
  print("loaded")
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

  
  save(file = file_bf, list = c("data_valid", "data_test", "list_subsets_valid", "list_subsets_test"))
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
  
  stop("CAREFUL - CHECK IF BIAS IS ADDED PROPERLY WHEN USING WHOLE TRAINING IN LEARNING")
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

P_BU <- Matrix(0, nrow = my_bights$nbts, ncol = my_bights$nts, sparse = T)
P_BU[cbind(seq(my_bights$nbts), seq(my_bights$nts - my_bights$nbts + 1, my_bights$nts)) ] <- 1
P_OLS <- solve(t(S) %*% S) %*% t(S)
P_0 <- Matrix(0, nrow = my_bights$nbts, ncol = my_bights$nts, sparse = T)

config <- list(glmnet = config_main, cvglmnet = config_cvglmnet)

h <-1 
Yhat_valid_h <- t(Yhat_valid_allh[h, , ])
Y_valid_h    <- t(Y_valid_allh[h, , ])
Yhat_test_h  <- t(Yhat_test_allh[h, , ])

if(do.usetrainingdata)
{
  stopifnot(h == 1)
  
  Yhat_valid_h <- rbind(data_valid$IN_yhat, Yhat_valid_h)
  Y_valid_h    <- rbind(data_valid$IN_y, Y_valid_h)
}

objreg_valid <- list(y = makey(Y_valid_h), X = makeX(Yhat_valid_h, my_bights), Y = Y_valid_h, Yhat = Yhat_valid_h)
objreg_test <- list(X = makeX(Yhat_test_h, my_bights), Yhat = Yhat_test_h)

nValid <- nrow(Yhat_valid_h)


foldid_joint <- make_foldid(nValid, my_bights$nts, config$cvglmnet$nfolds)
config_joint <- config
config_joint$cvglmnet$foldid <- foldid_joint


nfolds <- config$cvglmnet$nfolds
nbperfold <- floor(nValid/nfolds)
foldid_reg    <- rep(seq(nfolds), each = nbperfold)
remainder <- nValid %% nfolds
foldid_reg    <- c(foldid_reg, rep(nfolds, remainder))
config_reg <- config
config_reg$cvglmnet$foldid <- foldid_reg


results <- mclapply(name_methods, function(current_method){
  #print(paste(base::date(), " - ", current_method, sep = ""))
  obj_method <- NULL
  if(current_method == "ERMregbu"){
    obj_method <- list(algo = "ERMreg", Ptowards = pbu(my_bights), config = config_reg, selection = lambda_selection)
  }else if(current_method == "ERMreg"){
    obj_method <- list(algo = "ERMreg", Ptowards = P_0, config = config_reg, selection = lambda_selection)
  }else if(current_method == "BU"){
    obj_learn <- bu(my_bights)
  }else if(grepl("MINT", current_method)){
    wmethod <- substr(current_method, 5, 7)
    cov_method <-  switch(wmethod, ols = "ols", shr = "shrink", sam = "sample", NA)
    e_residuals <- switch(wmethod, ols = NULL, shr = Eresiduals, sam = Eresiduals, NA)
    obj_learn <- mint(my_bights, method = cov_method, e_residuals = e_residuals , h = h) # J U and W
  }else if(current_method == "ERM"){
    obj_learn <- erm(objreg_valid, my_bights) 
  }else if(current_method == "L1"){
    obj_method <- list(algo = "LASSO", Ptowards = NULL, config = config_joint, selection = lambda_selection)
  }else if(current_method == "L1-POLS"){
    obj_method <- list(algo = "LASSO", Ptowards = POLS, config = config_joint, selection = lambda_selection)
  }else if(current_method == "L1-PBU"){
    obj_method <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config_joint, selection = lambda_selection)
  }else if(current_method == "L2" || current_method == "L2-PBU"){
    config_ridge <- config_joint
    config_ridge$glmnet$alpha <- 0
    config_ridge$cvglmnet$alpha <- 0 # !!!!!!!
    if(current_method == "L2"){
      obj_method <- list(algo = "LASSO", Ptowards = NULL, config = config_ridge, selection = lambda_selection)
    }else{
      obj_method <- list(algo = "LASSO", Ptowards = pbu(my_bights), config = config_ridge, selection = lambda_selection)
    }
  }else if(current_method == "MRCE"){
    config_ridge <- config_joint
    config_ridge$glmnet$alpha <- 0
    config_ridge$cvglmnet$alpha <- 0 # !!!!!!!
    obj_method <- list(algo = "MRCE", Ptowards = pbu(my_bights), config = config_ridge, selection = lambda_selection)
  }else{
    obj_learn <- NULL
  }
  #}else if(current_method == "G-L1-PBU"){
  #  obj_method <- list(algo = "GGLASSO", Ptowards = pbu(my_bights), config = config_gglasso, selection = lambda_selection)
  #}
  
  if(!is.null(obj_method)){
      obj_learn <- new_learnreg(objreg_valid, my_bights, obj_method, TRUE, TRUE)
      print("learning finished.")
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


# results_allh[[h]] <-  results ONLY NEED PREDICTIONS + SAVE SPACE
results_h1 <- lapply(results , "[[", "predictions")
allP <- lapply(lapply(results , "[[", "obj_learn"), "[[", "P")

print(names(results_h1))

myfile <- file.path(results.folder, paste("resultsicml_", experiment, "_", lambda_selection, "_", 
                                          add.bias, "_", tag, "_", algobf, ".Rdata", sep = ""))
Y_test_h <- t(Y_test_allh[h, , ])
save(file = myfile, list = c("results_h1", "Y_test_h", "allP"))


print(paste("Finished ", base::date(), sep = ""))



