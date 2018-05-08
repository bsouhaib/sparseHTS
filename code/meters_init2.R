rm(list = ls())

sparsehts.work.folder <- file.path("/home/rstudio/sparseHTS", "work")
load(file.path(sparsehts.work.folder, "myinfo.Rdata"))

source("../../PROJ/code/config_paths.R")
source("../../PROJ/code/config_general.R")
source("../../PROJ/code/config_splitting.R")
#source("../../PROJ/code/jasa_utils.R")
source("../../PROJ/code/utils.R")



#stop("RUN AGAIN WITH RECENT FORECASTS")

res_file <- file.path(work.folder, "revisedf", paste("revised_meanf_Xmatrix_", algo.agg, "_", algo.bottom, ".Rdata", sep = ""))
load(res_file)
# c("Xhat_learn", "Xhat_test", "pday_learn", "pday_test")

Xhat_learn <- Xhat_learn[, -which(colnames(Xhat_learn) %in% c("C1", "S1"))]
Xhat_test  <- Xhat_test[, -which(colnames(Xhat_test) %in% c("C1", "S1"))]

Xhat_learn <- Xhat_learn[, which(colnames(Xhat_learn) %in% c(aggSeries, bottomSeries))]
Xhat_test  <- Xhat_test[, which(colnames(Xhat_test) %in% c(aggSeries, bottomSeries))]


# length(train$id)
# length(valid$id)
# length(learn$id)
# length(test$id)

######
n_series <- n_agg + n_bottom
Ylearn <- matrix(NA, nrow = length(learn$id), ncol = n_series)
Ytest <- matrix(NA, nrow = length(test$id), ncol = n_series)
for(j in seq(n_series)){
  if(j%%100 ==0)
    print(j)
  if(j <= n_agg){
    idseries <- aggSeries[j]
    load(file.path(aggseries.folder, paste("series-", idseries, ".Rdata", sep = "")))
  }else{
    idseries <- bottomSeries[j - n_agg]
    load(file.path(mymeters.folder, paste("mymeter-", idseries, ".Rdata", sep = "")))
  }
  Ylearn[, j] <- demand[learn$id]
  Ytest[, j]  <- demand[test$id]
}

#pday <- 1
pdays <- seq(10)
list_yhat_valid <- list_yhat_test <- vector("list", length(pdays))
list_y_valid <- list_y_test <- vector("list", length(pdays))

for(pday in pdays){
  # Ylearn AND Xhat_learn (learn)
  myids <- c(head(validation$id, 1) - seq(48 * 60 * 2, 1), validation$id)
  ids <- which((learn$id %in% myids) & (pday_learn == pday))
  Yforecast <- Xhat_learn[ids, ]
  Ytrue <- Ylearn[ids, ]
  Yforecast <- t(Yforecast)
  Ytrue <- t(Ytrue)
  
  #Yhat <- Y <- array(NA, c(1, dim(Yforecast)))
  #Yhat[1, , ] <- Yforecast
  #Y[1, , ] <- Ytrue
  
  list_yhat_valid[[pday]] <- Yforecast
  list_y_valid[[pday]] <- Ytrue
  
  # Ytest AND Xhat_test (test)
  ids <- which(pday_test == pday)
  Yforecast <- Xhat_test[ids, ]
  Ytrue <- Ytest[ids, ]
  Yforecast <- t(Yforecast)
  Ytrue <- t(Ytrue)
  
  #Yhat <- Y <- array(NA, c(1, dim(Yforecast)))
  #Yhat[1, , ] <- Yforecast
  #Y[1, , ] <- Ytrue
  
  list_yhat_test[[pday]] <- Yforecast
  list_y_test[[pday]] <- Ytrue
}

data_valid <- NULL
data_valid$Yhat <- aperm(simplify2array(list_yhat_valid), c(3, 1, 2))
data_valid$Y    <- aperm(simplify2array(list_y_valid), c(3, 1, 2))

data_test <- NULL
data_test$Yhat <- aperm(simplify2array(list_yhat_test), c(3, 1, 2))
data_test$Y    <- aperm(simplify2array(list_y_test), c(3, 1, 2))


# Eresiduals
n_past_obs_kd    <- 60 *48
n_total <- n_series
R_onestep <- matrix(NA, nrow = length(learn$id) - n_past_obs_kd, ncol = n_total)

for(do.agg in c(TRUE, FALSE)){
  if(do.agg){
    set_series <- aggSeries
    algo <- algo.agg
  }else{
    set_series <- bottomSeries
    algo <- algo.bottom
  }
  
  mat_residuals <- sapply(seq_along(set_series), function(j){
    if(j%%100 == 0)
      print(j)
    
    idseries <- set_series[j]
    if(algo == "KD-IC-NML"){
      resid_MINT_file <- file.path(insample.folder, algo, paste("residuals_MINT_", idseries, "_", algo, ".Rdata", sep = "")) 
      load(resid_MINT_file) # residuals_MINT
      e_vec <- c(rep(NA, n_past_obs_kd), residuals_MINT)
    }else if(algo == "DYNREG" || algo == "DETS"){
      resid_MINT_file <- file.path(insample.folder, algo, paste("residuals_MINT_", idseries, "_", algo, "_", 1, ".Rdata", sep = "")) 
      load(resid_MINT_file)
      e_vec <- residuals_MINT
    }
    e_vec
  })
  print(dim(mat_residuals))
  mat_residuals <- tail(mat_residuals, -n_past_obs_kd)
  if(do.agg){
    R_onestep[, seq(n_agg)] <- mat_residuals
  }else{
    R_onestep[, seq(n_agg + 1, n_total)] <- mat_residuals
  }
}

data_test$Eresiduals <- R_onestep

source("config_paths.R")
file_bf <- file.path(bf.folder, 
                     paste("bf_", "meters", ".Rdata", sep = "")) 
save(file = file_bf, list = c("data_valid", "data_test"))
