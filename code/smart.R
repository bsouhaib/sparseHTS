rm(list = ls())

source("../../PROJ/code/config_paths.R")
source("../../PROJ/code/config_general.R")
source("../../PROJ/code/config_splitting.R")
#source("../../PROJ/code/jasa_utils.R")
source("../../PROJ/code/utils.R")

load(file.path(work.folder, "myinfo.Rdata"))

#stop("RUN AGAIN WITH RECENT FORECASTS")

res_file <- file.path(work.folder, "revisedf", paste("revised_meanf_Xmatrix_", algo.agg, "_", algo.bottom, ".Rdata", sep = ""))
load(res_file)
# c("Xhat_learn", "Xhat_test", "pday_learn", "pday_test")

Xhat_learn <- Xhat_learn[, -which(colnames(Xhat_learn) %in% c("C1", "S1"))]
Xhat_test  <- Xhat_test[, -which(colnames(Xhat_test) %in% c("C1", "S1"))]

# length(train$id)
# length(valid$id)
# length(learn$id)
# length(test$id)

pday <- 1

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

# Ylearn AND Xhat_learn (learn)
myids <- c(head(validation$id, 1) - seq(48 * 60, 1), validation$id)
ids <- which((learn$id %in% myids) & (pday_learn == pday))
Yforecast <- Xhat_learn[ids, ]
Ytrue <- Ylearn[ids, ]
Yforecast <- t(Yforecast)
Ytrue <- t(Ytrue)

Yhat <- Y <- array(NA, c(1, dim(Yforecast)))
Yhat[1, , ] <- Yforecast
Y[1, , ] <- Ytrue
data_valid <- NULL
data_valid$Yhat <- Yhat
data_valid$Y <- Y

# Ytest AND Xhat_test (test)
ids <- which(pday_test == pday)
Yforecast <- Xhat_test[ids, ]
Ytrue <- Ytest[ids, ]

Yhat <- Y <- array(NA, c(1, dim(Yforecast)))
Yhat[1, , ] <- Yforecast
Y[1, , ] <- Ytrue
data_test <- NULL
data_test$Yhat <- Yhat
data_test$Y <- Y

# Eresiduals
n_past_obs_kd <- 
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


