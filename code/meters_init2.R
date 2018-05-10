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
  }else{
    idseries <- bottomSeries[j - n_agg]
  }
  
  myfile <- file.path(aggseries.folder, paste("series-", idseries, ".Rdata", sep = ""))
  if(!file.exists(myfile)){
    myfile <- file.path(mymeters.folder, paste("mymeter-", idseries, ".Rdata", sep = ""))
  }
  load(myfile)
  
  Ylearn[, j] <- demand[learn$id]
  Ytest[, j]  <- demand[test$id]
}

#pday <- 1
pdays <- seq(48)
list_yhat_valid <- list_yhat_test <- vector("list", length(pdays))
list_y_valid <- list_y_test <- vector("list", length(pdays))

Yhat <- rbind(Xhat_learn, Xhat_test)
Y    <- rbind(Ylearn, Ytest)
ypday <- c(pday_learn, pday_test)

###
id_learn <- c(validation$id, head(test$id, 48 * 60))
id_test <- tail(test$id, 48 * 32)
#id_learn <- c(head(test$id, 48 * 60))
#id_test <- tail(test$id, 48 * 32)
###  

Yhat_learn <- Yhat[id_learn, ]
Y_learn    <- Y[id_learn, ]
hday_learn <- ypday[id_learn]
  
Yhat_test <- Yhat[id_test, ]
Y_test    <- Y[id_test, ]
hday_test <- ypday[id_test]

for(pday in pdays){
  
  if(pday == 48){
    local_pday <- c(pday - 1, pday)
  }else{
    local_pday <- c(pday, pday +1)
  }
  # local_pday <- pday
  
  # 
  ids <- which(hday_learn %in% local_pday)
  Yforecast_learn <- Yhat_learn[ids, ]
  Ytrue_learn <- Y_learn[ids, ]

  #
  ids <- which(hday_test == pday)
  Yforecast_test <- Yhat_test[ids, ]
  Ytrue_test <- Y_test[ids, ]

  Yforecast_learn <- t(Yforecast_learn)
  Ytrue_learn <- t(Ytrue_learn)
  
  Yforecast_test <- t(Yforecast_test)
  Ytrue_test     <- t(Ytrue_test)
  
  list_yhat_valid[[pday]] <- Yforecast_learn
  list_y_valid[[pday]] <- Ytrue_learn
  
  list_yhat_test[[pday]] <- Yforecast_test
  list_y_test[[pday]] <- Ytrue_test
}

data_valid <- NULL
data_valid$Yhat <- aperm(simplify2array(list_yhat_valid), c(3, 1, 2))
data_valid$Y    <- aperm(simplify2array(list_y_valid), c(3, 1, 2))

data_test <- NULL
data_test$Yhat <- aperm(simplify2array(list_yhat_test), c(3, 1, 2))
data_test$Y    <- aperm(simplify2array(list_y_test), c(3, 1, 2))


# Eresiduals
n_past_obs_kd    <- 60 *48
# n_total <- n_series
# R_onestep <- matrix(NA, nrow = length(learn$id) - n_past_obs_kd, ncol = n_total)

mat_residuals <- sapply(c(aggSeries, bottomSeries), function(idseries){
  print(idseries)
  myfile <-  resid_MINT_file <- file.path(insample.folder, "KD-IC-NML", paste("residuals_MINT_", idseries, "_", "KD-IC-NML", ".Rdata", sep = "")) 
  filekd <- file.exists(myfile)	
  if(!filekd){
    resid_MINT_file <- file.path(insample.folder, "DETS", paste("residuals_MINT_", idseries, "_", "DETS", "_", 1, ".Rdata", sep = "")) 
  }
  load(resid_MINT_file)
  
  if(filekd){
    e_vec <- c(rep(NA, n_past_obs_kd), residuals_MINT)
  }else{
    e_vec <- residuals_MINT
  }
  e_vec
})
R_onestep <- tail(mat_residuals, -n_past_obs_kd)

data_test$Eresiduals <- R_onestep

# Tv <- ncol(data_valid$Yhat[1, , ])
# weights_glmnet <- rep(1/apply(R_onestep^2, 2, mean), each = Tv)


source("config_paths.R")
file_bf <- file.path(bf.folder, 
                     paste("bf_", "meters", ".Rdata", sep = "")) 
save(file = file_bf, list = c("data_valid", "data_test"))
