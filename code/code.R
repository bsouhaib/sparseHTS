
makeMatrices <- function(obj_bights, list_subsets, H, fmethod_agg, fmethod_bot, refit_step, mc.cores = 1){
  n <- obj_bights$nts
  m <- obj_bights$nbts
  results <- vector("list", n)
  
  algos <- c(rep(fmethod_agg, n-m), rep(fmethod_agg, m))
  #for(j in seq(n)){
  #  results[[j]] <- rolling.forecast(obj_bights$yts[, j], list_subsets, H, algos[j], refit_step = refit_step)
  #}
  #results <- lapply(seq(n), function(j){
  #  rolling.forecast(obj_bights$yts[, j], list_subsets, H, algos[j], refit_step = refit_step)
  #})
  results <- mclapply(seq(n), function(j){
    rolling.forecast(obj_bights$yts[, j], list_subsets, H, algos[j], refit_step = refit_step)
  }, mc.cores = mc.cores)
  
  Yhat <- simplify2array(lapply(results, "[[", "predictions"))
  Y <- simplify2array(lapply(results, "[[", "future"))
  Eresiduals <- simplify2array(lapply(results, "[[", "e_residuals"))
  
  Yhat <- aperm(Yhat, c(2, 3, 1))
  Y    <- aperm(Y, c(2, 3, 1))

  list(Yhat = Yhat, Y = Y, Eresiduals = Eresiduals)
}  


rolling.forecast <- function(series, list_subsets, H, fmethod = c("AR1", "ARIMA", "ETS"), refit_step){

  n_subsets <- length(list_subsets)
  predictions <- future <- matrix(NA, nrow = n_subsets, ncol = H)  
  e_residuals <- NULL
    
  if(fmethod == "ETS"){
    fit_fct <- ets
    forecast_fct <- ets
    list_param <- NULL
  }else if(fmethod == "ARIMA"){
    fit_fct <- auto.arima
    forecast_fct <- Arima
    list_param <- list(seasonal = FALSE, ic = "aic", max.p = 2, max.q = 2,  approximation = TRUE, stationary = FALSE)
  }else if(fmethod == "AR1"){
    fit_fct <- auto.arima
    forecast_fct <- Arima
    list_param <- list(seasonal = FALSE, ic = "aic", max.p = 1, max.q = 0,  approximation = TRUE, stationary = FALSE)
  }
  
  for(i in seq(n_subsets)){
    ts_split <- list_subsets[[i]]
    learn_series <- series[seq(ts_split[1], ts_split[2])]

    if( (i-1) %% refit_step == 0){
      #model <- fit_fct(learn_series, fmethod)
      model <- do.call(fit_fct, c(list(y = learn_series), list_param))
    }else{
      model <- forecast_fct(learn_series, model = model, use.initial.values = TRUE)
    }
    
    if(i == 1){
      e_residuals <- as.numeric(resid(model))
    }
    
    
    predictions[i, ] <- forecast(model, h = H)$mean
    future[i, ] <- series[ts_split[2] + seq(1, H)]
  }
  
  output <- list(future = future, predictions = predictions, e_residuals = e_residuals)
}

#fit_fct <- function(series, fmethod = c("ARIMA", "ETS")){
#  match.arg(fmethod)
#  if(fmethod == "ARIMA"){
#    model <- auto.arima(series, seasonal = FALSE, 
#                        ic = "aic", max.p = 2, max.q = 2, 
#                        approximation = TRUE, stationary = FALSE)
#  }else if(fmethod == "ETS"){
#    model <- ets(series)
#  }
#  model
#}

