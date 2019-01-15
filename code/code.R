
makeMatrices <- function(obj_bights, list_subsets, H, config_forecast_agg, config_forecast_bot, refit_step, mc.cores = 1){
  n <- obj_bights$nts
  m <- obj_bights$nbts
  results <- vector("list", n)
  
  #algos <- c(rep(fmethod_agg, n-m), rep(fmethod_agg, m))
  #for(j in seq(n)){
  #  results[[j]] <- rolling.forecast(obj_bights$yts[, j], list_subsets, H, algos[j], refit_step = refit_step)
  #}
  #results <- lapply(seq(n), function(j){
  #  rolling.forecast(obj_bights$yts[, j], list_subsets, H, algos[j], refit_step = refit_step)
  #})
  #print(mc.cores)
  
  results <- mclapply(seq(n), function(j){
    if(j <= n-m){
      config_forecast <- config_forecast_agg
    }else{
      config_forecast <- config_forecast_bot
    }
    rolling.forecast(obj_bights$yts[, j], list_subsets, H, config_forecast, refit_step = refit_step)
  }, mc.cores = mc.cores)
  
  Yhat <- simplify2array(lapply(results, "[[", "predictions"))
  Y <- simplify2array(lapply(results, "[[", "future"))
  Eresiduals <- simplify2array(lapply(results, "[[", "e_residuals"))
  
  IN_y    <- simplify2array(lapply(results, "[[", "insample_y"))
  IN_yhat <- simplify2array(lapply(results, "[[", "insample_yhat"))
  
  Yhat <- aperm(Yhat, c(2, 3, 1))
  Y    <- aperm(Y, c(2, 3, 1))

  list(Yhat = Yhat, Y = Y, Eresiduals = Eresiduals, IN_y = IN_y, IN_yhat = IN_yhat)
}  


rolling.forecast <- function(series, list_subsets, H, config_forecast, refit_step){
  n_subsets <- length(list_subsets)
  predictions <- future <- matrix(NA, nrow = n_subsets, ncol = H)  
  e_residuals <- NULL
  insample_yhat <- insample_y <- NULL
  
  fit_fct <- config_forecast$fit_fct  
  forecast_fct <- config_forecast$forecast_fct
  param_fit_fct <- config_forecast$param_fit_fct
  param_forecast_fct <- config_forecast$param_forecast_fct
  
  #if(fmethod == "ETS"){
  #  fit_fct <- ets
  #  forecast_fct <- ets
  #  list_param <- NULL
  #}else if(fmethod == "ARIMA"){
  #  fit_fct <- auto.arima
  #  forecast_fct <- Arima
  #  list_param <- list(seasonal = FALSE, ic = "aic", max.p = 2, max.q = 2,  approximation = TRUE, stationary = FALSE)
  #}else if(fmethod == "AR1"){
  #  fit_fct <- auto.arima
  #  forecast_fct <- Arima
  #  list_param <- list(seasonal = FALSE, ic = "aic", max.p = 1, max.q = 0,  approximation = TRUE, stationary = FALSE)
  #}
  
  for(i in seq(n_subsets)){
    ts_split <- list_subsets[[i]]
    if(is.ts(series)){
      learn_series <- subset(series, start = ts_split[1], end = ts_split[2])
    }else{
      learn_series <- series[seq(ts_split[1], ts_split[2])]
    }
    
    #browser()

    if( (i-1) %% refit_step == 0){
      #model <- fit_fct(learn_series, fmethod)
      model <- do.call(fit_fct, c(list(y = learn_series), param_fit_fct))
    }else{
      #model <- forecast_fct(learn_series, model = model, use.initial.values = TRUE)
      model <- do.call(forecast_fct, c(list(y = learn_series, model = model), param_forecast_fct))
    }
    
    if(i == 1){
      e_residuals <- as.numeric(resid(model))
      insample_y <- as.numeric(learn_series)
      insample_yhat <- as.numeric(fitted(model))
    }
    
    predictions[i, ] <- forecast(model, h = H)$mean
    future[i, ] <- series[ts_split[2] + seq(1, H)]
  }
  
  output <- list(future = future, predictions = predictions, e_residuals = e_residuals, insample_y = insample_y, insample_yhat = insample_yhat)
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

