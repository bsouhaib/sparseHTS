
makeMatrices <- function(obj_bights, list_subsets, H, fmethod, refit_step){
  n <- obj_bights$nts
  results <- vector("list", n)
  for(j in seq(n)){
    #print(j)
    series <- obj_bights$yts[, j]
    results[[j]] <- rolling.forecast(series, list_subsets, H, fmethod, refit_step = refit_step)
  }
  Yhat <- simplify2array(lapply(results, "[[", "predictions"))
  Y <- simplify2array(lapply(results, "[[", "future"))
  
  Yhat <- aperm(Yhat, c(2, 3, 1))
  Y    <- aperm(Y, c(2, 3, 1))

  list(Yhat = Yhat, Y = Y)
}  


rolling.forecast <- function(series, list_subsets, H, fmethod, refit_step){

  n_subsets <- length(list_subsets)
  predictions <- future <- matrix(NA, nrow = n_subsets, ncol = H)  
  
  if(fmethod == "ETS"){
    forecast_function <- ets
  }else if(fmethod == "ARIMA"){
    forecast_function <- Arima
  }
  
  for(i in seq(n_subsets)){
    ts_split <- list_subsets[[i]]
    learn_series <- series[seq(ts_split[1], ts_split[2])]

    if( (i-1) %% refit_step == 0){
      model <- fit_fct(learn_series, fmethod)
    }else{
      model <- forecast_function(learn_series, model = model, use.initial.values = TRUE)
    }
    
    predictions[i, ] <- forecast(model, h = H)$mean
    future[i, ] <- series[ts_split[2] + seq(1, H)]
  }
  
  output <- list(future = future, predictions = predictions)
}
