
simulate_hts <- function(n_simul){
  p <- 2
  d <- 0
  q <- 1
  if(FALSE){
    phi_2 <- runif(1, min = 0.5, max = 0.7)
    phi_1 <- runif(1, min = phi_2 - 1, max = 1 - phi_2)
    theta_1 <- runif(1, min = 0.5, max = 0.7)
  }else{
    phi_2 <- runif(1, min = 0.5, max = 0.7)
    phi_1 <- runif(1, min = phi_2 - 0.71, max = 0.71 - phi_2)
    theta_1 <- runif(1, min = 0.5, max = 0.7)
  }
  
  PHI <- c(phi_1, phi_2)
  THETA <- theta_1
  
  mus <- rep(0, 4)
  
  if(FALSE){
    Sigma <- rbind(c(5, 3, 2, 1), c(3, 4, 2, 1), c(2, 2, 5, 3), c(1, 1, 3, 4))
  }else{
    #varVec <- rep(1, 4) 
    varVec <- rep(2, 4)
    corMatB <- rbind(c(1, 0.7, 0.2, 0.3), 
                   c(0.7, 1, 0.3, 0.2),
                   c(0.2, 0.3, 1, 0.6),
                   c(0.3, 0.2, 0.6, 1))
    Sigma <- as.matrix(Diagonal(x = sqrt(varVec)) %*% corMatB %*% Diagonal(x = sqrt(varVec)))
  }
  
  Ematrix_insample <- mvrnorm(n = n_simul, mus, Sigma = Sigma)
  Ematrix_start <- mvrnorm(n = n_warm, mus, Sigma = Sigma)
  
  bts <- matrix(NA, nrow = n_simul, ncol = 4)
  for(j in seq(4)){
    bts[, j] <- arima.sim(n = n_simul, list(order = c(p, d, q), ar = PHI, ma = THETA), 
    n.start = n_warm, start.innov = Ematrix_start[, j], innov = Ematrix_insample[, j])
  }

  #bts <- tail(bts, -n_warm)
  list(bts = bts, param = list(phi_1 = phi_1, phi_2 = phi_2, theta_1 = theta_1))
}
