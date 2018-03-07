
simulte_hts <- function(n_simul){
  p <- 2
  d <- 0
  q <- 1
  phi_2 <- runif(1, min = 0.5, max = 0.7)
  phi_1 <- runif(1, min = phi_2 - 1, max = 1 - phi_2)
  PHI <- c(phi_1, phi_2)
  theta_1 <- runif(1, min = 0.5, max = 0.7)
  THETA <- theta_1
  
  mus <- rep(0, 4)
  Sigma <- rbind(c(5, 3, 2, 1), c(3, 4, 2, 1), c(2, 2, 5, 3), c(1, 1, 3, 4))
  Ematrix <- mvrnorm(n = n_simul, mus, Sigma = Sigma)
  
  
  bts <- NULL
  for(j in seq(4)){
    series <- arima.sim(n = n_simul, list(order = c(p, d, q), ar = PHI, ma = THETA), innov = Ematrix[, j])
    bts <- cbind(bts, series)
  }
  bts <- tail(bts, -n_warm)
}
