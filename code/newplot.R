rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(ggplot2)
library(reshape2)
plotmat <- function(M, option = c(1, 2)){
  X <- as.matrix(M)
  X <- t(apply(X, 2, rev))
  Xmelted <- melt(X)
  if(option == 1){
    g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= factor(value) )) + geom_tile() + 
      scale_fill_manual(values = c('white','blue')) 
  }else{
    g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= value )) + 
      geom_tile() + 
      #scale_fill_gradient(low="white",high="black")
      scale_fill_gradient()
  }
  
  g + theme_bw()
}

makePlot <- function(M, title = ""){
  Mprime <- M
  Mprime[, ] <- NA
  x <- as.numeric(M)
  if(any(x == 0)){
    x <- x[which(x != 0)]
  }
  mybreaks <- as.numeric(quantile(x, c(seq(0.01, 0.99, 0.2), 1)  ))
  mycols <-  rainbow(length(mybreaks))
  
  for(i in seq_along(mybreaks) ){
    if(i == 1){
      id <- which(M != 0 & M <= mybreaks[i], arr.ind = T)
    }else{
      id <- which(M <= mybreaks[i] &  M > mybreaks[i - 1], arr.ind = T)
    }
    Mprime[id] <- i
  }
  # == 0
  id_zero <- which(M == 0, arr.ind = T)
  if(length(id_zero) > 0)
  {
    Mprime[id_zero] <- 0
    mycols <- c("white", mycols)
  }
  X <- Mprime
  X <- t(apply(X, 2, rev))
  Xmelted <- melt(X)
  g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= factor(value) )) + geom_tile() + 
    scale_fill_brewer(direction = 1, labels = format(mybreaks, scientific = T, digits = 2), 
                      name = "") +
    theme_minimal() +
    xlab("") +
    ylab("") +
    ggtitle(title) + 
    theme(legend.position="none", axis.text = element_text(size=7), 
          plot.title = element_text(size=5))
  g
}

experiment <- "elec1"
algobf <- "arima"
refit_step <- ifelse(experiment == "tourism", 12, 14)

lambda_selection <- "min"
add.bias <- FALSE
do.log <- TRUE
do.deseasonalization <- TRUE
do.scaling <- FALSE
tag <- paste(do.log, "_", do.deseasonalization, "_", do.scaling, sep = "")


if(grepl("elec", experiment) ){
  info_file <- file.path(results.folder, 
                         paste("info_", "elec1" , ".Rdata", sep = ""))
}
load(info_file)

#

nbts <- ncol(A)
naggts <- nrow(A)
nts <- nbts + naggts
S <- rbind(A, diag(nbts))
P_BU <- matrix(0, nrow = nbts, ncol = nts)
P_BU[cbind(seq(nbts), seq(nts - nbts + 1, nts)) ] <- 1

myfile <- file.path(results.folder, paste("resultsicml_", experiment, "_", lambda_selection, "_", 
                                          add.bias, "_", tag, "_", algobf, ".Rdata", sep = ""))
load(myfile)


######
yhat <- results_h1$BASE
ahat <- t(yhat[, seq(naggts)])
bhat <- P_BU %*% t(yhat)
plot.ts(ahat[1, ])
lines(abu[1, ], col = "red")
###########

for(exp in c("small-1", "road_traffic1", "wikipedia1", "elec1")){
  savepdf(paste("tsplot-", exp, sep = ''))
  par(mfrow = c(3, 1),
      oma = c(1.5,1.5,0,0) + 0.1,
      mar = c(0,1.5,1,1) + 0.1)
  #for(j in seq(3)){
    myfile <- file.path(rdata.folder, paste(exp, ".Rdata", sep = ""))
    load(myfile)
    #if(j !=3){
    #  plot.ts(Z[, j], main = "", ylab = "", xaxt = "n", xlab = "")
    #}else{
    #  plot.ts(Z[, j], main = "", ylab = "")
    #}
    if(exp == "small-1")
    {
      X <- Z[seq(300), seq(3)]
    }else{
      X <- Z[, seq(3)]
    }
    colnames(X) <- rep("", ncol(X))
    plot.ts(X, main = "", ylab = "", lwd = 2)
  #}
  dev.off()
}

stop("done")


savepdf(paste(experiment, "_", "allP", sep = ''), width = 21 * 0.95, height = 29.7 * 0.2)
g1 <- makePlot(abs(as.matrix(allP[["MINTshr"]])), "MINTshr")
g2 <- makePlot(abs(as.matrix(allP[["ERM"]])), "ERM")
g3 <- makePlot(abs(as.matrix(allP[["REG"]])), "REG")
g4 <- makePlot(abs(as.matrix(allP[["REGBU"]])), "REGBU")
grid.arrange(g1, g2, g3, g4, nrow = 1)
dev.off()

stop("done")


savepdf(paste(experiment, "_", "allP", sep = ''))
print(plotmat(P_BU, 2))
print(plotmat(abs(allP[["MINTshr"]]), 2))
print(plotmat(abs(allP[["REG"]]), 2))
print(plotmat(abs(allP[["REGBU"]]), 2))
dev.off()

stop("done")


savepdf(paste(experiment, "_", "PMINT", sep = ''))
PMINT <- allP[["MINTshr"]]
print(plotmat(abs(PMINT), 2))
dev.off()

savepdf(paste(experiment, "_", "C-REGBU", sep = ''))
Pmatrix <- allP[["REGBU"]]
C <- Pmatrix - P_BU
print(plotmat( (C != 0)+0, 1))
dev.off()

savepdf(paste(experiment, "_", "C-REG", sep = ''))
PREG <- allP[["REG"]]
C <- PREG - P_BU
print(plotmat( (C != 0)+0, 1))
dev.off()



M <- as.matrix(Pmatrix %*% S - diag(nbts))
plotmat(M, 2)

#heatmap(C, Colv = NA, Rowv = NA, scale="column", 
#        col = colorRampPalette(brewer.pal(8, "Blues"))(25))



stop("done")

# ETS and ARMA bias

algos <- c("arima", "ets")

matbias <- sapply(algos, function(algo){
  file_bf <- file.path(bf.folder, 
                       paste("bf_", experiment, "_", refit_step, "_", tag, "_", algo, ".Rdata", sep = ""))  
  load(file_bf)
  Yhat_test_allh  <- data_test$Yhat
  Y_test_allh    <- data_test$Y
  yhat  <- t(Yhat_test_allh[1, , ])
  ytrue <- t(Y_test_allh[1, , ])
  
  apply((ytrue - yhat)/ytrue, 2, mean) 
})
matbias[258, ] <- c(0, 0)
matbias[379, ] <- c(0, 0)

matplot(matbias, type = 'l')




myindex <- index_bottom