

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
    ggtitle(title)
  g
}

g1 <- makePlot(abs(as.matrix(allP[["MINTshr"]])), "MINTshr")
g2 <- makePlot(abs(as.matrix(allP[["ERM"]])), "ERM")
g3 <- makePlot(abs(as.matrix(allP[["REG"]])), "REG")
g4 <- makePlot(abs(as.matrix(allP[["REGBU"]])), "REGBU")

gird.arrange(g1, g2, g3, g4, nrow = 1)

#M1 <- abs(as.matrix(allP[["MINTshr"]]))
M1 <- abs(as.matrix(allP[["REGBU"]]))
M1p <- M1 
M1p[, ] <- NA
x <- as.numeric(M1)
if(any(x == 0)){
  x <- x[which(x != 0)]
}
mybreaks <- as.numeric(quantile(x, c(seq(0.01, 0.99, 0.2), 1)  ))
mycols <-  rainbow(length(mybreaks))

for(i in seq_along(mybreaks) ){
  if(i == 1){
    id <- which(M1 != 0 & M1 <= mybreaks[i], arr.ind = T)
  }else{
    id <- which(M1 <= mybreaks[i] &  M1 > mybreaks[i - 1], arr.ind = T)
  }

  #print(i)
  #M1p[id] <- mycols[i]
  M1p[id] <- i
  
  #print(table(M1p))
  
}
# == 0
id_zero <- which(M1 == 0, arr.ind = T)
if(length(id_zero) > 0)
{
  #M1p[id_zero] <- "white"
  M1p[id_zero] <- 0
  mycols <- c("white", mycols)
}

X <- M1p
X <- t(apply(X, 2, rev))
Xmelted <- melt(X)
g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= factor(value) )) + geom_tile() + 
  scale_fill_brewer(direction = 1, labels = format(mybreaks, scientific = T, digits = 2), 
                    name = "") +
  theme_minimal()
plot(g)


M2 <- abs(as.matrix(allP[["REG"]]))



M <- M1
#x <- c(as.numeric(M1), as.numeric(M2))
#mybreaks <- as.numeric(quantile(x, seq(0, 1, by = .25)))
X <- as.matrix(M)
X <- t(apply(X, 2, rev))
Xmelted <- melt(X)

#colorRampPalette(c("black", "white"))(5)


g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= value )) + geom_tile() + 
  scale_fill_gradient2(low="black", mid = "white", high="blue", midpoint = 0)
plot(g)

g <- ggplot(data = Xmelted, aes(x=Var1, y=Var2, fill= value )) + geom_tile() + 
  scale_fill_gradient(low="white",high="black", breaks = mybreaks)
plot(g)
