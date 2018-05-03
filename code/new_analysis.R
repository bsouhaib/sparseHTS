rm(list = ls())
source("config_paths.R")
source("nicefigs.R")

#######
color_methods <- c("orange", "yellowgreen", "purple", "deeppink", "pink", "yellow",
                   "darkseagreen1", "darkblue", "darkgreen", "red", "cyan", "blue", "aquamarine", "darkseagreen1", "darkseagreen1", "darkseagreen1", "darkseagreen1", 
                   "slategray4", "slategray", "aquamarine")
nb_methods <- length(color_methods)
#######

#methods_toprint <- c("BU", "BASE", "BASE-NOB", "MINT", "MINTOLS",  "LSglmnet", "RID", "LAS", "RID_PBU","LAS_PBU")
methods_toprint <- c("BU", "BASE", "BASE2", "MINTSHR", "MINTOLS", "MINTSAM",  "OLS", "L2", "L1", "L2-PBU", "L1-PBU")

do.ratio <- TRUE
tag <- "nips"
experiment <- "large-biased"
H <- 1
id_jobs <- 1 # seq(2000, 2060) #1986 #seq(400, 450) #420 #seq(200, 210) #420 #c(200, 210)
nb_simulations <- 100 #500
ids_simulations <- seq(nb_simulations) # 37
fmethod_agg <- "ARIMA"
fmethod_bot <- "ARIMA"
lambda_selection <- "1se" # "min"

nbfiles <- 0
results <- vector("list", length(id_jobs) * length(ids_simulations))

for(idjob in id_jobs){
  for(i in ids_simulations){
    if(i %% 50 == 0) print(i)
    myfile <- file.path(rdata.folder, paste("results_", experiment, "_", fmethod_agg, "_", fmethod_bot, "_", idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
    if(file.exists(myfile)){
      nbfiles <- nbfiles + 1
      load(myfile) 
      #print(dim(Ytilde_test_allh))
      #browser()
      res <- sapply(seq(nb_methods), function(imethod){
        apply((Ytilde_test_allh[seq(H), , , imethod] - Y_test_allh[seq(H), , ])^2, c(1, 2), mean)
      }, simplify = "array")
      dimnames(res)[3] <- dimnames(Ytilde_test_allh)[4]
      results[[nbfiles]] <-  res
      
      #i_base <- match("BASE2", dimnames(res)[[3]])
      #res_ratio <- sapply(seq(nb_methods), function(imethod){
      #  (res[, , imethod] - res[, , i_base] )/res[, , i_base]
      #}, simplify = "array")
      #dimnames(res_ratio)[3] <- dimnames(res)[3]
      #results_ratio[[nbfiles]] <-  res_ratio
      
    }
  }
}

err <- sapply(seq(nbfiles), function(ifile){
  mean(results[[ifile]], na.rm = T)
})

res_sort <- sort(err, decreasing = T, index = T)
print(res_sort)
#myresults <- results[-res_sort$ix[seq(8)]]


myresults <- results[seq(nbfiles)] ## seq() !!!!!!!!!!!!!


id.keep <- match(methods_toprint, dimnames(myresults[[1]])[[3]])
id.base <- match("BASE2", dimnames(myresults[[1]])[[3]])

do.log <- FALSE
myfile <- paste(tag, "_", ifelse(do.ratio, "ratio", "absolute"), "_", experiment, "_", fmethod_agg, "_", fmethod_bot, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.9)
par(mfrow = c(4, 2))
par(cex.axis=.4)

v <- sapply(myresults, function(mat){
  if(do.ratio){
    (apply(mat[seq(H), , id.keep], 2,  sum) - sum(mat[seq(H), , id.base]))/sum(mat[seq(H), , id.base])
  }else{
    apply(mat[seq(H), , id.keep], 2,  sum)
  }
})
v <- t(v)
#colnames(v) <- name_methods[id.keep]
final_v <- v
final_medians <- apply(v, 2, median)
if(do.log){
  final_v <- log(final_v)
  final_medians <- log(final_medians)
}
mymin <- min(final_medians)
boxplot(final_v, outline = F, 
        main = paste(fmethod_agg, "-", fmethod_bot, " - ", "Total - ", nbfiles), cex = .5, col = color_methods[id.keep])
abline(h = mymin, col = "red")


for(j in seq(7)){
  v <- sapply(myresults, function(mat){
    if(do.ratio){
      (t(mat[seq(H), j, id.keep]) - mat[seq(H), j, id.base])/mat[seq(H), j, id.base]
    }else{
      t(mat[seq(H), j, id.keep])
    }
  })
  v <- t(v)
  #colnames(v) <- name_methods[id.keep]
  
  final_v <- v
  final_medians <- apply(v, 2, median)
  if(do.log){
    final_v <- log(final_v)
    final_medians <- log(final_medians)
  }
  mymin <- min(final_medians)
  
  boxplot(final_v, outline = F, main = j, col = color_methods[id.keep])
  abline(h = mymin, col = "red")
  #boxplot(v - v[, "BASE"], outline = F)
  #stop("done")
}
dev.off()

stop("done")








###############
res <- simplify2array(results)
res <- apply(res, c(1, 2, 3), mean)
par(mfrow = c(3, 3))
for(j in seq(nb_methods)){
  #matplot(res[, j, ], lty = 1, type = "l")
  matplot(res[, j, ])
}

res_overall <- apply(res, c(1, 3), mean)
matplot(res_overall)

#####
#if(do.ratio){
#  myresults <- results_ratio[seq(nbfiles)]
#}

err_all <- Reduce("+", myresults)/length(myresults)
err_sd <- sqrt(Reduce("+", lapply(myresults, function(mat){ (mat - err_all)^2}))/length(myresults))
err_avgnodes <- apply(err_all, c(1, 3), mean)
#colnames(err_avgnodes) <- name_methods

v <- lapply(myresults, function(mat){ apply(mat, c(1, 3), mean) })
v_mean  <- Reduce("+", v)/length(v)
v_sd <- sqrt(Reduce("+", lapply(v, function(mat){ (mat - v_mean)^2}))/length(v))
####

htoplot <- 3
plot.splitted <- TRUE
if(FALSE){
  if(TRUE){
    par(mfrow = c(3, 3))
    
    barplot(err_avgnodes[1, id.keep], cex.names = .5, main = "average over hts (h =1)")
    
    matplot(log(err_avgnodes[seq(htoplot), id.keep]), type = 'p', pch = 21, col = color_methods[id.keep], main = paste(lambda_selection, " - average over hts"), ylab = "MSE (log)")
    #matpoints(err_avgnodes[, id.keep], col = color_methods[id.keep], type = 'p')
    legend("bottomright", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
    
    for(j in seq(7)){
      
      matplot(log(err_all[seq(htoplot), j, id.keep]), type = 'p', pch = 21, col = color_methods[id.keep], main = j, ylab = "MSE (log)")
      if(j == 1){
        legend("bottomright", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
      }
      #barplot(err_all[1, j, id.keep])
    }
  }else{
    #####
    matplot(err_avgnodes[seq(10), id.keep], type = 'p', pch = 21, col = color_methods[id.keep], main = paste(lambda_selection))
    #matpoints(err_avgnodes[, id.keep], col = color_methods[id.keep], type = 'p')
    legend("topleft", name_methods[id.keep], col = color_methods[id.keep], lty = 1, cex = .5)
  }
  #u <- apply(results[[1]], c(1, 3), mean)
  #matplot(u[seq(3), c(4, 5, 7, 8)])
  # 
}
