rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
#######
color_methods <- c("orange", "yellowgreen", "purple", "deeppink", "pink", "yellow",
                   "darkseagreen1", "darkblue", "darkgreen", "red", "cyan", "blue", "aquamarine", "darkseagreen1", "darkseagreen1", "darkseagreen1", "darkseagreen1", 
                   "slategray4", "slategray", "aquamarine")
nb_methods <- length(color_methods)
#######

#methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "MINTsam",  "ERM", "L2-PBU", "L1-PBU", "G-L1-PBU")
#methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "MINTsam",  "ERM", "L2-PBU", "L1-PBU", "L1", "L2")

methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "MINTsam",  "ERM", "L2-PBU", "L1-PBU")
#methods_toprint <- c("BU", "BASE", "BASE2", "MINTshr", "MINTols", "ERM", "L2-PBU", "L1-PBU")


do.ratio <- FALSE
tag <- "NIPS"
experiment <- "large-biased"
id_jobs <- 1 # seq(2000, 2060) #1986 #seq(400, 450) #420 #seq(200, 210) #420 #c(200, 210)

#nb_simulations <- 100 #500
# ids_simulations <- seq(nb_simulations) # 37

ids_simulations <- seq(50)

lambda_selection <- "1se" # "min"

res_experiment <- unlist(strsplit(experiment, "-"))
DGP <- res_experiment[1]

info_file <- file.path(results.folder, paste("info_", DGP, "_", id_jobs[1], "_", lambda_selection, ".Rdata", sep = ""))
load(info_file)

nbfiles <- 0
errors_simulations <- vector("list", length(id_jobs) * length(ids_simulations))

for(idjob in id_jobs){
  files_available <- NULL
  for(i in ids_simulations){
    if(i %% 50 == 0) print(i)
    myfile <- file.path(results.folder, paste("results_", experiment, "_", idjob, ".", i, "_", lambda_selection, ".Rdata", sep = ""))
    
    if(file.exists(myfile)){
      load(myfile) 
      H <- length(results_allh)
      
      do.proceed <- !any(sapply(seq(H), function(h){
        any(sapply(results_allh[[h]], class) == "try-error")
      }))
      
      if(do.proceed){
        nbfiles <- nbfiles + 1
        files_available <- c(files_available, i)
        
        errors <- lapply(seq(H), function(h){
          #print(h)
          FUTURE <- t(Y_test_allh[h, , ])
          res <- sapply(results_allh[[h]], function(results_method){
            #(results_method$predictions - FUTURE)^2
            (results_method - FUTURE)^2
          }, simplify = "array")
          
          res
          #PRED <- simplify2array(lapply(results_allh[[h]], "[[", "predictions"))
        })
        
        errors_simulations[[nbfiles]] <-  errors
      }
    }
  }
}

print(nbfiles)

errors_all <- errors_simulations[seq(nbfiles)] ## seq() !!!!!!!!!!!!!

errors <- sapply(errors_all, function(errors_simulation){
  errors_h <- sapply(errors_simulation, function(obj){
    apply(obj, c(2, 3), mean)
  }, simplify = "array") # 5 horizons
  errors_h
}, simplify = "array")

par(mfrow = c(1, 2))
plot(apply(errors, c(2, 4), mean)[1, ])
#which(apply(errors, c(2, 4), mean)[1, ] > 60)
#browser()

if(experiment == "small-unbiased"){
  simulation_outliers <- c(49, 145)
  errors <- errors[, , , -simulation_outliers]
}
plot(apply(errors, c(2, 4), mean)[1, ])



#v <- apply(res, c(1, 2, 3), mean)
#err <- t(apply(v, c(2, 3), sum))



name_methods <- dimnames(errors)[[2]]
id.keep <- match(methods_toprint, name_methods)
id.base <- match("BASE2", name_methods)

errors <- errors[, id.keep , , ]


myfile <- paste(tag, "_", ifelse(do.ratio, "ratio", "absolute"), "_", experiment, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.9, width = 21 * 0.9)
par(mfrow = c(5, 3))
par(cex.axis=1)


naggts <- nrow(A)
nbts <- ncol(A)
nts <- naggts + nbts
groupings <- list(all = rep(1, nts), agg = c(rep(1, naggts), rep(0, nbts)), bot = c(rep(0, naggts), rep(1, nbts)))

for(h in seq(H)){
  for(i in seq_along(groupings)){
    mygroup <- which(groupings[[i]] == 1)
    
    errors_grouping <- apply(errors[mygroup, , , ], c(2, 3, 4), sum)
    
    mse_mean <- t(apply(errors_grouping, c(1, 2), mean))
    mse_sd   <- t(apply(errors_grouping, c(1, 2), std.error))
    
    mse_mean_h <- mse_mean[h, ]
    mse_sd_h <- mse_sd[h, ]
    
    
    #browser()
    plotCI(mse_mean_h, uiw = mse_sd_h, liw = mse_sd_h, main = paste("Horizon h = ", h, sep = ""), ylab = "MSE", xaxt = "n", xlab = "")
    axis(1, at = seq(length(mse_mean_h)) , labels = methods_toprint, col.axis="blue", las=2)
    #browser()
    
    #boxplot(err_toplot, outline = T, main = paste("Total - ", nbfiles), cex = .5, col = color_methods[id.keep])
  }
}
dev.off()
