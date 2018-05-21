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


tag <- "SIM-NIPS"
id_jobs <- 1 # seq(2000, 2060) #1986 #seq(400, 450) #420 #seq(200, 210) #420 #c(200, 210)

#nb_simulations <- 100 #500
# ids_simulations <- seq(nb_simulations) # 37

ids_simulations <- seq(100)

lambda_selection <- "1se" # "min"

do.file <- TRUE

if(do.file){
myfile <- paste(tag, sep = "")
savepdf(file.path(pdf.folder, myfile), height = 26 * 0.2, width = 21 * 0.9)
par(mfrow = c(1, 4))
par(cex.axis=.7)
}

#myexp <- c("small-unbiased", "small-biased", "large-unbiased", "large-biased")
#myexp <- c("small-unbiased", "small-biased")
myexp <- c("small-unbiased", "large-unbiased", "small-biased", "large-biased")
for(experiment in myexp){
  
  if(experiment == "small-unbiased"){
    methods_toprint <- c("BU", "BASE", "MINTshr", "MINTols", "MINTsam",  "ERM", "L1-PBU")
    mymain <- "UNBIASED/SMALL"
  }else if(experiment == "small-biased"){
    methods_toprint <- c("BASE2", "BU", "BASE", "MINTshr", "MINTols", "MINTsam",  "ERM", "L1-PBU")
    mymain <- "BIASED/SMALL"
  }else if(experiment == "large-unbiased"){
    methods_toprint <- c("BU", "BASE", "MINTshr", "MINTols", "L1-PBU")
    mymain <- "UNBIASED/LARGE"
  }else if(experiment == "large-biased"){
    methods_toprint <- c("BASE2", "BU", "BASE", "MINTshr", "MINTols", "L1-PBU")
    mymain <- "BIASED/LARGE"
  }


  print(experiment)
  
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
            #if(DGP == "small"){
            #  (results_method$predictions - FUTURE)^2
            #}else{
              (results_method - FUTURE)^2
            #}
            
            
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

#par(mfrow = c(1, 2))
#plot(apply(errors, c(2, 4), mean)[1, ])
#which(apply(errors, c(2, 4), mean)[1, ] > 60)
#browser()



#plot(apply(errors, c(2, 4), mean)[1, ])
#stop("done")
print(sort(apply(errors, c(2, 4), mean)[1, ], index = T, decreasing = T))

if(experiment == "small-biased"){
  simulation_outliers <- c(10)
  errors <- errors[, , , -simulation_outliers]
}

#v <- apply(res, c(1, 2, 3), mean)
#err <- t(apply(v, c(2, 3), sum))



name_methods <- dimnames(errors)[[2]]
id.keep <- match(methods_toprint, name_methods)
id.base <- match("BASE2", name_methods)

errors <- errors[, id.keep , , ]


naggts <- nrow(A)
nbts <- ncol(A)
nts <- naggts + nbts
#groupings <- list(all = rep(1, nts), agg = c(rep(1, naggts), rep(0, nbts)), bot = c(rep(0, naggts), rep(1, nbts)))
#groupings <- list(agg = c(rep(1, naggts), rep(0, nbts)), bot = c(rep(0, naggts), rep(1, nbts)))
groupings <- list(all = rep(1, nts))



  for(i in seq_along(groupings)){
    mygroup <- which(groupings[[i]] == 1)
    
    errors_grouping <- apply(errors[mygroup, , , ], c(2, 3, 4), sum)
    errors_grouping <- apply(errors_grouping, c(1, 3), mean)
    
    mse_mean <- t(apply(errors_grouping, c(1), mean))
    mse_std   <- t(apply(errors_grouping, c(1), std.error))
    
    
    #browser()
    plotCI(mse_mean, uiw = mse_std, liw = mse_std, main = mymain, ylab = "MSE", xaxt = "n", xlab = "")
    axis(1, at = seq(length(mse_mean)) , labels = methods_toprint, col.axis="blue", las=2)
    #browser()
    
    #boxplot(err_toplot, outline = T, main = paste("Total - ", nbfiles), cex = .5, col = color_methods[id.keep])
  }

}


if(do.file){
dev.off()
}