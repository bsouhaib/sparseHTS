rm(list = ls())
source("config_paths.R")
source("nicefigs.R")
library(plotrix)
#######

hierarchy_name <- "UKF16"
filehtsplot <- paste(hierarchy_name, "-graph")
#savepdf(file.path(pdf.folder, filehtsplot), height = 26 * 0.2, width = 21 * 0.5)
savepdf(file.path(pdf.folder, filehtsplot))
sparsehts.work.folder <- file.path("/home/rstudio/sparseHTS", "work")
load(file.path(sparsehts.work.folder, paste("myinfo_", hierarchy_name,".Rdata", sep = "")))
plot(itree, layout = layout.reingold.tilford(itree, root=1, circular=T), vertex.size=0, vertex.label=NA, edge.arrow.size=0)
dev.off()

# hierarchy_name <- "UKF13"
# hierarchy_name <- "UKF15"
# hierarchy_name <- "UKF16"
# hierarchy_name <- "UKF21"
# hierarchy_name <- "UKF23"

inithierarchy_names <- c("UKF21", "UKF13", "UKF16", "UKF15", "UKF23")

#list_plots <- list(first = seq(5), second = 4)
list_plots <- list(first = seq(5))

#
hselected <- seq(48)

for(ilist in seq_along(list_plots)){

  hierarchy_names <- inithierarchy_names[list_plots[[ilist]]]
  
  
tag <- "METERS-NIPS"
myfile <- paste(tag, "_", ilist, sep = "")

if(ilist == 1){
  savepdf(file.path(pdf.folder, myfile), height = 26 * 0.5, width = 21 * 0.9)
  par(mfrow = c(2, 3))
}else{
  savepdf(file.path(pdf.folder, myfile))
  par(mfrow = c(1, 1))
}

mymat <- NULL
for(hierarchy_name in hierarchy_names){
  print(hierarchy_name)
  
  methods_toprint <- c("BU", "BASE", "MINTshr", "MINTols", "L1-PBU")
  
  
  lambda_selection <- "1se"
  #lambda_selection <- "min"
  
  
  info_file <- file.path(results.folder, paste("info_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
  load(info_file)
  
  
  myfile <- file.path(results.folder, paste("results_", hierarchy_name, "_", lambda_selection, ".Rdata", sep = ""))
  if(!file.exists(myfile)){
    stop("FILE DOES NOT EXIST")
  }
  load(myfile) 
  
  H <- length(results_allh)
  errors <- lapply(seq(H), function(h){
    FUTURE <-t(Y_test_allh[h, , ])
    
    res <- sapply(results_allh[[h]], function(results_method){
      (results_method - FUTURE)^2
    }, simplify = "array")
    
    res
  })
  
  errors_all <- errors
  
  name_methods <- dimnames(errors_all[[1]])[[3]]
  id.keep <- match(methods_toprint, name_methods)
  
  naggts <- nrow(A)
  nbts <- ncol(A)
  nts <- naggts + nbts
  
  
  #listhorizons <- list(night = seq(1, 12), day = seq(13, 48))
  #listhorizons <- list(night = seq(1, 12), day = seq(13, max(hselected) ))
  #listhorizons <- list(a = seq(1, 12), b = seq(13, 24) , c = seq(25, 36), d = seq(37, max(hselected)))
  listhorizons <- list(all = seq(1, max(hselected)))
  
  
  #listgroups <- list(agg = seq(naggts) , bot = seq(naggts + 1, nts))
  #listgroups <- list(all = seq(nts), agg = seq(naggts) , bot = seq(naggts + 1, nts))
  listgroups <- list(all = seq(nts))
  
  for(ihgroup in seq_along(listhorizons)){
    for(i in seq_along(listgroups)){
      mygroup <- listgroups[[i]]
      myhorizons <- listhorizons[[ihgroup]]
      
      m <- sapply(errors_all[hselected], function(v){ apply(v[, mygroup, id.keep], c(1, 3), sum)   }, simplify = "array")
      
      #res <- apply(m[, , myhorizons], c(1, 2), mean) # avg horizon
      #mse_mean <- apply(res, 2, mean)
      #mse_std <- apply(res, 2, std.error)
      #plotCI(mse_mean, uiw = mse_std, liw = mse_std, ylab = "MSE", xaxt = "n", xlab = "")
      #axis(1, at = seq(length(mse_mean)) , labels = name_methods, col.axis="blue", las=2)
     
      mse_mean <- apply(m[, , myhorizons], 2, mean)
      mse_std <- apply(m[, , myhorizons], 2, function(obj){std.error( c(obj) )})
      plotCI(mse_mean, uiw = mse_std, liw = mse_std, ylab = "MSE", xaxt = "n", xlab = "", 
             main = paste("Hierarchy ", 
                          which(hierarchy_name == inithierarchy_names), " (", nbts, ", ", naggts - 1,  ", ",  1, ")", sep = ""))
      axis(1, at = seq(length(mse_mean)) , labels = methods_toprint, col.axis="blue", las=2)
      
      vec <- paste(format(mse_mean, digits = 3), " (", format(mse_std, digits = 3), ")" , sep = "")
      vec <- c(nbts, naggts, vec)
      mymat <- rbind(mymat, vec)
    }
  }
 
}
dev.off()
colnames(mymat) <- c("$m$", "$k$", methods_toprint)
row.names(mymat) <- hierarchy_names
stargazer(mymat, title = "Results", out = file.path(pdf.folder, "table-meters.tex"))
}


