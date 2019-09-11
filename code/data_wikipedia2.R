rm(list = ls())
source("config_paths.R")
library(readr)
library(dplyr)
source("nicefigs.R")

myfile <- file.path(rdata.folder, paste("wikipediaDT", ".Rdata", sep = ""))
load(myfile)

# SUBSET
set_series <- seq(71, 100)
mylist <- vector("list", length(set_series))
for(myhts in set_series){
  set.seed(500 + myhts)
  
  #nb <- 50
  nb <- 25
  i_en <- sample(which(meta_DT$lang == "en"), nb)
  i_fr <- sample(which(meta_DT$lang == "fr"), nb)
  i_de <- sample(which(meta_DT$lang == "de"), nb)
  i_ja <- sample(which(meta_DT$lang == "ja"), nb)
  i_ru <- sample(which(meta_DT$lang == "ru"), nb)
  i_zh <- sample(which(meta_DT$lang == "zh"), nb)
  i_all <- c(i_en, i_fr, i_de, i_ja, i_ru, i_zh)
  mymetaDT <- meta_DT[i_all, ]
  
  mylist[[myhts]] <- mymetaDT %>% group_by(lang, access) %>% count() %>% .$n
  
  
  tempDT <- DT[i_all, ]
  myDT <- t(tempDT[, -1])
  
  res <- sapply(seq(nrow(mymetaDT)), function(irow){
    paste(mymetaDT[irow, "lang"], 
          mymetaDT[irow, "access"], 
          mymetaDT[irow, "agent"], sprintf("%04d", irow), sep = "")
  })
  
  colnames(myDT) <- res
  Z <- myDT
  Z <- tail(Z, -184)
  
  # with set.seed(1466)
  #set.outliers <- c(185L, 277L, 4L, 182L, 65L, 77L, 138L, 128L, 30L, 274L)
  #Z <- Z[, -set.outliers]
  
  if(FALSE){
  myfile <- file.path(rdata.folder, paste("wikipedia-", myhts, ".Rdata", sep = ""))
  save(file = myfile, list = c("Z", "i_all"))
  }
  
}

#stopifnot(length(table(sapply(mylist, length))) == 1)
if(length(table(sapply(mylist, length))) != 1){
  print("some splits are not available")
}


library(igraph)
MAT <- cbind(substring(res, 1, 2), substring(res, 3, 5), substring(res, 6, 8), substring(res, 9, 12))
MAT <- cbind("T", MAT)

tags <- cbind(MAT[, 1], sapply(seq(2, ncol(MAT)), function(j)(
  apply(MAT[, seq(j)], 1, paste, collapse = "")
)))

myedges <- data.frame(do.call(rbind, lapply(seq(ncol(tags) - 1), function(j){
  cbind(tags[, j], tags[, j+1])
})))

itree <- graph.data.frame(myedges)
g <- simplify(itree, remove.loops = F)

savepdf("tree-wikipedia")
plot(g, layout = layout.reingold.tilford(g, root=1, circular=T), 
     vertex.size=1, edge.arrow.size=0, vertex.label= NA, vertex.color = "white")
dev.off()

plot(g, layout = layout_as_tree(g), vertex.size=1, edge.arrow.size=0, vertex.label= NA, vertex.color = "white")
