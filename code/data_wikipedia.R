rm(list = ls())
source("config_paths.R")
library(readr)
library(dplyr)

DT <- read_csv("../data/train_1.csv")
DT <- dplyr::filter(DT, grepl(".wikipedia.org", Page))
name_pages <- DT %>% .$Page

meta_DT <- t(sapply(name_pages, function(name_page){
  first_split <- unlist(strsplit(name_page, "\\.wikipedia\\.org"))
  lang <- tail(unlist(strsplit(first_split[1], "_")), 1)
  
  second_split <- unlist(strsplit(first_split[2], "_"))
  access <- second_split[2]
  agent <- second_split[3]
  c(lang, access, agent)
}))
colnames(meta_DT) <- c("lang", "access", "agent")
meta_DT <- tbl_df(meta_DT)

i <- which(meta_DT[, "access"] == "all-access")
meta_DT[i, "access"] <- "AAC"
i <- which(meta_DT[, "access"] == "desktop")
meta_DT[i, "access"] <- "DES"
i <- which(meta_DT[, "access"] == "mobile-web")
meta_DT[i, "access"] <- "MOB"

i <- which(meta_DT[, "agent"] == "all-agents")
meta_DT[i, "agent"] <- "AAG"
i <- which(meta_DT[, "agent"] == "spider")
meta_DT[i, "agent"] <- "SPD"



# No missing values
nb_nas  <- apply(is.na(DT[, -1]) , 1, sum)
id_nonas <- which(nb_nas == 0)

# Few zeroes
nb_zeroes  <- apply(DT[id_nonas, -1] == 0 , 1, sum)
id_fewzeroes <- which(nb_zeroes < 10)

id_selected <- id_nonas[id_fewzeroes]

DT <- DT[id_selected, ]
meta_DT <- meta_DT[id_selected, ]

myfile <- file.path(rdata.folder, paste("wikipediaDT", ".Rdata", sep = ""))
save(file = myfile, list = c("DT", "meta_DT"))

# SUBSET
for(myhts in seq(1, 70)){
  set.seed(myhts)
  
  i_en <- sample(which(meta_DT$lang == "en"), 50)
  i_fr <- sample(which(meta_DT$lang == "fr"), 50)
  i_de <- sample(which(meta_DT$lang == "de"), 50)
  i_ja <- sample(which(meta_DT$lang == "ja"), 50)
  i_ru <- sample(which(meta_DT$lang == "ru"), 50)
  i_zh <- sample(which(meta_DT$lang == "zh"), 50)
  i_all <- c(i_en, i_fr, i_de, i_ja, i_ru, i_zh)
  mymetaDT <- meta_DT[i_all, ]
  
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
  
  myfile <- file.path(rdata.folder, paste("wikipedia", myhts, ".Rdata", sep = ""))
  save(file = myfile, list = c("Z", "i_all"))
}

stop("done")



imat <- NULL
for(myhts in seq(5)){
  myfile <- file.path(rdata.folder, paste("wikipedia", myhts, ".Rdata", sep = ""))
  load(myfile)
  imat <- cbind(imat, i_all)
}
