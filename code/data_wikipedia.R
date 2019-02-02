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