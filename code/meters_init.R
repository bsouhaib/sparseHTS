rm(list = ls())
library(lubridate)
library(gdata)
library(dplyr)
library(igraph)

meters.work.folder <- file.path("/home/rstudio/PROJ", "work")
sparsehts.work.folder <- file.path("/home/rstudio/sparseHTS", "work")
initmeters.folder  <- file.path(file.path("/home/rstudio/PROJ", "procdata"), "initmeters")

source(file.path("/home/rstudio/PROJ/code", "config_splitting.R"))

makelist <- function(vecintervals){
  sol <- lapply(seq(length(vecintervals)), function(i){list(vecintervals[i])})
  mylist <- lapply(sol, "[[", 1)
}

load(file.path(meters.work.folder, "info.Rdata"))

# "UKG" "UKL" "UKJ" "UKI" "UKM" "UKF" "---" "UKK" "UKD"
myregion <- "UKF"

subInfoDT <- infoDT %>%
  filter(NUTS1 %in% myregion)

allintervals <- subInfoDT %>%
  transmute(interval = lubridate::interval(firstAdvance, lastAdvance)) %>%
  .$interval

listintervals <- makelist(allintervals)

myinterval <- interval(startObs, endObs)
seq_myinterval <- seq(startObs, endObs, by = "30 min")

matches <- lapply(listintervals, function(oneinterval){ lubridate::intersect(oneinterval, myinterval) == myinterval })
metersInInterval <- subInfoDT[which(unlist(matches)), ] %>% .$IDMETER
print(length(metersInInterval))


pctFound <- n <- n_na <- n_expected <- n_avail <- numeric(length(metersInInterval)) + NA
listmissing <- vector("list", length(metersInInterval))
for(i in seq_along(metersInInterval)){
  print(i)
  meter <- metersInInterval[i]
  
  infoMeter <- subInfoDT %>% filter(IDMETER == meter) %>% dplyr::select(firstAdvance, lastAdvance) 
  firstAdvance <- infoMeter %>% .$firstAdvance
  lastAdvance  <- infoMeter %>% .$lastAdvance
  alldates <- seq(firstAdvance, lastAdvance, by = "30 min")
  ids <- match(seq_myinterval, alldates) 
  stopifnot(all(!is.na(ids)))
  
  load(file.path(initmeters.folder, paste("meter-", meter, ".Rdata", sep = "")))
  n[i] <- nrow(dataset)
  n_expected[i] <- length(alldates)
  n_na[i] <- length(which(is.na(dataset[ids, 1])))
  n_avail[i] <- n_expected[i] - n_na[i]
  pctFound[i] <- 1 - (n_na[i]/n[i])
}

id <- which(pctFound > 0.99)
finalmeters <- metersInInterval[id]

# Keep meters with few missing values
x <- subInfoDT %>% filter(IDMETER %in% finalmeters)

# Keep meters with complete NUTS information
x <- x %>% filter(!grepl("-", NUTS4))

# Remove few weird meters
x <- x %>% filter(!IDMETER %in% c(6228, 13154, 9503))

# Some meter with high consumption during the night 
x <- x %>% filter(!IDMETER %in% c(12874L, 6951L, 14738L, 925L, 8255L))

# Each node must have at least two children nodes (NUTS HIERARCHY)
idset <- which(x[, "NUTS4"] == "UKF2100")
res <- split(idset, c(1,2))
x[res[[1]], "NUTS4"] <- "UKF2100"
x[res[[2]], "NUTS4"] <- "UKF2101"

idset <- which(x[, "NUTS4"] == "UKF2202")
res <- split(idset, c(1,2))
x[res[[1]], "NUTS4"] <- "UKF2202"
x[res[[2]], "NUTS4"] <- "UKF2209"

idset <- which(x[, "NUTS4"] == "UKF1100")
res <- split(idset, c(1,2))
x[res[[1]], "NUTS4"] <- "UKF1100"
x[res[[2]], "NUTS4"] <- "UKF1101"

idset <- which(x[, "NUTS4"] == "UKF1400")
res <- split(idset, c(1,2))
x[res[[1]], "NUTS4"] <- "UKF1400"
x[res[[2]], "NUTS4"] <- "UKF1401"

for(mynuts4 in c("UKF3004", "UKF3006")){
  idset <- which(x[, "NUTS4"] == mynuts4)
  x[idset, "NUTS4"] <- paste("UKF31", substr(mynuts4, 6,7), sep = "")
  x[idset, "NUTS3"] <- "UKF31"
}

# Remove some branches in DEMO HIERARCHY
DT <- x %>% filter(DEMO2 != "D517", DEMO1 != "D2")


print("NEW SELECTION !!")
 
hierarchy_names <- c("UKF23", "UKF16", "UKF21", "UKF13", "UKF15")

for(hierarchy_name in hierarchy_names){
  
  x <- DT
  print(x[which(grepl(hierarchy_name, x$NUTS4)), ] %>% count(NUTS4))
  x <- x[which(grepl(hierarchy_name, x$NUTS4)), ]
  print(dim(x))
  
  # Save myinfo.Rdata
  myinfoDT <- x
  
  lname <- nchar(hierarchy_name)
  if(lname == 3){
    myedges <- data.frame(rbind(cbind(myinfoDT$NUTS1, myinfoDT$NUTS2), cbind(myinfoDT$NUTS2, myinfoDT$NUTS3),
                                cbind(myinfoDT$NUTS3, myinfoDT$NUTS4) ))
  }else if(lname == 4){
    myedges <- data.frame(rbind(cbind(myinfoDT$NUTS2, myinfoDT$NUTS3), cbind(myinfoDT$NUTS3, myinfoDT$NUTS4), cbind(myinfoDT$NUTS4, myinfoDT$IDMETER)))
  }else if(lname == 5){
    myedges <- data.frame(rbind(cbind(myinfoDT$NUTS3, myinfoDT$NUTS4), cbind(myinfoDT$NUTS4, myinfoDT$IDMETER)))
  }else{
    stop("ERROR !!")
  }
  
  
  itree <- graph.data.frame(myedges)
  itree <- simplify(itree, remove.loops = F)
  # plot(itree, layout = layout.reingold.tilford(itree, root=1, circular=T), vertex.label.cex = 0.4, vertex.size = 1, vertex.label.dist = .2)
  # MUCH BETTER: plot(itree, layout = layout.reingold.tilford(itree, root=1, circular=T), vertex.size=0, vertex.label=NA, edge.arrow.size=0)
  #browser()
  
  all.nodes.names <- V(itree)$name
  agg.nodes.names <- aggSeries <- all.nodes.names[which(degree(itree, V(itree), "out")!=0)]
  n_agg <- length(agg.nodes.names)
  
  bottom.nodes.names <- bottomSeries <- all.nodes.names[which(degree(itree, V(itree), "out")==0)]
  n_bottom <- length(bottom.nodes.names)
  
  Sagg <- matrix(0, nrow = n_agg, ncol = n_bottom)
  for(i in seq_along(agg.nodes.names)){
    agg.node.name <- agg.nodes.names[i]
    reachable <- which(shortest.paths(itree, agg.node.name, mode="out") != Inf)
    terminal.nodes <- reachable[which(degree(itree, reachable, mode="out") == 0)]
    #print(terminal.nodes)
    terminal.nodes.names <- all.nodes.names[terminal.nodes]
    ids <- match(terminal.nodes.names, bottomSeries)
    stopifnot(all(!is.na(ids)))
    Sagg[i, ids] <- 1
  }
  
  #file.remove(file.path(sparsehts.work.folder, "myinfo.Rdata"))
  save(file = file.path(sparsehts.work.folder, paste("myinfo_", hierarchy_name,".Rdata", sep = "")) , list = c("bottomSeries", "itree", "Sagg", "aggSeries", "n_agg", "n_bottom"))
}

