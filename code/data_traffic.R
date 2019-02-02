rm(list = ls())
source("config_paths.R")
library(readr)

#DT <- read_delim("LD2011_2014.txt", delim = ";")
#ft <- factor(rep(seq(974), each = 144))
#res <- tapply(DT[, seq(2)], ft, sum)
#https://robjhyndman.com/hyndsight/seasonal-periods/

myfct <- function(filepath){
  con = file(filepath, "r")
  x = readLines(con, n = 1)
  close(con)
  x <- gsub("\\[", "", x)
  x <- gsub("\\]", "", x)
  x <- unlist(sapply(x, strsplit, " "))
  x <- as.numeric(x)
}


trainlabels <- myfct(file.path(main.folder, "data", "PEMS-SF/PEMS_trainlabels"))
testlabels <- myfct(file.path(main.folder, "data", "PEMS-SF/PEMS_testlabels"))
labels <- c(trainlabels, testlabels)

randperm <- myfct(file.path(main.folder, "data", "PEMS-SF/randperm"))

i <- 0
list_dt <- vector("list", 440)
for(filepath in c("PEMS-SF/PEMS_train", "PEMS-SF/PEMS_test")){
  print(filepath)
  con = file(file.path(main.folder, "data", filepath), "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    i <- i + 1
    print(i)
    x <- line
    res <- unlist(strsplit(x, ";"))
      
    res2 <- gsub("\\[", "", res)
    res3 <- gsub("\\]", "", res2)
    m <- lapply(res3, strsplit, " ")
    m <- lapply(m, unlist)
    m <- lapply(m, as.numeric)
    vars <- sapply(m, sum)
    
    list_dt[[i]] <- m
      
    # 263 lines
    # m has 963 component (variables/sensors) and each component has 144 entries (one day)
    #sapply(m, length)
    #  263 + 173
    
  }
  close(con)
}

correct_order <- sort(randperm, index = T)$ix
dw <- labels[correct_order]
final_list <- list_dt[correct_order]

# length(final_list[[1]])
# 963 = 6 x 24 heures

res <- sapply(final_list, function(mylist){
  sapply(mylist, sum)
})

myres <- matrix(NA, nrow = nrow(res), ncol = 65*7)
mylist <- list(c(1, 2, 3),NA,NA, 2,NA,NA , NA ,2,NA,  NA,1,NA,NA , 
               NA ,NA,NA , NA  ,NA,NA , NA  ,NA ,2 ,NA,NA , NA  ,NA ,
               6,NA, NA,NA , NA  ,NA,NA , NA  ,NA,2,NA,  NA , NA,NA , 
               NA  ,NA,2,NA, NA , NA,2,NA,  NA,1,NA, NA,NA , NA  ,NA,NA,
               NA  ,NA,NA , NA,2,NA,NA , NA,c(6, 7))


mystatus <- unlist(lapply(mylist, function(missing_in_day){
  if(length(missing_in_day) == 1 && is.na(missing_in_day)){
    x <- rep("OK", 7)
  }else{
    x <- rep("OK", 7)
    x[missing_in_day] <- "KO"
  }
  x
}))
myres[seq(nrow(res)), which(mystatus == "OK")] <- res

m <- seq(dmy("01/01/2008"), dmy("30/03/2009"), by = "day")
id <- which(year(m) == 2008)
DT <- myres[, id]
save(file = file.path(rdata.folder, "road_trafficDT.Rdata"), list = c("DT"))

stop("done")


y <- ts(myres[40, ], freq = 7)
y <- tsclean(y)
plot(y)

obj_stl <- stl(y, s.window = "periodic")
plot(obj_stl)

ydetrend <- y - obj_stl$time.series[, "trend"]
plot(ydetrend)
lines(obj_stl$time.series[, "seasonal"], col = "red")

plot.ts(obj_stl$time.series[, "remainder"])

x <- obj_stl$time.series[, "remainder"]
plot(x)
model <- ets(x)
fit <- fitted(model)
lines(fit, col = "red")

model2 <- auto.arima(x)
fit2 <- fitted(model2)
lines(fit2, col = "blue")









