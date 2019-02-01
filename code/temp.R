myhts <- 1
myfile <- file.path(rdata.folder, paste("wikipedia", myhts, ".Rdata", sep = ""))
load(myfile)
Z <- ts(Z, freq = 7)

x <- Z[, 20]
plot(x)
x <- log(x)
plot(x)
obj_outliers <- tsoutliers(x)
x[obj_outliers$index] <- obj_outliers$replacements
lines(x, col = "red")
plot.ts(x, col = "red")
plot(stl(x, s.window = "period"))
