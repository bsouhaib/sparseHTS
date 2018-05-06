main.folder   <- "/Users/souhaibt/Dropbox/Code/hts_sparse"

code.folder <- file.path(main.folder, "code")
work.folder <- file.path(main.folder, "work")


rdata.folder  <- file.path(work.folder , "rdata")
bf.folder  <- file.path(rdata.folder , "bf")
results.folder  <- file.path(rdata.folder , "results")
pdf.folder <- file.path(work.folder , "pdfs")

x <- c(work.folder, rdata.folder, bf.folder, results.folder, pdf.folder)
sapply(x, dir.create, showWarnings = FALSE, recursive = TRUE)