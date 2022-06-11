source("./Scripts/common.R")

prefix <- "03"
device <- "pdf"
# device <- NULL

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("resample")) {
  install.packages("resample")
  library(resample)
}
if (!require("matrixStats")) {
  install.packages("matrixStats")
  library(matrixStats)
}

seurat <- readRDS(file = paste0(checkpoint_folder, "2-loaded.rds"))

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 1e4)

loginfo("calculating variances")
variances <- rowVars(as.matrix(GetAssayData(object = seurat, slot = "data")),
                         useNames = "TRUE")

my_hist <- function(x,
                    breaks = 100,
                    xlab = "Variance",
                    ylab = "Number of cells",
                    main = "Histogram of variances") {
  hist(x,
       breaks = breaks,
       main = main,
       xlab = xlab,
       ylab = ylab)
}

loginfo("printing histogram of variances")
print_img(variances,
          fun = my_hist,
          prefix = prefix,
          title = "variances",
          device = device)

variances <- as.data.frame(variances)
variances <- data.frame("H"=rownames(variances), variances)

# wektor wariancji eksportujemy (wyliczony ze znormalizowanych danych)
loginfo("exporting table of row variances to file")
write.table(variances,
            paste0(export_folder, "variances.data.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

saveRDS(seurat, file = paste0(checkpoint_folder, "3-normalized.rds"))
