source("./Scripts/common.R")

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
seurat <- readRDS(file = paste0(checkpoint_folder, "2-loaded.rds"))

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 1e4)

# loginfo("exporting pre-normalization data to file")
# write.table(as.matrix(GetAssayData(object = seurat, slot = "data")),
#             paste0(export_folder, "pre-normalization.data.txt"),
#             sep = '\t',
#             row.names = TRUE,
#             col.names = TRUE,
#             quote = FALSE)

loginfo("calculating variances")
column_variances <- colVars(as.matrix(GetAssayData(object = seurat, slot = "data")))

loginfo("printing histogram of variances")
print_img(column_variances,
          fun = hist,
          title = "02-column-variances",
          device = device)

# wektor wariancji eksportujemy (wyliczony ze znormalizowanych danych)
loginfo("exporting table of column variances to file")
write.table(column_variances,
            paste0(export_folder, "variances.data.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

saveRDS(seurat, file = paste0(checkpoint_folder, "3-normalized.rds"))
