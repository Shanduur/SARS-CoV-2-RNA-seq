source("./Scripts/common.R")

prefix <- "06"
device <- "pdf"
# device <- NULL

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
}

BiocManager::install("SingleR")

seurat <- readRDS(file = paste0(checkpoint_folder, "6-biomarkers.rds"))

saveRDS(seurat, file = paste0(checkpoint_folder, "7-annotate.rds"))
