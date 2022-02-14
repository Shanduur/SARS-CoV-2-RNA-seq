source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

seurat <- readRDS(file = paste(output_folder, "3-scaled-normalized.rds", sep = ""))

# Performing linear dimensional reduction
seurat <- RunPCA(seurat, npcs = 100)

seurat <- ProjectDim(object = seurat)

# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

DimPlot(seurat, reduction = "pca")

DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:20)

pc.pval <- seurat@reductions$pca@jackstraw@overall.p.values; # get p-value for each PC
write.table(x = pc.pval,
            file = paste(output_folder, "PCA_jackstraw_scores.xls"),
            quote = FALSE,
            sep = '\t',
            col.names = TRUE);

JackStrawPlot(seurat, dims = 1:15)

saveRDS(seurat, file = paste(output_folder, "4-pca-jackstraw.rds", sep = ""))
