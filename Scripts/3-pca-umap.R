source("./Scripts/common.R")

if (!require("Seurat")) { 
  install.packages("Seurat")
  library(Seurat)
}

seurat <- readRDS(file=paste(path.data.out, '2-cleaned.rds', sep=''))

# linear transformation 
# ScaleData():
# - Shifts the expression of each gene, so that the mean expression across cells is 0
# - Scales the expression of each gene, so that the variance across cells is 1
#   - This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

# Performing linear dimensional reduction
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

DimPlot(seurat, reduction = "pca")

DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:20)

JackStrawPlot(seurat, dims = 1:15)

saveRDS(seurat, file = paste(path.data.out, "3-scaled-pca.rds", sep=""))
