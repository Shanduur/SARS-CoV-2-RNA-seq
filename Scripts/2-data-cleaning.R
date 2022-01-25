source("./Scripts/common.R")

if (!require("Seurat")) { install.packages("Seurat") }
if (!require("hdf5r")) { install.packages("hdf5r") }

path.data = path.data.raw

h5_files <- list.files(path.data, pattern = "*.h5")
h5_files <- mapply(paste, path.data, h5_files, sep='')

# file = 'all'
file = h5_files[1]

library(hdf5r)
library(Seurat)

if (file == 'all') {
  pbmc.data <- lapply(h5_files, Read10X_h5)
  pbmc <- lapply(pbmc.data, CreateSeuratObject)
  rm(pbmc.data)
  pbmc  <- merge(pbmc[[1]], y = pbmc[2:length(pbmc)])
} else {
  pbmc.data <- Read10X_h5(file)
  pbmc <- CreateSeuratObject(pbmc.data)
  rm(pbmc.data)
}

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.mt"]]

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# calculating a subset of features that exhibit high cell-to-cell variation in the dataset
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels to identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

saveRDS(pbmc, file=paste(path.data.out, '2-cleaned.rds', sep=''))
