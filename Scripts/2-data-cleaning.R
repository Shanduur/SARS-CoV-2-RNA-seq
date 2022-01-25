source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("hdf5r")) {
  install.packages("hdf5r")
  library(hdf5r)
}

path.data = path.data.raw

h5_files <- list.files(path.data, pattern = "*.h5")
h5_files <- mapply(paste, path.data, h5_files, sep='')

# file = 'all'
file = h5_files[1]

if (file == 'all') {
  seurat.data <- lapply(h5_files, Read10X_h5)
  seurat <- lapply(seurat.data, CreateSeuratObject)
  rm(seurat.data)
  seurat  <- merge(seurat[[1]], y = seurat[2:length(seurat)])
} else {
  seurat.data <- Read10X_h5(file)
  seurat <- CreateSeuratObject(seurat.data)
  rm(seurat.data)
}

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.mt"]]

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# TODO: histogram

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# calculating a subset of features that exhibit high cell-to-cell variation in the dataset
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
# TODO: histogram

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and without labels to identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

saveRDS(seurat, file=paste(path.data.out, '2-cleaned.rds', sep=''))
