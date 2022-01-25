source("./Scripts/common.R")

if (!require("Seurat")) { 
  install.packages("Seurat")
  library(Seurat)
}
if (!require("dplyr")) { 
  install.packages("dplyr")
  library(dplyr)
}

seurat <- readRDS(file=paste(path.data.out, '3-scaled-pca.rds', sep=''))

# Clustering the cells
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)

head(Idents(seurat), 5)

seurat <- RunUMAP(seurat, dims = 1:10)

# The goal of the following algorithms is to learn the underlying manifold of the data 
# in order to place similar cells together in low-dimensional space.
DimPlot(seurat, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.markers <- FindAllMarkers(seurat, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

head(seurat.markers, n=5)

print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

VlnPlot(seurat, features = c("SERPINA1", "CTSH"))

VlnPlot(seurat, features = c("SPARCL1", "CLDN5"), slot = "counts", log = TRUE)

FeaturePlot(seurat, features = c(
  "SFTPC", "SLPI", "C1QA", "C1QB", "IL7R", "KLRB1", "MSR1", "APOC1", "EMP2", "AGER"
))

seurat.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()

