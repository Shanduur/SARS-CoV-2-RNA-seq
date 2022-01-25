source("./Scripts/common.R")

if (!require("Seurat")) { 
  install.packages("Seurat")
  library(Seurat)
}
if (!require("dplyr")) { 
  install.packages("dplyr")
  library(dplyr)
}

pbmc <- readRDS(file=paste(path.data.out, '3-scaled-pca.rds', sep=''))

# Clustering the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

# The goal of the following algorithms is to learn the underlying manifold of the data 
# in order to place similar cells together in low-dimensional space.
DimPlot(pbmc, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

head(pbmc.markers, n=5)

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VlnPlot(pbmc, features = c("SERPINA1", "CTSH"))

VlnPlot(pbmc, features = c("SPARCL1", "CLDN5"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c(
  "SFTPC", "SLPI", "C1QA", "C1QB", "IL7R", "KLRB1", "MSR1", "APOC1", "EMP2", "AGER"
))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

