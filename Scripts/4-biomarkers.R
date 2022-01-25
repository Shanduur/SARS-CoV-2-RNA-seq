source("./Scripts/common.R")

if (!require("dplyr")) { 
  install.packages("dplyr")
}
library(dplyr)

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
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, 
                                ident.1 = 2, 
                                min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, 
                                ident.1 = 5, 
                                ident.2 = c(0, 3), 
                                min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, 
                                ident.1 = 0, 
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VlnPlot(pbmc, features = c("SERPINA1", "CTSH"))

VlnPlot(pbmc, features = c("SPARCL1", "CLDN5"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c(
  "SERPINA1", "CTSH", "APOC1", "NPC2"
))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

