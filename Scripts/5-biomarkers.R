source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require("cowplot")) {
  install.packages("cowplot")
  library("cowplot")
}

seurat <- readRDS(file = paste0(output_folder, "4-pca-jackstraw.rds"))

# Clustering the cells
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat,
                       resolution = 0.5) # use a value above (below) 1.0 if
                                         # you want to obtain a larger (smaller)
                                         # number of communities

head(Idents(seurat), 5)

seurat <- RunUMAP(object = seurat, dims = 1:10)
# duplicates are cells with the same PCA scores for the specified dimensions
seurat <- RunTSNE(object = seurat,
                  reduction = "pca",
                  dims = 1:10,
                  check_duplicates = FALSE);


# The goal of the following algorithms is to learn the underlying manifold of the data
# in order to place similar cells together in low-dimensional space.
p1 <- DimPlot(object = seurat,
              reduction = "umap",
              split.by = "DataSet",
              pt.size=0.1)
p2 <- DimPlot(object = seurat,
              reduction = "tsne",
              split.by = "DataSet",
              pt.size=0.1)
p3 <- DimPlot(object = seurat,
              reduction = "pca",
              split.by = "DataSet",
              pt.size=0.1)
print(p1)
print(p2)
print(p3)

plot_grid(p1, p2)

# Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_markers <- FindAllMarkers(seurat,
  only.pos = TRUE,
  min.pct = 0.1, # only test genes that are detected in a minimum 
                 # fraction of min.pct cells in either of the two populations
  logfc.threshold = 0.25 # default value
)
markers <- (seurat_markers %>%
              group_by(cluster) %>%
              slice_max(n = 2, order_by = avg_log2FC))

head(seurat_markers, n = 5)

print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

VlnPlot(seurat,
        features = markers$gene[1:2],
        split.by = "DataSet")

VlnPlot(seurat,
        features = markers$gene[1:2],
        slot = "counts",
        log = TRUE,
        split.by = "DataSet")

for (i in markers$gene) {
  print(i)
  px <- FeaturePlot(seurat, 
              features = i,
              split.by = "DataSet")
  print(px)
}

seurat_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()
