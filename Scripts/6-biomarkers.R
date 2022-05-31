source("./Scripts/common.R")

prefix <- "06"
device <- "pdf"
# device <- NULL

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

seurat <- readRDS(file = paste0(checkpoint_folder, "5-dim-reduction.rds"))

# The goal of the following algorithms is to learn the underlying manifold of the data
# in order to place similar cells together in low-dimensional space.
loginfo(paste("dim plot - umap"))
dp1 <- DimPlot(object = seurat,
              reduction = "umap",
              split.by = "smoking",
              pt.size = 0.1)
print_img(dp1,
          prefix = prefix,
          title = "dim-UMAP",
          device = device,
          # width = 25
          )

loginfo(paste("dim plot - tsne"))
dp2 <- DimPlot(object = seurat,
              reduction = "tsne",
              split.by = "smoking",
              pt.size = 0.1)
print_img(dp2,
          prefix = prefix,
          title = "dim-TSNE",
          device = device,
          # width = 25
          )

loginfo(paste("dim plot - pca"))
dp3 <- DimPlot(object = seurat,
              reduction = "pca",
              split.by = "smoking",
              pt.size = 0.1)
print_img(dp3,
          prefix = prefix,
          title = "dim-PCA",
          device = device,
          # width = 25
          )

# loginfo(paste("dim plot - ica"))
# dp4 <- DimPlot(object = seurat,
#               reduction = "ica",
#               split.by = "smoking",
#               pt.size = 0.1)
# print_img(dp4,
#           prefix = prefix,
#           title = "dim-ICA",
#           device = device,
#           # width = 25
#           )

# loginfo(paste("dim plot - mds"))
# dp5 <- DimPlot(object = seurat,
#               reduction = "mds",
#               split.by = "DataSet",
#               pt.size = 0.1)
# print_img(dp5,
#           title = "dim-MDS",
#           device = device)

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

loginfo(paste("violin plot - 1"))
vln1 <- VlnPlot(seurat,
        features = unique(markers$gene),
        split.by = "smoking"
        )
print_img(vln1,
          prefix = prefix,
          title = "violin-1",
          device = device,
          # width = 25
          )

loginfo(paste("violin plot - 2"))
vln2 <- VlnPlot(seurat,
        features = unique(markers$gene),
        slot = "counts",
        log = TRUE,
        split.by = "smoking"
        )
print_img(vln2,
          prefix = prefix,
          title = "violin-2",
          device = device,
          width = 25,
          height = 18
          )

loginfo(paste("feature plots - 1"))
for (i in markers$gene) {
  print(i)
  px <- FeaturePlot(seurat,
              features = i,
              split.by = "smoking"
              )
  print_img(px,
            prefix = prefix,
            title = paste0("features-plot-", i),
            device = device,
            # width = 25
            )
}

seurat_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
hm1 <- DoHeatmap(seurat, features = unique(markers$gene)) + NoLegend()
print_img(hm1,
          prefix = prefix,
          title = "heatmap",
          device = device)

loginfo("exporting table of row variances to file")
write.table(markers,
            paste0(export_folder, "markers.data.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

saveRDS(seurat, file = paste0(checkpoint_folder, "6-biomarkers.rds"))
