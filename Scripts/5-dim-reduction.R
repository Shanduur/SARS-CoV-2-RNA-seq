source("./Scripts/common.R")

prefix <- "05"
device <- "pdf"
# device <- NULL

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("parallelDist")) {
  install.packages("parallelDist")
  library(parallelDist)
}

seurat <- readRDS(file = paste0(checkpoint_folder, "4-scaled.rds"))

# Gamred
# SNN
# mclust

# Performing linear dimensional reduction
loginfo(paste("run PCA and ProjectDim"))
seurat <- RunPCA(seurat, npcs = 40)
seurat <- ProjectDim(object = seurat)

# loginfo(paste("run ICA and ProjectDim"))
# seurat <- RunICA(seurat, nics = 40)
# seurat <- ProjectDim(object = seurat)

# loginfo(paste("calculate distances between data"))
# d <- parDist(t(GetAssayData(seurat, slot = "scale.data")))
#
# loginfo(paste("run CMD scaling"))
# mds <- cmdscale(d = d, k = 2)
#
# loginfo(paste("create subobject with reduced dimensions"))
# colnames(mds) <- paste0("MDS_", 1:2)
# seurat[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(seurat))

# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
# print(seurat[["ica"]], dims = 1:5, nfeatures = 5)

loginfo(paste("visialize Dim reduction - pca"))
vdl1 <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
print_img(vdl1,
          prefix = prefix,
          title = "VizDim-PCA",
          device = device)

# loginfo(paste("visialize Dim reduction - ica"))
# vdl2 <- VizDimLoadings(seurat, dims = 1:2, reduction = "ica")
# print_img(vdl2,
#           prefix = prefix,
#           title = "VizDim-ICA",
#           device = device)

# loginfo(paste("visialize Dim reduction - mds"))
# vdl3 <- VizDimLoadings(seurat, dims = 1:2, reduction = "mds")
# print_img(vdl3,
#           title = "VizDim-MDS",
#           device = device)

loginfo(paste("dim plot - pca"))
dm1 <- DimPlot(seurat,
               reduction = "pca",
               label.box = TRUE,
               group.by = "smoking")
print_img(dm1,
          prefix = prefix,
          title = "DimPlot-PCA",
          device = device)

# loginfo(paste("dim plot - ica"))
# dm2 <- DimPlot(seurat, reduction = "ica")
# print_img(dm2,
#           prefix = prefix,
#           title = "DimPlot-ICA",
#           device = device)

# loginfo(paste("dim plot - mds"))
# dm3 <- DimPlot(seurat, reduction = "mds")
# print_img(dm3,
#           title = "DimPlot-MDS",
#           device = device)

loginfo(paste("dim heatmap 1"))
hm1 <- DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
print_img(hm1,
          prefix = prefix,
          title = "DimHeatmap-1",
          device = device)

loginfo(paste("dim heatmap 2"))
hm2 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
print_img(hm2,
          prefix = prefix,
          title = "DimHeatmap-1-15",
          device = device)

dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
}

loginfo("elbow plot")
elb1 <- ElbowPlot(seurat)

distances <- c()
for (i in 1:length(elb1$data$dims)) {
  distances[i] <- dist2d(c(elb1$data$dims[i], elb1$data$stdev[i]),
                                c(elb1$data$dims[1], elb1$data$stdev[1]),
                                c(elb1$data$dims[length(elb1$data$dims)], 
                                  elb1$data$stdev[length(elb1$data$dims)])
                                )
}

inflection_point <- match(max(distances), distances)
loginfo(paste("inflection point is equal to", inflection_point))

print_img(elb1 + geom_vline(xintercept = inflection_point, color = "red"),
          prefix = prefix,
          title = "ElbowPlot",
          device = device)

# # Determine the ‘dimensionality’ of the dataset
# seurat <- JackStraw(seurat, num.replicate = 100)
# seurat <- ScoreJackStraw(seurat, dims = 1:inflection_point)
# 
# pc_pval <- seurat@reductions$pca@jackstraw@overall.p.values; # get p-value for each PC
# write.table(x = pc_pval,
#             file = paste0(output_folder, "PCA_jackstraw_scores.xls"),
#             quote = FALSE,
#             sep = "\t",
#             col.names = TRUE);
# 
# loginfo(paste("jackstraw plot"))
# jsp1 <- JackStrawPlot(seurat, dims = 1:inflection_point)
# print_img(jsp1,
#           prefix = prefix,
#           title = "JackStrawPlot",
#           device = device)

# Clustering the cells
seurat <- FindNeighbors(seurat, dims = 1:inflection_point)
seurat <- FindClusters(seurat,
                       resolution = 0.5) # use a value above (below) 1.0 if
                                         # you want to obtain a larger (smaller)
                                         # number of communities

head(Idents(seurat), 5)

seurat <- RunUMAP(object = seurat, dims = 1:inflection_point)
loginfo(paste("dim plot - umap"))
dp1 <- DimPlot(object = seurat,
               reduction = "umap",
               label.box = TRUE,
               pt.size = 0.1,
               group.by = "smoking")
print_img(dp1,
          prefix = prefix,
          title = "DimPlot-UMAP",
          device = device)

# duplicates are cells with the same PCA scores for the specified dimensions
seurat <- RunTSNE(object = seurat,
                  reduction = "pca",
                  dims = 1:inflection_point,
                  check_duplicates = FALSE);
loginfo(paste("dim plot - tsne"))
dp1 <- DimPlot(object = seurat,
               reduction = "tsne",
               pt.size = 0.1,
               label.box = TRUE,
               group.by = "smoking")
print_img(dp1,
          prefix = prefix,
          title = "DimPlot-TSNE",
          device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "5-dim-reduction.rds"))
