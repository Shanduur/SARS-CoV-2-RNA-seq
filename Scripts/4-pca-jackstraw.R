source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("parallelDist")) {
  install.packages("parallelDist")
  library(parallelDist)
}

seurat <- readRDS(file = paste0(output_folder, "3-scaled-normalized.rds"))

# Performing linear dimensional reduction
print(paste("run PCA and ProjectDim"))
seurat <- RunPCA(seurat, npcs = 40)
seurat <- ProjectDim(object = seurat)

print(paste("run ICA and ProjectDim"))
seurat <- RunICA(seurat, nics = 40)
seurat <- ProjectDim(object = seurat)


print(paste("calculate distances between data"))
d <- parDist(t(GetAssayData(seurat, slot = "scale.data")))

print(paste("run CMD scaling"))
mds <- cmdscale(d = d, k = 2)

print(paste("create subobject with reduced dimensions"))
colnames(mds) <- paste0("MDS_", 1:2)
seurat[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(seurat))

# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
print(seurat[["ica"]], dims = 1:5, nfeatures = 5)

print(paste("visialize Dim reduction - pca"))
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

print(paste("visialize Dim reduction - ica"))
VizDimLoadings(seurat, dims = 1:2, reduction = "ica")

print(paste("dim plot - pca"))
DimPlot(seurat, reduction = "pca")

print(paste("dim plot - ica"))
DimPlot(seurat, reduction = "ica")

print(paste("dim plot - mds"))
DimPlot(seurat, reduction = "mds")

print(paste("dim heatmap 1"))
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)

print(paste("dim heatmap 2"))
DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:20)

# smn - algorytmy grupowania grafowe, nieliniowe

pc_pval <- seurat@reductions$pca@jackstraw@overall.p.values; # get p-value for each PC
write.table(x = pc_pval,
            file = paste0(output_folder, "PCA_jackstraw_scores.xls"),
            quote = FALSE,
            sep = "\t",
            col.names = TRUE);

print(paste("jackstraw plot"))
JackStrawPlot(seurat, dims = 1:15)

saveRDS(seurat, file = paste0(output_folder, "4-pca-jackstraw.rds"))
