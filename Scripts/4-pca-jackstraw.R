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

# selekcja cech ze względu na wariancje
# model składowych gausowskich
# pierwszy komponent
# Gamred
# SNN
# mclust

# Performing linear dimensional reduction
loginfo(paste("run PCA and ProjectDim"))
seurat <- RunPCA(seurat, npcs = 40)
seurat <- ProjectDim(object = seurat)

loginfo(paste("run ICA and ProjectDim"))
seurat <- RunICA(seurat, nics = 40)
seurat <- ProjectDim(object = seurat)

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
print(seurat[["ica"]], dims = 1:5, nfeatures = 5)

loginfo(paste("visialize Dim reduction - pca"))
vdl1 <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
print_img(vdl1)

loginfo(paste("visialize Dim reduction - ica"))
vdl2 <- VizDimLoadings(seurat, dims = 1:2, reduction = "ica")
print_img(vdl2)

loginfo(paste("visialize Dim reduction - mds"))
vdl3 <- VizDimLoadings(seurat, dims = 1:2, reduction = "mds")
print_img(vdl3)

loginfo(paste("dim plot - pca"))
dm1 <- DimPlot(seurat, reduction = "pca")
print_img(dm1)

loginfo(paste("dim plot - ica"))
dm2 <- DimPlot(seurat, reduction = "ica")
print_img(dm2)

loginfo(paste("dim plot - mds"))
dm3 <- DimPlot(seurat, reduction = "mds")
print_img(dm3)

loginfo(paste("dim heatmap 1"))
hm1 <- DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
print_img(hm1)

loginfo(paste("dim heatmap 2"))
hm2 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
print_img(hm2)

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

loginfo(paste("jackstraw plot"))
jsp1 <- JackStrawPlot(seurat, dims = 1:15)
print_img(jsp1)

saveRDS(seurat, file = paste0(output_folder, "4-pca-jackstraw.rds"))
