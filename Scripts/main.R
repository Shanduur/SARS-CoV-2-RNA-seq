source("./Scripts/common.R")

prefix <- "02"
device <- "pdf"

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("mclust")) {
  install.packages("mclust")
  library(mclust)
}
if (!require("resample")) {
  install.packages("resample")
  library(resample)
}
if (!require("matrixStats")) {
  install.packages("matrixStats")
  library(matrixStats)
}
if (!require("parallelDist")) {
  install.packages("parallelDist")
  library(parallelDist)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require("cowplot")) {
  install.packages("cowplot")
  library("cowplot")
}

samples <- c(
  # "fibrosis-01",
  # "fibrosis-02",
  "cap-ctrl",
  "covid-ctrl",
  "covid-cap"
  # "smokers",
  # "non-smokers"
  )

files <- c(
  "./Data/Pneumonia/GSE164948_covid_control_RNA_counts.csv"
)

for (i in 1:length(files)) {
  if (!file.exists(files[i])) {
    logerror(paste("data file", files[i], "does not exist!"))
  } else {
    loginfo(paste("data file location", files[i], "ok"))
  }
}

meta <- c(
  "./Data/Pneumonia/GSE164948_covid_control_count_metadata_added.csv"
)

for (i in 1:length(meta)) {
  if (!is.null(meta[i]) && !file.exists(meta[i])) {
    logerror(paste("meta file", meta[i], "does not exist!"))
  } else if (!is.null(meta[i]) && file.exists(meta[i])) {
    loginfo(paste("meta file location", files[i], "ok"))
  } else if (is.null(meta[i])) {
    logwarn("meta file location not privided")
  }
}

seurat_list <- list()
for (i in 1:length(files)) {
  loginfo(paste("loading file", paste0(i, ":"), files[i], "|", meta[i]))
  seurat_list[[i]] <- load_seurat(
    filename = files[i],
    project = samples[i],
    meta = meta[i],
    meta_separator = "\t"
  )
  seurat_list[[i]][["DataSet"]] <- samples[i]
}

if (length(files) >= 5) {
  loginfo("saving full seurat list")
  saveRDS(seurat_list, file = paste0(checkpoint_folder, "2-seurat-list-full.rds"))
}

loginfo(paste("building seurat object"))
seurat <- merge(
  x = seurat_list[[1]],
  y = unlist(seurat_list)[2:length(seurat_list)],
  add.cell.ids = samples,
  project = "SARS-CoV-2-RNA-seq"
)
str(seurat@meta.data)

seurat@meta.data$original <- seurat@meta.data$orig.ident
seurat@meta.data$orig.ident <- NULL

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RB-")
head(seurat[["percent.mt"]], n = 20)

loginfo(paste("histogram 1: nFeature_RNA"))
print_img(table(seurat[["nFeature_RNA"]]),
          fun = hist,
          prefix = prefix,
          title = "nFeature_RNA",
          device = device)

loginfo(paste("histogram 2: percent.mt"))
print_img(table(seurat[["percent.mt"]]),
          fun = hist,
          prefix = prefix,
          title = "percent-mt",
          device = device)

loginfo(paste("histogram 3: percent.rb"))
print_img(table(seurat[["percent.rb"]]),
          fun = hist,
          prefix = prefix,
          title = "percent-rb",
          device = device)

loginfo(paste("histogram 4: nCount_RNA"))
print_img(table(seurat[["nCount_RNA"]]),
          fun = hist,
          prefix = prefix,
          title = "nCount_RNA",
          device = device)

# Visualize QC metrics as a violin plot
loginfo(paste("violin plot 1"))
vln1 <- VlnPlot(seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 2,
  pt.size = 0.001
)
print_img(vln1,
          prefix = prefix,
          title = "violin-qc-metrics-1",
          device = device)

loginfo(paste("violin plot 2"))
vln2 <- VlnPlot(seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 2,
  pt.size = 0
)
print_img(vln2,
          prefix = prefix,
          title = "violin-qc-metrics-2",
          device = device)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") # , group.by = "DataSet")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.rb") # , group.by = "DataSet")
plot3 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # , group.by = "DataSet")

feature_min <- min(seurat@meta.data$nFeature_RNA)
feature_m <- median(seurat@meta.data$nFeature_RNA)
feature_max <- max(seurat@meta.data$nFeature_RNA)
feature_s <- sd(seurat@meta.data$nFeature_RNA)
count_min <- min(seurat@meta.data$nCount_RNA)
count_max <- max(seurat@meta.data$nCount_RNA)
count_m <- mean(seurat@meta.data$nCount_RNA)
count_s <- sd(seurat@meta.data$nCount_RNA)
count_q <- quantile(seurat@meta.data$nCount_RNA, 0.95)

loginfo(paste("plot1 - percent.mt"))
print_img(plot1 + geom_hline(yintercept = 5, color = "red"),
          prefix = prefix,
          title = "plot1-percent.mt",
          device = device)

loginfo(paste("plot2 - percent.rb"))
print_img(plot2 + geom_hline(yintercept = 5, color = "red"),
          prefix = prefix,
          title = "plot2-percent.rb",
          device = device)

loginfo(paste("plot3 - nFeature_RNA is the number of genes detected in each cell, nCount_RNA is the total number of molecules detected within a cell"))
print_img(plot3 +
            geom_hline(yintercept = 500, color = "red") +
            geom_vline(xintercept = count_q, color = "red"),
          prefix = prefix,
          title = "plot3-nFeature_RNA",
          device = device)

loginfo(paste("Feature stats:", feature_min, feature_m, feature_max, feature_s))
loginfo(paste("UMI stats:", count_min, count_m, count_max, count_s, count_q))
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nCount_RNA < count_q & percent.mt < 5 & smoking != -1)

saveRDS(seurat, file = paste0(checkpoint_folder, "2-loaded.rds"))
# rm(seurat_list)

prefix <- "03"

loginfo("performing log normalization")
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 1e4)

loginfo("calculating variances")
variances <- rowVars(as.matrix(GetAssayData(object = seurat, slot = "data")),
                         useNames = "TRUE")
variances <- as.table(variances)

loginfo("printing histogram of variances")
print_img(variances,
          fun = hist,
          prefix = prefix,
          title = "variances",
          device = device)

# wektor wariancji eksportujemy (wyliczony ze znormalizowanych danych)
loginfo("exporting table of row variances to file")
write.table(variances,
            paste0(export_folder, "variances.data.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

saveRDS(seurat, file = paste0(checkpoint_folder, "3-normalized.rds"))

prefix <- "04"

threshold <- 0.181

selected <- variances > threshold
loginfo(paste("number of selected features:", sum(selected)))

# seurat <- subset(x=seurat, features=rownames(selected)[selected])
VariableFeatures(seurat) <- rownames(selected)[selected]

# selekcja cech ze względu na wariancje
# model składowych gausowskich
# pierwszy komponent
# calculating a subset of features that exhibit high cell-to-cell variation in the dataset
# próg powinien zostać dobrany w inny sposób, np. modelowanie (mieszanie rozkładów normalnych)

loginfo("scaling data")
seurat <- ScaleData(object = seurat,
                    features = rownames(x = seurat));

saveRDS(seurat, file = paste0(checkpoint_folder, "4-scaled.rds"))

prefix <- "05"

# Performing linear dimensional reduction
loginfo(paste("run PCA and ProjectDim"))
seurat <- RunPCA(seurat, npcs = 40)
seurat <- ProjectDim(object = seurat)

loginfo(paste("visialize Dim reduction - pca"))
vdl1 <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
print_img(vdl1,
          prefix = prefix,
          title = "viz-dim-PCA",
          device = device)

loginfo(paste("dim plot - pca"))
dm1 <- DimPlot(seurat, reduction = "pca")
print_img(dm1 + NoLegend(),
          prefix = prefix,
          title = "dim-plot-PCA",
          device = device)

loginfo(paste("dim heatmap 1"))
hm1 <- DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
print_img(hm1,
          prefix = prefix,
          title = "dim-heatmap-1",
          device = device)

loginfo(paste("dim heatmap 2"))
hm2 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
print_img(hm2,
          prefix = prefix,
          title = "dim-heatmap-1-15",
          device = device,
          width = 30,
          height = 20)

loginfo("elbow plot")
elb1 <- ElbowPlot(seurat)
print_img(elb1,
          prefix = prefix,
          title = "elbow-plot",
          device = device)

dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
}

distances <- c()
for (i in 1:length(elb1$data$dims)) {
  distances[i] <- dist2d(c(elb1$data$dims[i], elb1$data$stdev[i]),
                                c(elb1$data$dims[1], elb1$data$stdev[1]),
                                c(elb1$data$dims[length(elb1$data$dims)], 
                                  elb1$data$stdev[length(elb1$data$dims)])
                                )
}

inflection_point <- match(max(distances), distances)

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
               group.by = NULL,
               pt.size = 0.1)
print_img(dp1 + NoLegend(),
          prefix = prefix,
          title = "dim-plot-UMAP",
          device = device)

# duplicates are cells with the same PCA scores for the specified dimensions
seurat <- RunTSNE(object = seurat,
                  reduction = "pca",
                  dims = 1:inflection_point,
                  check_duplicates = FALSE);
loginfo(paste("dim plot - tsne"))
dp1 <- DimPlot(object = seurat,
               reduction = "tsne",
               group.by = NULL,
               pt.size = 0.1)
print_img(dp1 + NoLegend(),
          prefix = prefix,
          title = "dim-plot-TSNE",
          device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "5-dim-reduction.rds"))

prefix <- "06"

# The goal of the following algorithms is to learn the underlying manifold of the data
# in order to place similar cells together in low-dimensional space.
loginfo(paste("dim plot - umap"))
dp1 <- DimPlot(object = seurat,
              reduction = "umap",
              # split.by = "DataSet",
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
              # split.by = "DataSet",
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
              # split.by = "DataSet",
              pt.size = 0.1)
print_img(dp3,
          prefix = prefix,
          title = "dim-PCA",
          device = device,
          # width = 25
          )

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
        features = markers$gene[1:2],
        # split.by = "DataSet"
        )
print_img(vln1,
          prefix = prefix,
          title = "violin-1",
          device = device,
          # width = 25
          )

loginfo(paste("violin plot - 2"))
vln2 <- VlnPlot(seurat,
        features = markers$gene[1:2],
        slot = "counts",
        log = TRUE,
        # split.by = "DataSet"
        )
print_img(vln2,
          prefix = prefix,
          title = "violin-2",
          device = device,
          # width = 25
          )

loginfo(paste("feature plots - 1"))
for (i in markers$gene) {
  print(i)
  px <- FeaturePlot(seurat,
              features = i,
              # split.by = "DataSet"
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
hm1 <- DoHeatmap(seurat, features = top10$gene) + NoLegend()
print_img(hm1,
          prefix = prefix,
          title = "heatmap",
          device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "6-biomarkers.rds"))
