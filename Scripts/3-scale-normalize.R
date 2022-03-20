source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

seurat <- readRDS(file = paste0(output_folder, "2-loaded.rds"))

# Normalizing the data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 1e4)

# calculating a subset of features that exhibit high cell-to-cell variation in the dataset
# próg powinien zostać dobrany w inny sposób, np. modelowanie (mieszanie rozkładów normalnych)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
loginfo(paste("Number of Variable Features:", length(x = VariableFeatures(object = seurat))));

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)
loginfo(paste("Top 10 most highly variable genes:", toString(top10)))

# linear transformation
# ScaleData():
# - Shifts the expression of each gene, so that the mean expression across cells is 0
# - Scales the expression of each gene, so that the variance across cells is 1
#   - This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
seurat <- ScaleData(object = seurat,
                    features = rownames(x = seurat),
                    verbose = FALSE);

# plot variable features with and without labels to identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

loginfo(paste("plot1+plot2 - VariableFeaturePlot"))
print_img(plot1 + plot2)

loginfo(paste("plot2 - VariableFeaturePlot labeled"))
print_img(plot2)

saveRDS(seurat, file = paste0(output_folder, "3-scaled-normalized.rds"))
