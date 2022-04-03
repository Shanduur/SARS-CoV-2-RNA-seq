source("./Scripts/common.R")

device <- "pdf"
# device <- NULL

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

seurat <- readRDS(file = paste0(checkpoint_folder, "3-normalized.rds"))

# loginfo("running GMM fitting")
# fit <- Mclust(seurat,
#               G=4, # TODO: how to select this | expectation maximization
#               model="V")
# summary(fit)

# selekcja cech ze względu na wariancje
# model składowych gausowskich
# pierwszy komponent
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
loginfo("scaling data")
seurat <- ScaleData(object = seurat,
                    features = rownames(x = seurat),
                    verbose = FALSE);

# plot variable features with and without labels to identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

loginfo(paste("plot1+plot2 - VariableFeaturePlot"))
print_img(plot1 + plot2,
          title = "VariableFeatures",
          device = device)

loginfo(paste("plot2 - VariableFeaturePlot labeled"))
print_img(plot2,
          title = "VariableFeatures-labeled",
          device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "4-scaled.rds"))
