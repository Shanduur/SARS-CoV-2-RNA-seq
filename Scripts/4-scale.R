source("./Scripts/common.R")

prefix <- "04"
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

variances <- read.table(paste0(export_folder, "variances.data.txt"), header = TRUE)
threshold_log <- 0.181
selected <- variances$variances > threshold_log
loginfo(paste("number of selected features:", sum(selected)))
VariableFeatures(seurat) <- rownames(seurat)[selected]
seurat <- subset(x=seurat, features=rownames(seurat)[selected])

var_mean <- mean(variances$variances)
var_med <- median(variances$variances)

loginfo(paste("mean of variances =", var_mean))
loginfo(paste("median of variances =", var_med))

print_img(ggplot(variances, aes(x=variances)) +
            geom_histogram(bins=100, color="black", fill="grey") +
            theme_bw() +
            scale_x_continuous(trans='log2') +
            geom_vline(xintercept=threshold_log, color = "red") +
            geom_vline(xintercept=var_mean, color = "green") +
            geom_vline(xintercept=var_med, color = "blue") +
            ggtitle("Histogram of annotation scores") +
            xlab("Variance") +
            ylab("Number of cells"),
          prefix = prefix,
          title = "variances-selected",
          device = device)
 
# selekcja cech ze względu na wariancje
# model składowych gausowskich
# pierwszy komponent
# calculating a subset of features that exhibit high cell-to-cell variation in the dataset
# próg powinien zostać dobrany w inny sposób, np. modelowanie (mieszanie rozkładów normalnych)
# seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = length(selected))
# loginfo(paste("Number of Variable Features:", length(x = VariableFeatures(object = seurat))));

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(seurat), 10)
# loginfo(paste("Top 10 most highly variable genes:", toString(top10)))

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
# plot1 <- VariableFeaturePlot(seurat)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# 
# loginfo(paste("plot1+plot2 - VariableFeaturePlot"))
# print_img(plot1 + plot2,
#           prefix = prefix,
#           title = "VariableFeatures",
#           device = device)
# 
# loginfo(paste("plot2 - VariableFeaturePlot labeled"))
# print_img(plot2,
#           prefix = prefix,
#           title = "VariableFeatures-labeled",
#           device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "4-scaled.rds"))
