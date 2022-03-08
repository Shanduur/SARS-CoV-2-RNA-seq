source("./Scripts/common.R")

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

samples <- c("cap-ctrl", "covid-ctrl", "covid-cap")
files <- c("./Data/Pneumonia/GSE164948_cap_control_RNA_counts.csv",
          "./Data/Pneumonia/GSE164948_covid_control_RNA_counts.csv",
          "./Data/Pneumonia/GSE164948_covid_cap_RNA_counts.csv")

seurat_list <- list()

for (i in 1:length(files)) {
  seurat_list[[i]] <- load_seurat(filename = files[i],
                                  project = samples[i])
  seurat_list[[i]][["DataSet"]] <- samples[i];
}

seurat <- merge(x = seurat_list[[1]],
                y = c(seurat_list[[2]], seurat_list[[3]]),
                add.cell.ids = samples,
                project = "SARS-CoV-2-RNA-seq");

str(seurat@meta.data)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RB-")
head(seurat[["percent.mt"]], n = 20)
hist(table(seurat[["nFeature_RNA"]]))
hist(table(seurat[["percent.mt"]]))
hist(table(seurat[["nCount_RNA"]]))

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot3 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

feature_min <- min(seurat@meta.data$nFeature_RNA);
feature_m <- median(seurat@meta.data$nFeature_RNA)
feature_max <- max(seurat@meta.data$nFeature_RNA)
feature_s <- sd(seurat@meta.data$nFeature_RNA)
count_min <- min(seurat@meta.data$nCount_RNA)
count_max <- max(seurat@meta.data$nCount_RNA)
count_m <- mean(seurat@meta.data$nCount_RNA)
count_s <- sd(seurat@meta.data$nCount_RNA)
count_q <- quantile(seurat@meta.data$nCount_RNA, 0.95)

(plot1
  + geom_hline(yintercept = 5, color = "red"))

(plot2
  + geom_hline(yintercept = 5, color = "red"))

(plot3
  + geom_hline(yintercept = 500, color = "red")
  + geom_vline(xintercept = count_q, color = "red"))

print(paste0("Feature stats:", feature_min, feature_m, feature_max, feature_s));
print(paste0("UMI stats:", count_min, count_m, count_max, count_s, count_q));

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nCount_RNA < count_q & percent.mt < 5)

saveRDS(seurat, file = paste0(output_folder, "2-loaded.rds"))
