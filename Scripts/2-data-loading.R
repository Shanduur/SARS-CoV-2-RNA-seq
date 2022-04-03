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

# samples <- c("cap-ctrl", "covid-ctrl", "covid-cap")
samples <- c("smokers", "non-smokers")
files <- c(
  "./Data/SARS-COV-2/NonSmokers/internal_nonsmokerslung.expression.txt",
  "./Data/SARS-COV-2/Smokers/internal_smokerslung.expression.txt"
  # "./Data/Pneumonia/GSE164948_cap_control_RNA_counts.csv",
  # "./Data/Pneumonia/GSE164948_covid_control_RNA_counts.csv",
  # "./Data/Pneumonia/GSE164948_covid_cap_RNA_counts.csv"
)

meta <- c(
  "./Data/SARS-COV-2/NonSmokers/internal_nonsmokerslung.meta.txt",
  "./Data/SARS-COV-2/Smokers/internal_smokerslung.meta.txt"
  # "./Data/Pneumonia/GSE164948_cap_control_count_metadata.csv",
  # "./Data/Pneumonia/GSE164948_covid_control_count_metadata.csv",
  # "./Data/Pneumonia/GSE164948_covid_cap_count_metadata.csv"
)

seurat_list <- list()

for (i in 1:length(files)) {
  loginfo(paste("loading file:", files[i], "|", meta[i]))
  seurat_list[[i]] <- load_seurat(
    filename = files[i],
    project = samples[i],
    meta = meta[i],
    meta_separator = '\t'
  )
  seurat_list[[i]][["DataSet"]] <- samples[i]
}

loginfo(paste("building seurat object"))
seurat <- merge(
  x = seurat_list[[1]],
  y = unlist(seurat_list)[2:length(seurat_list)],
  add.cell.ids = samples,
  project = "SARS-CoV-2-RNA-seq"
)
str(seurat@meta.data)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RB-")
head(seurat[["percent.mt"]], n = 20)

loginfo(paste("histogram 1: nFeature_RNA"))
print_img(table(seurat[["nFeature_RNA"]]),
          fun = hist,
          title = "nFeature_RNA",
          device = device)

loginfo(paste("histogram 2: percent.mt"))
print_img(table(seurat[["percent.mt"]]),
          fun = hist,
          title = "percent-mt",
          device = device)

loginfo(paste("histogram 3: percent.rb"))
print_img(table(seurat[["percent.rb"]]),
          fun = hist,
          title = "percent-rb",
          device = device)

loginfo(paste("histogram 4: nCount_RNA"))
print_img(table(seurat[["nCount_RNA"]]),
          fun = hist,
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
          title = "violin-qc-metrics-1",
          device = device)

loginfo(paste("violin plot 2"))
vln2 <- VlnPlot(seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 2,
  pt.size = 0
)
print_img(vln2,
          title = "violin-qc-metrics-2",
          device = device)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot3 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

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
          title = "plot1-percent.mt",
          device = device)

loginfo(paste("plot2 - percent.rb"))
print_img(plot2 + geom_hline(yintercept = 5, color = "red"),
          title = "plot2 - percent.rb",
          device = device)

loginfo(paste("plot3 - nFeature_RNA"))
print_img(plot3 +
            geom_hline(yintercept = 500, color = "red") +
            geom_vline(xintercept = count_q, color = "red"),
          title = "plot3-nFeature_RNA",
          device = device)

loginfo(paste("Feature stats:", feature_min, feature_m, feature_max, feature_s))
loginfo(paste("UMI stats:", count_min, count_m, count_max, count_s, count_q))
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nCount_RNA < count_q & percent.mt < 5)

saveRDS(seurat, file = paste0(checkpoint_folder, "2-loaded.rds"))
rm(seurat_list)
