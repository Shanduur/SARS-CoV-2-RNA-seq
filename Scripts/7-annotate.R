source("./Scripts/common.R")

prefix <- "07"
device <- "pdf"
# device <- NULL

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
}

BiocManager::install("SingleR")

if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("satijalab/azimuth", ref = "master")

library(Azimuth)

seurat <- readRDS(file = paste0(checkpoint_folder, "6-biomarkers.rds"))

seurat <- RunAzimuth(seurat, reference = "./Data/Azimuth/")

annotationLevels <- c(
  "predicted.ann_level_1",
  "predicted.ann_level_2",
  "predicted.ann_level_3",
  "predicted.ann_level_4",
  "predicted.ann_level_5",
  "predicted.ann_finest_level"
)

for (level in annotationLevels) {
  loginfo(paste("predicted Azimuth annotations on dim plot for [", paste(intersection, collapse = ', '), "]"))
  dpann <- DimPlot(seurat, group.by = level, label = TRUE, label.size = 3) + NoLegend()
  print_img(dpann,
            prefix = prefix,
            title = paste0("dim_plot", level, sep = "-"),
            device = device)
}

annotationLevels <- c("1", "2", "3", "4", "5", "Finest")

for (level in annotationLevels) {
  annotationFile <- paste0("./Data/Azimuth/Azimuth_Lung_Annotations_-_Level_", level ,".csv")
  
  loginfo(paste("loading data from file:", annotationFile))
  annos <- read.table(file = annotationFile, sep = ',', header = TRUE)
  annos$Markers <- strsplit(annos$Markers, split = ',')

  i <- 0
  for (id in annos$Label) {
    features <- unlist(annos$Markers[annos$Label == id])
    
    intersection <- intersect(rownames(seurat), features)
    
    dir.create(paste(output_folder, level, id, sep = "/"), recursive = TRUE)
    
    if (length(intersection) > 0) {
      loginfo(paste("Violin plot for [", paste(intersection, collapse = ', '), "]"))
      vpl <- VlnPlot(seurat, features = intersection, log = TRUE, pt.size = 0)
      print_img(vpl,
                prefix = paste(level, id, prefix, sep = "/"),
                title = paste0("violin_plot", level, id, i, sep = "-"),
                device = device)
      
      loginfo(paste("Feature plot for [", paste(intersection, collapse = ', '), "]"))
      fpl <- FeaturePlot(seurat, features = intersection)
      print_img(fpl,
                prefix = paste(level, id, prefix, sep = "/"),
                title = paste("feature_plot", level, i, sep = "-"),
                device = device)
    } else {
      logerror(paste("unable to create plots - no genes found in dataset for", id))
    }
    
    i <- i+1
  }
}

saveRDS(seurat, file = paste0(checkpoint_folder, "7-annotate.rds"))

