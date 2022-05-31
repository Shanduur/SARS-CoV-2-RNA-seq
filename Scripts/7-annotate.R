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

BiocManager::insestall("SingleR")

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

replace <- Vectorize(function(x) {
  gsub("_", " ", x)
})

for (level in annotationLevels) {
  loginfo(paste("predicted Azimuth annotations on dim plot for", level))
  dpann <- DimPlot(seurat, group.by = level, label = TRUE, label.size = 3, split.by = "smoking") + NoLegend()
  print_img(dpann,
            prefix = prefix,
            title = paste0("dim_plot", level, sep = "-"),
            device = device)

  
  loginfo(paste("predicted Azimuth annotation scores", level))
  vlnScore <- VlnPlot(seurat, split.by = "smoking", features = paste0(level, ".score")) + NoLegend()
  print_img(vlnScore,
            prefix = prefix,
            title = paste0("vln_plt_score", level, sep = "-"),
            device = device)
}

annotationLevels <- c("1", "2", "3", "4", "5", "Finest")

loginfo("creating smokers object")
smokers_cells <- rownames(seurat@meta.data)[seurat@meta.data$smoking == "smokers"]
seurat_smokers <- seurat[,colnames(seurat) %in% smokers_cells]

loginfo("creating non-smokers object")
non_smokers_cells <- rownames(seurat@meta.data)[seurat@meta.data$smoking == "non-smokers"]
seurat_non_smokers <- seurat[,colnames(seurat) %in% non_smokers_cells]

genes_more_expr_smokers <- rowSums(seurat_smokers@assays$RNA@data) > rowSums(seurat_non_smokers@assays$RNA@data)
genes_more_expr_smokers_cnt <- sum(genes_more_expr_smokers)
genes_more_expr_smokers_names <- rownames(seurat)[genes_more_expr_smokers]

genes_less_expr_smokers <- rowSums(seurat_smokers@assays$RNA@data) < rowSums(seurat_non_smokers@assays$RNA@data)
genes_less_expr_smokers_cnt <- sum(genes_less_expr_smokers)
genes_less_expr_smokers_names <- rownames(seurat)[genes_less_expr_smokers]

for (level in annotationLevels) {
  annotationFile <- paste0("./Data/Azimuth/Azimuth_Lung_Annotations_-_Level_", level , ".csv")
  
  loginfo(paste("loading data from file:", annotationFile))
  annos <- read.table(file = annotationFile, sep = ",", header = TRUE)
  annos$Markers <- strsplit(annos$Markers, split = ",")
  
  genes_less_expr_smokers_names_intersection <-c()
  genes_more_expr_smokers_names_intersection <-c()

  i <- 1
  for (id in annos$Label) {
    features <- unlist(annos$Markers[annos$Label == id])
    
    intersection <- intersect(rownames(seurat), features)

    genes_less_expr_smokers_names_intersection[i] <- paste(intersect(genes_less_expr_smokers_names, features), collapse = ",")
    genes_more_expr_smokers_names_intersection[i] <- paste(intersect(genes_more_expr_smokers_names, features), collapse = ",")

    dir.create(paste(output_folder, level, id, sep = "/"), recursive = TRUE)
    
    if (length(intersection) > 0) {
      loginfo(paste("Violin plot for [", paste(intersection, collapse = ", "), "]"))
      vpl <- VlnPlot(seurat,
                     features = intersection,
                     log = TRUE,
                     pt.size = 0,
                     split.by = "smoking")
      print_img(vpl,
                prefix = paste(level, id, prefix, sep = "/"),
                title = paste0("violin_plot", level, id, i, sep = "-"),
                device = device)

      loginfo(paste("Feature plot for [", paste(intersection, collapse = ", "), "]"))
      fpl <- FeaturePlot(seurat,
                         features = intersection,
                         split.by = "smoking")
      print_img(fpl,
                prefix = paste(level, id, prefix, sep = "/"),
                title = paste("feature_plot", level, i, sep = "-"),
                device = device)
    } else {
      logerror(paste("unable to create plots - no genes found in dataset for", id))
    }
    
    i <- i+1
  }
  
  df <- data.frame(
    labels = annos$Label,
    expressed_less_in_smokers = genes_less_expr_smokers_names_intersection,
    expressed_more_in_smokers = genes_more_expr_smokers_names_intersection
  )
  
  loginfo("writting expressed markers for each cluster type based on smoking status")
  write.csv2(df,
              paste0(export_folder, paste0(level, ".csv")),
              sep = ";",
              row.names = FALSE,
              col.names = TRUE,
              quote = TRUE)
}

loginfo("violin plot comparing counts of two groups: cell annotations vs predicted annotations")
vln_P <- VlnPlot(seurat,
                 features = "nCount_RNA",
                 pt.size = 0.001,
                 cols = 2,
                 group.by = "predicted.ann_finest_level",
                 split.by = "smoking") + 
          ggtitle("The total number of molecules detected within a cell per predicted cluster")
vln_A <- VlnPlot(seurat,
                 features = "nCount_RNA",
                 pt.size = 0.001,
                 cols = 2,
                 group.by = "cluster",
                 split.by = "smoking") +
          ggtitle("The total number of molecules detected within a cell per annotated cluster")
print_img(vln_P + vln_A,
          height = 11,
          prefix = prefix,
          title = "violin_predicted_clusters",
          device = device)

loginfo("dim plot comparing two groups: cell annotations vs predicted annotations")
dp_P <- DimPlot(seurat,
                 group.by = c("predicted.ann_finest_level", "cluster"),
                 label = TRUE,
                 label.size = 3,
                 split.by = "smoking")
print_img(dp_P,
          height = 11,
          prefix = prefix,
          title = "dim_plot_predicted_clusters",
          device = device)

dotplt1 <- DotPlot(seurat_smokers,
                  features = sample(rownames(seurat), 20),
                  cols = c("lightgrey", "blue"),
                  group.by = "predicted.ann_finest_level") + RotatedAxis()
dotplt2 <- DotPlot(seurat_non_smokers,
                   features = sample(rownames(seurat), 20),
                   cols = c("lightgrey", "red"),
                   group.by = "predicted.ann_finest_level") + RotatedAxis()
print_img(dotplt1 + dotplt2,
          prefix = prefix,
          title = "dot_plot",
          device = device)


markers <- read.table(paste0(export_folder, "markers.data.txt"), header = TRUE)

loginfo("DotPlot")
dotplt <- DotPlot(seurat,
                  features = unique(markers$gene),
                  group.by= "predicted.ann_finest_level",
                  # split.by = "smoking",
                  dot.scale = 16) +
  RotatedAxis() +
  labs(y="Class", x = "Genes") +
  theme(axis.text = element_text(size=20)) +
  scale_y_discrete(labels = replace) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray",
                                        linetype = "dashed", size=0.35),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  facet_grid(~ unique(seurat@meta.data$smoking)) +
  theme(strip.text = element_text(size = 30))
print_img(dotplt,
          width = 22,
          height = 16,
          prefix = prefix,
          title = "dot_plot",
          device = device)

loginfo(paste("violin plot"))
vln <- VlnPlot(seurat,
                features = unique(markers$gene),
                group.by = "predicted.ann_finest_level",
                slot = "counts",
                log = TRUE,
                split.by = "smoking"
)
print_img(vln,
          prefix = prefix,
          title = "violin",
          device = device,
          width = 25,
          height = 18
)

hm1 <- DoHeatmap(seurat, features = unique(markers$gene)) + NoLegend()
print_img(hm1,
          prefix = prefix,
          title = "heatmap",
          device = device)

saveRDS(seurat, file = paste0(checkpoint_folder, "7-annotate.rds"))
