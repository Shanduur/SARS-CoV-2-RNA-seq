source("./Scripts/common.R")

prefix <- "07"
device <- "jpeg"
# device <- NULL

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("plyr")) {
  install.packages("plyr")
  library(plyr)
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
}

BiocManager::install("SingleR", update = FALSE, ask = FALSE)

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
  dpann <- DimPlot(seurat,
                   group.by = level,
                   split.by = "smoking")
  dpann <- LabelClusters(plot = dpann, id = level, color = "black", box = TRUE)
  print_img(dpann,
            prefix = prefix,
            title = paste0("dim_plot", level, sep = "-"),
            device = device)

  
  loginfo(paste("predicted Azimuth annotation scores", level))
  vlnScore <- VlnPlot(seurat,
                      split.by = "smoking",
                      group.by = level,
                      features = paste0(level, ".score"))
  print_img(vlnScore,
            prefix = prefix,
            title = paste0("vln_plt_score", level, sep = "-"),
            device = device,
            res = 200)
 
  df = as.data.frame(FetchData(seurat, paste0(level, ".score"))[[1]])
  colnames(df) = c('score')
  loginfo(paste("predicted Azimuth annotation scores on histogram", level)) 
  print_img(ggplot(df, aes(x=score)) +
              geom_histogram(bins=10, color="black", fill="grey") +
              theme_bw() +
              scale_y_continuous(trans='log10') +
              ggtitle("Histogram of annotation scores") +
              xlab("Annotation score") +
              ylab("Number of cells"),
            prefix = prefix,
            title = paste0("histogram-score", level, sep = "-"),
            device = device)
}

annotationLevels <- c("1", "2", "3", "4", "5", "Finest")

loginfo("creating smokers object")
smokers_cells <- rownames(seurat@meta.data)[seurat@meta.data$smoking == "smokers"]
seurat_smokers <- seurat[,colnames(seurat) %in% smokers_cells]

loginfo("creating non-smokers object")
non_smokers_cells <- rownames(seurat@meta.data)[seurat@meta.data$smoking == "non-smokers"]
seurat_non_smokers <- seurat[,colnames(seurat) %in% non_smokers_cells]

genes_more_expr_smokers <- rowMeans(seurat_smokers@assays$RNA@data) > rowMeans(seurat_non_smokers@assays$RNA@data)
genes_more_expr_smokers_cnt <- sum(genes_more_expr_smokers)
genes_more_expr_smokers_names <- rownames(seurat)[genes_more_expr_smokers]

genes_less_expr_smokers <- rowMeans(seurat_smokers@assays$RNA@data) < rowMeans(seurat_non_smokers@assays$RNA@data)
genes_less_expr_smokers_cnt <- sum(genes_less_expr_smokers)
genes_less_expr_smokers_names <- rownames(seurat)[genes_less_expr_smokers]

expression_diff <- rowMeans(seurat_smokers@assays$RNA@data) - rowMeans(seurat_non_smokers@assays$RNA@data)

ted <- tail(sort(expression_diff))
ted <- as.data.frame(ted)
ted$ted <- abs(ted$ted)
hed <- head(sort(expression_diff))
hed <- as.data.frame(hed)
hed$hed <- abs(hed$hed)
print(ted)
print(hed)

thed_genes <- unique(c(rownames(ted), rownames(hed)))

expression_diff_df <- as.data.frame(expression_diff)
markers_diff <- as.data.frame(expression_diff_df[unique(markers$gene),])
rownames(markers_diff) <- unique(markers$gene)

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
                device = device,
                res = 200)

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
          device = device,
          res = 200)

loginfo("dim plot comparing two groups: cell annotations vs predicted annotations")
dp_P_1 <- DimPlot(seurat,
                group.by = "predicted.ann_finest_level",
                split.by = "smoking")
dp_P_1 <- LabelClusters(plot = dp_P_1, id = "predicted.ann_finest_level", color = "black", box = TRUE, size = 2)
dp_P_2 <- DimPlot(seurat,
                  group.by = "cluster",
                  split.by = "smoking")
dp_P_2 <- LabelClusters(plot = dp_P_2, id = "cluster", color = "black", box = TRUE, size = 2)
print_img(dp_P_1 + dp_P_2,
          height = 11,
          prefix = prefix,
          title = "dim_plot_predicted_clusters",
          device = device)

pred_anno <- unique(seurat@meta.data$predicted.ann_finest_level)
clus_anno <- unique(seurat@meta.data$cluster)

df <- data.frame(matrix(NA,
                        nrow = 0,
                        ncol = 3))
colnames(df) <- c("Annotated", "Predicted", "count")

for (i in clus_anno) {
  for (j in pred_anno) {
    count <- sum(seurat@meta.data$cluster == i & seurat@meta.data$predicted.ann_finest_level == j)
    df[nrow(df) +1,] <- c(
      i,
      j,
      count
    )
  }
}
df$count <- as.integer(df$count)

lut <- ggplot(df, aes(x = Predicted, y = Annotated, fill = count)) +
  geom_tile(color = "black") +
  geom_text(aes(label = count), color = "black", size = 12) +
  theme(text = element_text(size=40),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(angle=45, hjust=1, size = 24),
        axis.text.y = element_text(hjust=1, size = 24)) +
  scale_fill_gradient(low="aliceblue", high="yellow") +
  coord_fixed()
print_img(lut,
          width = 10,
          height = 11,
          res = 100,
          prefix = prefix,
          title = "look_up_table",
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
          width = 11,
          height = 8,
          res = 100,
          prefix = prefix,
          title = "dot_plot",
          device = device)

loginfo("DotPlot 2")
dotplt2 <- DotPlot(seurat,
                   features = thed_genes,
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
print_img(dotplt2,
          width = 11,
          height = 8,
          res = 100,
          prefix = prefix,
          title = "dot_plot_2",
          device = device)

loginfo("creating covid object")
covid_cells <- rownames(seurat@meta.data)[seurat@meta.data$status == "COVID"]
seurat_covid <- seurat[,colnames(seurat) %in% covid_cells]

loginfo("creating control object")
control_cells <- rownames(seurat@meta.data)[seurat@meta.data$status == "control"]
seurat_control <- seurat[,colnames(seurat) %in% control_cells]

seurat_covid_means <- rowMeans(seurat_covid)
seurat_control_means <- rowMeans(seurat_control)

scmd <- data.frame(genes = rownames(as.data.frame(seurat_covid_means)),
                   avg.covid = as.vector(seurat_covid_means),
                   avg.cntrl = as.vector(seurat_control_means)
)

loginfo("DotPlot 3")
dotplt3 <- DotPlot(seurat,
                   features = thed_genes,
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
  facet_grid(~ unique(seurat@meta.data$status)) +
  theme(strip.text = element_text(size = 30))
print_img(dotplt3,
          width = 11,
          height = 8,
          res = 100,
          prefix = prefix,
          title = "dot_plot_3",
          device = device)

# ------------------------------------------------------------------------------

loginfo(paste("violin plot"))
vln <- VlnPlot(seurat,
              features = unique(markers$gene)[1:8],
              group.by = "predicted.ann_finest_level",
              slot = "counts",
              log = TRUE,
              ncol = 2,
              split.by = "smoking"
)
print_img(vln,
          prefix = prefix,
          title = "violin",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("violin plot cont"))
vln_cont <- VlnPlot(seurat,
                    features = unique(markers$gene)[9:length(unique(markers$gene))],
                    group.by = "predicted.ann_finest_level",
                    slot = "counts",
                    log = TRUE,
                    ncol = 2,
                    split.by = "smoking"
)
print_img(vln_cont,
          prefix = prefix,
          title = "violin-cont",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("violin plot 2"))
vln2 <- VlnPlot(seurat,
                features = unique(markers$gene)[1:8],
                group.by = "cluster",
                slot = "counts",
                log = TRUE,
                ncol = 2,
                split.by = "smoking"
)
print_img(vln2,
          prefix = prefix,
          title = "violin-2",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("violin plot 2 cont"))
vln2_cont <- VlnPlot(seurat,
                     features = unique(markers$gene)[9:length(unique(markers$gene))],
                     group.by = "cluster",
                     slot = "counts",
                     log = TRUE,
                     ncol = 2,
                     split.by = "smoking"
)
print_img(vln2_cont,
          prefix = prefix,
          title = "violin-2-cont",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("violin plot 3"))
vln3 <- VlnPlot(seurat,
                features = thed_genes[1:8],
                group.by = "cluster",
                slot = "counts",
                log = TRUE,
                ncol = 2,
                split.by = "smoking"
)
print_img(vln3,
          prefix = prefix,
          title = "violin-3",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("violin plot 3 cont"))
vln3_cont <- VlnPlot(seurat,
                     features = thed_genes[9:length(thed_genes)],
                     group.by = "cluster",
                     slot = "counts",
                     log = TRUE,
                     ncol = 2,
                     split.by = "smoking"
)
print_img(vln3_cont,
          prefix = prefix,
          title = "violin-3-cont",
          device = device,
          width = 9,
          height = 6,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("violin plot 4"))
vln4 <- VlnPlot(seurat,
                features = thed_genes[1:8],
                group.by = "predicted.ann_finest_level",
                slot = "counts",
                log = TRUE,
                ncol = 2,
                split.by = "smoking"
)
print_img(vln4,
          prefix = prefix,
          title = "violin-4",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("violin plot 4 cont"))
vln4_cont <- VlnPlot(seurat,
                     features = thed_genes[9:length(thed_genes)],
                     group.by = "predicted.ann_finest_level",
                     slot = "counts",
                     log = TRUE,
                     ncol = 2,
                     split.by = "smoking"
)
print_img(vln4_cont,
          prefix = prefix,
          title = "violin-4-cont",
          device = device,
          width = 9,
          height = 6,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("ridge plot"))
rdg <- RidgePlot(seurat,
                features = unique(markers$gene)[1:8],
                group.by = "predicted.ann_finest_level",
                slot = "counts",
                log = TRUE,
                ncol = 2
)
print_img(rdg,
          prefix = prefix,
          title = "ridge",
          device = device,
          width = 9,
          height = 12,
          res = 100
)


loginfo(paste("ridge plot cont"))
rdg_cont <- RidgePlot(seurat,
                      features = unique(markers$gene)[9:length(unique(markers$gene))],
                      group.by = "predicted.ann_finest_level",
                      slot = "counts",
                      log = TRUE,
                      ncol = 2
)
print_img(rdg_cont,
          prefix = prefix,
          title = "ridge-cont",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("ridge plot 2"))
rdg2 <- RidgePlot(seurat,
                  features = unique(markers$gene)[1:8],
                  group.by = "cluster",
                  slot = "counts",
                  log = TRUE,
                  ncol = 2
)
print_img(rdg2,
          prefix = prefix,
          title = "ridge-2",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("ridge plot 2 cont"))
rdg2_cont <- RidgePlot(seurat,
                       features = unique(markers$gene)[9:length(unique(markers$gene))],
                       group.by = "cluster",
                       slot = "counts",
                       log = TRUE,
                       ncol = 2
)
print_img(rdg2_cont,
          prefix = prefix,
          title = "ridge-2-cont",
          device = device,
          width = 9,
          height = 12,
          res = 100
)


# ------------------------------------------------------------------------------

loginfo(paste("ridge plot 3"))
rdg3 <- RidgePlot(seurat,
                  features = thed_genes[1:8],
                  group.by = "cluster",
                  slot = "counts",
                  log = TRUE,
                  ncol = 2
)
print_img(rdg3,
          prefix = prefix,
          title = "ridge-3",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("ridge plot 3 cont"))
rdg3_cont <- RidgePlot(seurat,
                       features = thed_genes[9:length(thed_genes)],
                       group.by = "cluster",
                       slot = "counts",
                       log = TRUE,
                       ncol = 2
)
print_img(rdg3_cont,
          prefix = prefix,
          title = "ridge-3-cont",
          device = device,
          width = 9,
          height = 6,
          res = 100
)

# ------------------------------------------------------------------------------

loginfo(paste("ridge plot 4"))
rdg4 <- RidgePlot(seurat,
                  features = thed_genes[1:8],
                  group.by = "predicted.ann_finest_level",
                  slot = "counts",
                  log = TRUE,
                  ncol = 2
)
print_img(rdg4,
          prefix = prefix,
          title = "ridge-4",
          device = device,
          width = 9,
          height = 12,
          res = 100
)

loginfo(paste("ridge plot 4"))
rdg4_cont <- RidgePlot(seurat,
                       features = thed_genes[9:length(thed_genes)],
                       group.by = "predicted.ann_finest_level",
                       slot = "counts",
                       log = TRUE,
                       ncol = 2
)
print_img(rdg4_cont,
          prefix = prefix,
          title = "ridge-4-cont",
          device = device,
          width = 9,
          height = 6,
          res = 100
)

# ------------------------------------------------------------------------------

hm1 <- DoHeatmap(seurat, features = unique(markers$gene))
print_img(hm1,
          prefix = prefix,
          title = "heatmap",
          device = device)

hm2 <- DoHeatmap(seurat, features = thed_genes, 
                 raster = FALSE, group.by = "smoking")
print_img(hm2,
          prefix = prefix,
          title = "heatmap_2",
          device = device)

covid_sns_all <- c("ACE2", "TMPRSS2", "IL6R", "IL6ST", "PCSK1", "IL6", "PCSK2", "CTSE", "MYRF", "MAG", "MOG", "MBP", "PLP1", "FURIN", "PCSK4", "PCSK5", "PCSK6", "PCSK7", "C1R", "C2", "C3", "C5", "CFI", "CTSS", "CTSL", "CTSB", "CTSC")
covid_sns <- intersect(rownames(seurat), covid_sns_all)

loginfo(paste("covid sns all genes:", length(covid_sns_all)))
loginfo(paste("covid sns expressed genes:", length(covid_sns)))
loginfo(paste("covid sns and thed genes:", length(intersect(thed_genes, covid_sns))))

loginfo(paste("violin plot 5"))
vln5 <- VlnPlot(seurat,
                features = covid_sns,
                group.by = "predicted.ann_finest_level",
                slot = "counts",
                log = TRUE,
                ncol = 2,
                split.by = "smoking"
)
print_img(vln5,
          prefix = prefix,
          title = "violin-5",
          device = device,
          res = 100,
          width = 9,
          height = 9
)

saveRDS(seurat, file = paste0(checkpoint_folder, "7-annotate.rds"))
