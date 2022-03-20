# rm(list = ls())

output_folder <- "./Data/Output/"

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("hdf5r")) {
  install.packages("hdf5r")
  library(hdf5r)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dyplr)
}
if (!require("Matrix")) {
  install.packages("Matrix")
  library(Matrix)
}
if (!require("logging")) {
  install.packages("logging")
  library(logging)
}
if (!require("stringr")) {
  install.packages("stringr")
  library(stringr)
}
if (!require("stringi")) {
  install.packages("stringi")
  library(stringi)
}

load_counts <- function(filename, separator, project, min_cells, min_features) {
  raw_data <- read.table(file = filename, sep = separator)
  print(head(raw_data))
  hist(colSums(raw_data),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression"
  )
  seurat_object <- CreateSeuratObject(counts = raw_data,
                                      min.cells = min_cells,
                                      min.features = min_features,
                                      project = project)
  return(seurat_object)
}

load_hdf5 <- function(filename,
                      project,
                      min_cells,
                      min_features) {
  raw_data <- Read10X_h5(file)
  print(head(raw_data))
  hist(colSums(raw_data),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression"
  )
  seurat_object <- CreateSeuratObject(counts = raw_data,
                                      min.cells = min_cells,
                                      min.features = min_features,
                                      project = project)
  return(seurat_object)
}

load_seurat <- function(filename,
                        separator = FALSE,
                        project = "seurat",
                        min_cells = 100,
                        min_features = 500) {
  if (grepl(".h5", filename, fixed = TRUE)) {
    loginfo(paste('loading hdf5 file'))
    
    return(load_hdf5(filename = filename,
                     project = project,
                     min_cells = min_cells,
                     min_features = min_features))
  } else {
    loginfo(paste('loading regular file'))

    if (!separator) {
      if (grepl(".txt", filename, fixed = TRUE) || grepl(".tsv", filename, fixed = TRUE)) {
        loginfo("assuming separator as = \\t")
        sep <- "\t"
      } else {
        loginfo("assuming separator as = ,")
        sep <- ","
      }
    } else {
      loginfo(paste("using provided separator =", separator))
      sep <- separator
    }

    return(load_counts(filename = filename,
                       separator = sep,
                       project = project,
                       min_cells = min_cells,
                       min_features = min_features))
  }
}

print_img <- function(x) {
  graphics.off()
  print(x)
  if ("rmote" %in% .packages(all.available = TRUE)) {
    if (rmote:::is_rmote_on()) {
      rmote::plot_done()
    }
  }
}

src <- function(stage = 1) {
  if (stage == 1) {
    return(source("Scripts/1-start-rmote.R"))
  } else if (stage == 2) {
    return(source("Scripts/2-data-loading.R"))
  } else if (stage == 3) {
    return(source("Scripts/3-scale-normalize.R"))
  } else if (stage == 4) {
    return(source("Scripts/4-pca-jackstraw.R"))
  } else if (stage == 5) {
    return(source("Scripts/5-biomarkers.R"))
  }
}

graphics.off()

loginfo("setup done!")
