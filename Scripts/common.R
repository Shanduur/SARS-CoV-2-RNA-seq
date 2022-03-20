# rm(list = ls())

output_folder <- "./Data/Output/"

device <- "pdf"

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

colSumHist <- function(x) {
  hist(colSums(x),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression"
  )
}

print_img <- function(x,
                      fun = NA,
                      width = 11,
                      height = 8,
                      title = NA,
                      device = NA,
                      output_folder = "./Data/Output/") {
  envRmote <- FALSE
  if ("rmote" %in% .packages(all.available = TRUE)) {
    if (rmote:::is_rmote_on()) {
      logwarn("rmote environment active, all printing is piped to rmote")

      envRmote <- TRUE
      device <- NA
    }
  }
  
  if (is.na(title)) {
    title <- format(Sys.time(), "%s")
  }
  title <- basename(title)
  title <- gsub(" ", "_", title)

  graphics.off()

  if (!is.na(device)) {
    if (device == "pdf") {
      grDev <- pdf
    } else if (device == "jpeg") {
      grDev <- jpeg
    } else if (device == "png") {
      grDev <- png
    }

    grDev(file = file.path(output_folder, paste0(title, ".", device)),
          width = width,
          height = height)
  }

  if (is.na(fun)) {
    print(x)
  } else {
    fun(x)
  }

  if (!is.na(device)) {
    dev.off()
  }

  if (envRmote) {
    rmote::plot_done()
  }
}

load_counts <- function(filename, separator, project, min_cells, min_features) {
  raw_data <- read.table(file = filename, sep = separator)
  print_img(raw_data, fun = colSumHist, device = device, title = paste0("histogram", filename))
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
  print_img(raw_data, fun = colSumHist, device = device, title = paste0("histogram", filename))
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
    loginfo(paste("loading hdf5 file"))

    return(load_hdf5(filename = filename,
                     project = project,
                     min_cells = min_cells,
                     min_features = min_features))
  } else {
    loginfo(paste("loading regular file"))

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
