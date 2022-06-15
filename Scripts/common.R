# rm(list = ls())

if (!require("munsell")) {
  install.packages("munsell")
}

if (!require("logging")) {
  install.packages("logging")
  library(logging)
}

output_folder <- "./Data/Output/"
if (!dir.exists(output_folder)) {
  loginfo(paste("creating", output_folder))
  dir.create(output_folder)
} else {
  loginfo(paste(output_folder, "already exists"))
}

export_folder <- "./Data/Exported/"
if (!dir.exists(export_folder)) {
  loginfo(paste("creating", export_folder))
  dir.create(export_folder)
} else {
  loginfo(paste(export_folder, "already exists"))
}

checkpoint_folder <- "./Data/Checkpoints/"
if (!dir.exists(checkpoint_folder)) {
  loginfo(paste("creating", checkpoint_folder))
  dir.create(checkpoint_folder)
} else {
  loginfo(paste(checkpoint_folder, "already exists"))
}

device <- "pdf"

if (!require("hms")) {
  install.packages("hms")
  library(hms)
}
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
if (!require("stringr")) {
  install.packages("stringr")
  library(stringr)
}
if (!require("stringi")) {
  install.packages("stringi")
  library(stringi)
}
if (!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install()
}

# BiocManager::install("limma", ask = FALSE, update = FALSE)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

col_sum_hist <- function(x) {
  hist(colSums(x),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression",
    ylab = "Number of cells"
  )
}

print_img <- function(x,
                      fun = NULL,
                      width = 11,
                      height = 8,
                      title = NULL,
                      device = NULL,
                      pointsize = 12,
                      res = 200,
                      prefix = "00",
                      output_folder = "./Data/Output/") {
  env_rmote <- FALSE
  if ("rmote" %in% .packages(all.available = TRUE)) {
    if (rmote:::is_rmote_on()) {
      logwarn("rmote environment active, all printing is piped to rmote")

      env_rmote <- TRUE
      device <- NA
    }
  }

  if (is.null(title)) {
    title <- format(Sys.time(), "%s")
  }
  title <- basename(title)
  title <- gsub(" ", "_", title)

  graphics.off()

  if (!is.null(device)) {
    if (device == "pdf") {
      graphics_device <- pdf
    } else if (device == "jpeg") {
      graphics_device <- jpeg
      width <- width * 200
      height <- height * 200
    } else if (device == "png") {
      graphics_device <- png
      width <- width * 200
      height <- height * 200
    }

    graphics_device(file = file.path(output_folder, paste0(prefix, "-", title, ".", device)),
          width = width,
          height = height,
          pointsize = pointsize,
          res = res,
          quality = 50)
  }

  if (is.null(fun)) {
    print(x)
  } else {
    fun(x)
  }

  if (!is.null(device)) {
    dev.off()
  }

  if (env_rmote) {
    rmote::plot_done()
  }
}

load_counts <- function(filename,
                        meta,
                        separator,
                        meta_separator,
                        project,
                        quote,
                        min_cells,
                        min_features) {
  loginfo("loading raw expression matrix")
  time_start <- Sys.time()
  raw_data <- read.table(file = filename, sep = separator, quote = quote)
  time_end <- Sys.time()
  loginfo(paste("loading expression matrix took", as_hms(time_end - time_start)))

  if (!is.null(meta) || meta == "") {
    loginfo("loading raw metadata")
    raw_meta <- read.table(meta, header = TRUE, sep = meta_separator, quote = quote)
    raw_meta <- as.data.frame(raw_meta)
    rownames(raw_meta) <- sub("-", ".", raw_meta$X) 
    raw_meta <- subset(raw_meta, select=-c(X))
  } else {
    raw_meta <- NULL
  }

  loginfo("printing histogram")
  print_img(raw_data, fun = col_sum_hist, device = device, title = paste0("00-histogram", filename))

  loginfo("creating seurat object")
  seurat_object <- CreateSeuratObject(counts = raw_data,
                                      min.cells = min_cells,
                                      min.features = min_features,
                                      meta.data = raw_meta,
                                      project = project)

  # seurat_object <- RenameCells(seurat_object, add.cell.id = project)

  loginfo("load_counts done!")
  return(seurat_object)
}

load_hdf5 <- function(filename,
                      project,
                      min_cells,
                      min_features) {
  loginfo(paste("loading expression matrix (h5)", filename))
  raw_data <- Read10X_h5(filename = filename)

  loginfo("printing histogram (h5)")
  print_img(raw_data, fun = col_sum_hist, device = device, title = paste0("00-histogram", filename))

  loginfo("creating seurat object (h5)")
  seurat_object <- CreateSeuratObject(counts = raw_data,
                                      min.cells = min_cells,
                                      min.features = min_features,
                                      project = project)
  loginfo("load_hdf5 done!")
  return(seurat_object)
}

load_seurat <- function(filename,
                        meta = NULL,
                        separator = NULL,
                        meta_separator = ";",
                        project = "seurat",
                        min_cells = 100,
                        min_features = 10,
                        quote = "\"'") {
  if (grepl(".h5", filename, fixed = TRUE)) {
    loginfo(paste("loading hdf5 file"))

    return(load_hdf5(filename = filename,
                     project = project,
                     min_cells = min_cells,
                     min_features = min_features))
  } else {
    loginfo(paste("loading regular file"))

    if (is.null(separator)) {
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
                       meta = meta,
                       separator = sep,
                       meta_separator = meta_separator,
                       project = project,
                       min_cells = min_cells,
                       min_features = min_features,
                       quote = quote))
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
