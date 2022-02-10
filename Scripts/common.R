rm(list = ls())

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

load_counts <- function(filename, separator) {
  raw_data <- read.table(file = filename, sep = separator)
  hist(colSums(raw_data),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression"
  )
  seurat_object <- CreateSeuratObject(raw_data)
  return(seurat_object)
}

load_hdf5 <- function(filename) {
  raw_data <- Read10X_h5(file)
  hist(colSums(raw_data),
    breaks = 100,
    main = "Expression sum per cell",
    xlab = "Sum expression"
  )
  seurat_object <- CreateSeuratObject(raw_data)
  return(seurat_object)
}

load_seurat <- function(filename, separator = ";") {
  if (grepl(filename, ".h5", fixed = TRUE)) {
    return(load_hdf5(filename))
  } else {
    return(load_counts(filename = filename, separator = separator))
  }
}