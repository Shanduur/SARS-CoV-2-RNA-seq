if (!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install()
}

BiocManager::install(c("GEOquery"))
library("GEOquery")

gse <- getGEO(filename = "GSE122960_series_matrix.txt",
        GSEMatrix = TRUE,
        getGPL = FALSE)

x <- phenoData(gse)
pData(x)

dir.create(file.path("RAW/"), showWarnings = TRUE)
dir.create(file.path("Filtred/"), showWarnings = TRUE)

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

for (i in 1:length(gse$title)) {
  a <- split_path(gse$supplementary_file_1[i])
  b <- split_path(gse$supplementary_file_2[i])
  print(a[1])

  f1 <- paste0("Filtred/", a[1])
  if (!file.exists(f1)) {
    download.file(gse$supplementary_file_1[i], f1)
  }

  print(b[1])
  f2 <- paste0("RAW/", b[1])
  if (!file.exists(f2)) {
    download.file(gse$supplementary_file_2[i], f2)
  }

  print("====================")
}