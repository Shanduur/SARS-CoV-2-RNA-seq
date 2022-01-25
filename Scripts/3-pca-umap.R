source("./Scripts/common.R")

pbmc <- readRDS(file=paste(path.data.out, '2-cleaned.rds', sep=''))

# linear transformation 
# ScaleData():
# - Shifts the expression of each gene, so that the mean expression across cells is 0
# - Scales the expression of each gene, so that the variance across cells is 1
#   - This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Performing linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

saveRDS(pbmc, file = paste(path.data.out, "3-scaled-pca.rds", sep=""))
