cnt.file <- "./Data/Pneumonia/GSE164948_covid_control_RNA_counts.csv"
meta.file <- "./Data/Pneumonia/GSE164948_covid_control_count_metadata.csv"

cnt.raw <- read.table(file = cnt.file, sep = ',')
meta.raw <- read.table(file = meta.file, sep = ';', header = TRUE)

head(colnames(cnt.raw))

meta.df <- as.data.frame(meta.raw)
rownames(meta.df) <-  sub("-", ".", meta.df$X) 
meta.df <- subset(meta.df, select=-c(X))

head(meta.df)

length(rownames(meta.df))
length(colnames(cnt.raw))

intersect(colnames(cnt.raw), rownames(meta.df))

colnames(cnt.raw)[1]
rownames(meta.df)[1]

so <- CreateSeuratObject(counts = cnt.raw,
                         min.cells = 500,
                         min.features = 1,
                         meta.data = meta.df,
                         project = 'project')

head(so@meta.data$cluster)