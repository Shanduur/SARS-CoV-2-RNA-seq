rm(list = ls())

print(Sys.getpid())

#source("./Scripts/1-start-rmote.R")
source("./Scripts/2-data-loading.R")
source("./Scripts/3-normalize.R")
source("./Scripts/4-scale.R")
source("./Scripts/5-dim-reduction.R")
source("./Scripts/6-biomarkers.R")
