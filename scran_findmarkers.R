#script to run scran's get markers

args <- commandArgs(TRUE)

countfile <- args[1] #raw data h5 file
cellinfofile <- args[2]

data <- Read10X_h5(countfile)
cellinfo <-read.csv(cellinfofile, row.names = 'id')

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data), colData=cellinfo)
