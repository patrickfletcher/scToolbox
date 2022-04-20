# command line script to run scran's scoreMarkers

# scoreMarkers(
#   x,
#   groups,
#   block = NULL,
#   pairings = NULL,
#   lfc = 0,
#   row.data = NULL,
#   full.stats = FALSE,
#   subset.row = NULL,
#   BPPARAM = SerialParam()
# )

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(R.matlab))

args <- commandArgs(TRUE)
parsmatfile <- args[1]
resultfile <- args[2]

mat <- readMat(parsmatfile)
# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellinfofile <- unlist(mat$cellinfofile)

groupvar <- unlist(mat$groupvar)

blockvar <- NULL
par <- unlist(mat$blockvar)
if (par!="NULL") {
  blockvar<-par
}

pairings <- NULL
par <- unlist(intpars$pairings)
if (par!="NULL") {
  pairings<-par
}

lfc <- as.numeric(mat$lfc)
row.data = NULL
full.stats = FALSE
subset.row = NULL

normpars=mat$normpars[,,1]


#testing: hard paths
datafile <- "D:/scRNAseq/human/h5/Aggr_all.h5"
# cellinfofile <- "C:/Users/fletcherpa/Box/scRNAseq/human_adenoma/matlab/tmp_sct_cellsub.csv"
cellinfofile <- "C:/Users/fletcherpa/Box/scRNAseq/human_adenoma/matlab/integrated_6P/p1p2p3_cellinfo_042021.csv"

cellinfo <-read.csv(cellinfofile, header = T)

keep <- cellinfo$cleaned_SL_AmbUnc
cellinfo$keep=keep
groupvar <- "type_cleaned_SL_AmbUnc"
blockvar <- "patient"
lfc <- 1.5


data <- Seurat::Read10X_h5(datafile)
sce <- SingleCellExperiment(assays=list(counts=data), colData=as(cellinfo, "DataFrame"))
sce <- sce[,sce$keep==1]

groups <- sce@colData[,groupvar]
block<-sce@colData[,blockvar]

#options? libsize
clust <- quickCluster(sce, block=block)
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- multiBatchNorm(sce, batch=block)

scores <- scoreMarkers(sce, groups, block, pairings, lfc, full.stats=T)

# how to get table into matlab? 
# - extract data as matrices, column names, row names. 
# - reconstruct on the matlab side...


# writeMat()
