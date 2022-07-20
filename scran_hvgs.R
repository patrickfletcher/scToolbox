# command line script to run scran sizefactor, multibatchnorm, modelgenevar
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
batch <- unlist(mat$batch)

normpars=mat$normpars[,,1]
# print(normpars)

min.mean = as.numeric(normpars$min.mean)
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)

hvgpars=mat$hvgpars[,,1]
# print(hvgpars)
n_hvg = as.numeric(hvgpars$n.features)

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")

data <- Seurat::Read10X_h5(datafile)
# data<-read10xCounts(datafile)
# data <- as.array(counts(data))
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

#blocking on most granular grouping...
block=sce@colData[,batch]

#options? libsize
# - seems like a different option for partitioning might work... 

if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  clust.pin <- quickCluster(sce, block=block)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min.mean)
}

if (do_multibatch==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=block)
} else
{
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}

blk <- modelGeneVarByPoisson(sce, block=block)
# blk <- modelGeneVar(sce, block=block)

chosen.hvgs <- getTopHVGs(blk, n=n_hvg)

writeMat(resultsmatfile, hvgs=chosen.hvgs, sizefactors=sizeFactors(sce))