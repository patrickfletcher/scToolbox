# command line script to run Seurat SCTransform normalization + fastMNN
suppressPackageStartupMessages(library(DropletUtils))
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
batchvar <- unlist(mat$batchvar)
# print(batchvar)

normpars=mat$normpars[,,1]
# print(normpars)

n_hvg = as.numeric(normpars$n.features)
min.mean = as.numeric(normpars$min.mean)
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch_norm = as.logical(normpars$do.multibatch.norm)

pcapars=mat$pcapars[,,1]
# print(pcapars)

d <- as.numeric(pcapars$d)
w <- as.numeric(pcapars$w)

w <- NULL
par <- unlist(pcapars$w)
if (par!="NULL") {
  w <- par
}

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")

# restrict=cellsubset$restrict

data <- Seurat::Read10X_h5(datafile)
# data<-read10xCounts(datafile)
# data <- as.array(counts(data))
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

# dim(sce)

batch=sce@colData[,batchvar]

# head(batch)

#options? libsize
# - seems like a different option for partitioning might work... 

if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  clust.pin <- quickCluster(sce, block=batch)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min.mean)
}

if (do_multibatch_norm==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=batch)
} else
{
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}

print('modelGeneVarByPoisson...')
blk <- modelGeneVarByPoisson(sce, block=batch)
# blk <- modelGeneVar(sce, block=batch)

chosen.hvgs <- getTopHVGs(blk, n=n_hvg)

print('multiBatchPCA...')
pca=multiBatchPCA(sce, batch=batch, d=d, weights=w, preserve.single=T)

pca_extra=pca@metadata

writeMat(resultsmatfile, pca=unlist(pca), pca_extra=pca_extra, hvgs=chosen.hvgs, sizefactors=sizeFactors(sce))
# writeMat(resultsmatfile, pca=unlist(pca), hvgs=chosen.hvgs, sizefactors=sizeFactors(sce))