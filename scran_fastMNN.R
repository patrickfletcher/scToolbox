# command line script to run Seurat SCTransform normalization + fastMNN
library(scran)
library(batchelor)
library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
splits <- unlist(mat$splitby)

normpars=mat$normpars[,,1]

n_hvg = normpars$n.features

mnnpars=mat$mnnpars[,,1]
# print(mnnpars)

k <- mnnpars$k
d <- mnnpars$d
ndist <- mnnpars$ndist
#auto.merge

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)

# restrict=cellsubset$restrict

data <- Seurat::Read10X_h5(datafile)
sce <- SingleCellExperiment(assays=list(counts=data), colData=as(cellinfo, "DataFrame"))
sce <- sce[,sce$keep==1]

split1<-sce@colData[,splits[1]]

merge.cats <- unique(split1)
if (length(splits)==2) {
  split2<-sce@colData[,splits[2]]
  merge.cats <- split(split1,split2)
  merge.cats <- lapply(merge.cats, unique)
}
print(merge.cats)

mergeix<-seq(from=1, to=length(merge.cats))
if (length(mnnpars$merge.order)!=0) {
  mergeix <- mnnpars$merge.order
}
print(mergeix)

merge.order<- merge.cats[c(mergeix)]
print(merge.order)

#blocking on most granular grouping...
block=split1

#options? libsize
clust.pin <- quickCluster(sce, block=block)
sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=0.1)
sce <- multiBatchNorm(sce, batch=block)

blk <- modelGeneVarByPoisson(sce, block=block)
# blk <- modelGeneVar(sce, block=block)

chosen.hvgs <- getTopHVGs(blk, n=n_hvg)

fp <-FastMnnParam(k=k, d=d, ndist=ndist, merge.order=merge.order)
sce.mnn <- correctExperiments(sce,
                              batch = block,
                              subset.row=chosen.hvgs,
                              PARAM = fp)

mnn <- reducedDim(sce.mnn, "corrected")
writeMat(resultsmatfile, mnn=mnn, hvgs=chosen.hvgs, sizefactors=sizeFactors(sce))