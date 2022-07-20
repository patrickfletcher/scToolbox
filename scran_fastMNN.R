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
splits <- unlist(mat$splitby)

normpars=mat$normpars[,,1]
# print(normpars)

n_hvg = as.numeric(normpars$n.features)
min.mean = as.numeric(normpars$min.mean)
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)

mnnpars=mat$mnnpars[,,1]
# print(mnnpars)

k <- as.numeric(mnnpars$k)
d <- as.numeric(mnnpars$d)
ndist <- as.numeric(mnnpars$ndist)
#auto.merge

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")

# restrict=cellsubset$restrict

data <- Seurat::Read10X_h5(datafile, use.names = F)
# data<-read10xCounts(datafile)
# data <- as.array(counts(data))
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

#figure out the split/merge order
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
# - seems like a different option for partitioning might work... 

if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  start_time = Sys.time()
  
  clust.pin <- quickCluster(sce, block=block)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min.mean)
  
  end_time = Sys.time()
  print(end_time - start_time)
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

print('Performing fastMNN correction...')
start_time = Sys.time()

fp <-FastMnnParam(k=k, d=d, ndist=ndist, merge.order=merge.order)
sce.mnn <- correctExperiments(sce,
                              batch = block,
                              subset.row=chosen.hvgs,
                              PARAM = fp)
end_time = Sys.time()
print(end_time - start_time)

mnn <- reducedDim(sce.mnn, "corrected")
rot <- as.matrix(rowData(sce.mnn))

print('Writing results to MAT file...')
start_time = Sys.time()

writeMat(resultsmatfile, mnn=mnn, hvgs=chosen.hvgs, rot=rot, sizefactors=sizeFactors(sce))

end_time = Sys.time()
print(end_time - start_time)
