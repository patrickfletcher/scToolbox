# command line script to run scran multibatchnorm and get batch-adjusted sizefactors
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

args <- commandArgs(TRUE)

tmp_path <- args[1]
tmp_fileroot <- args[2]

parsmatfile<-file.path(tmp_path,paste0(tmp_fileroot,"_pars.mat"))

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
batchvar <- unlist(mat$batchvar)
# clustvar =  unlist(mat$clustvar) #provide external clustering for pooledsfs

normpars=mat$normpars[,,1]
# print(normpars)

do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)
min_mean = as.numeric(normpars$min.mean)

hvgpars=mat$hvgpars[,,1]
do_poissonvar = as.logical(hvgpars$do.poissonvar)
do_densityweights = as.logical(hvgpars$do.densityweights)

var_thr = as.numeric(hvgpars$var.thr)
n_hvg = as.numeric(hvgpars$n.features)
if (n_hvg==0) { n_hvg = NULL}
fdr_thr = as.numeric(hvgpars$fdr.thr)
if (fdr_thr==1) { fdr_thr = NULL}

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")

data <- Seurat::Read10X_h5(datafile, use.names = F)
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

gene_sub = as.logical(normpars$gene.subset)
if (length(gene_sub)==dim(sce)[1]) {
  sce <- sce[gene_sub==T,]
}

print(dim(sce))

block<-sce@colData[,batchvar]

clust.pin <- NULL
if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  start_time = Sys.time()
  
  clust.pin <- quickCluster(sce, block=block)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min_mean)
  
  end_time = Sys.time()
  print(end_time - start_time)
}

if (do_multibatch==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=block, min.mean=min_mean)
} else
{
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}

if (do_poissonvar==T) {
  blk <- modelGeneVarByPoisson(sce, block=block, min.mean=min_mean)
} else
{
  blk <- modelGeneVar(sce, block=block, min.mean=min_mean, density.weights=do_densityweights)
}


chosen.hvgs <- getTopHVGs(blk, n=n_hvg, var.threshold = var_thr, fdr.threshold = fdr_thr)


print('Writing results to CSV files...')
start_time = Sys.time()

hvgfile<-file.path(tmp_path,paste0(tmp_fileroot,"_hvgs.txt"))
sfsfile<-file.path(tmp_path,paste0(tmp_fileroot,"_sfs.txt"))

writeLines(chosen.hvgs, con = hvgfile)
write.table(sizeFactors(sce), file=sfsfile, sep=',', row.names = F, col.names = F)

if (do_pooledsizefactors==T) {
  clustfile<-file.path(tmp_path,paste0(tmp_fileroot,"_qclust.txt"))
  write.table(clust.pin, file=clustfile, sep=',', row.names = F, col.names = F)
}

end_time = Sys.time()
print(end_time - start_time)
