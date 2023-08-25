# command line script to run multibatchnorm, modelGeneVar, and getTopHGVs
# suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

tmp_path <- "D:/tmp/tmp_scran_multibatch_hvg/"
tmp_fileroot <- "tmp"
args <- commandArgs(TRUE)
tmp_path <- args[1]
tmp_fileroot <- args[2]

print('Loading data...')
start_time = Sys.time()

parsmatfile<-file.path(tmp_path,paste0(tmp_fileroot,"_pars.mat"))
mat <- readMat(parsmatfile)
datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
batchvar <- unlist(mat$batchvar)

normpars=mat$normpars[,,1]
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)
min_mean = as.numeric(normpars$min.mean)

hvgpars=mat$hvgpars[,,1]
do_poissonvar = as.logical(hvgpars$do.poissonvar)
do_densityweights = as.logical(hvgpars$do.densityweights)
min_mean_hvg = as.numeric(hvgpars$min.mean.hvg)

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
block<-sce@colData[,batchvar]


print(Sys.time() - start_time)
print(dim(sce))


clust.pin <- NULL
if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  start_time = Sys.time()
  
  clust.pin <- quickCluster(sce, block=block)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min_mean)
  
  print(Sys.time() - start_time)
}


start_time = Sys.time()

if (do_multibatch==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=block, min.mean=min_mean)
} else
{
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}

print(Sys.time() - start_time)


print('modelGeneVar + getTopHVGs...')
start_time = Sys.time()

MEAN = matrix(data=NA_real_,nrow=n_hvg,ncol=length(min_mean_hvg))
VAR = matrix(data=NA_real_,nrow=n_hvg,ncol=length(min_mean_hvg))
BIO = matrix(data=NA_real_,nrow=n_hvg,ncol=length(min_mean_hvg))
TECH = matrix(data=NA_real_,nrow=n_hvg,ncol=length(min_mean_hvg))
for (i in seq(length(min_mean_hvg))) {
mm <- min_mean_hvg[[i]]
if (do_poissonvar==T) {
  blk <- modelGeneVarByPoisson(sce, block=block, min.mean=min_mean)
} else
{
  blk <- modelGeneVar(sce, block=block, min.mean=mm, density.weights=do_densityweights)
}

chosen.hvgs <- getTopHVGs(blk, n=n_hvg, var.threshold = var_thr, fdr.threshold = fdr_thr)

chosen.stats = blk[chosen.hvgs,]
MEAN[,i] = chosen.stats$mean
VAR[,i] = chosen.stats$total
BIO[,i] = chosen.stats$bio
TECH[,i] = chosen.stats$tech
}

print(Sys.time() - start_time)


print('Writing results to CSV files...')
start_time = Sys.time()

hvgfile<-file.path(tmp_path,paste0(tmp_fileroot,"_hvgs.txt"))
sfsfile<-file.path(tmp_path,paste0(tmp_fileroot,"_sfs.txt"))

writeLines(chosen.hvgs, con = hvgfile)
write.table(sizeFactors(sce), file=sfsfile, sep=',', row.names = F, col.names = F)

meanfile<-file.path(tmp_path,paste0(tmp_fileroot,"_mean.txt"))
varfile<-file.path(tmp_path,paste0(tmp_fileroot,"_var.txt"))
biofile<-file.path(tmp_path,paste0(tmp_fileroot,"_bio.txt"))
techfile<-file.path(tmp_path,paste0(tmp_fileroot,"_tech.txt"))
write.table(MEAN, file=meanfile, sep=',', row.names = F, col.names = F)
write.table(VAR, file=varfile, sep=',', row.names = F, col.names = F)
write.table(BIO, file=biofile, sep=',', row.names = F, col.names = F)
write.table(TECH, file=techfile, sep=',', row.names = F, col.names = F)

if (do_pooledsizefactors==T) {
  clustfile<-file.path(tmp_path,paste0(tmp_fileroot,"_qclust.txt"))
  write.table(clust.pin, file=clustfile, sep=',', row.names = F, col.names = F)
}

print(Sys.time() - start_time)
