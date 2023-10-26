# command line script to run Seurat scran normalization + fastMNN
# suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

tmp_path <- "D:/tmp/tmp_scran_fastMNN/"
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
splits <- unlist(mat$splitby)

# normalization parameters
normpars=mat$normpars[,,1]
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)
doCosNorm = !do_multibatch
min_mean = as.numeric(normpars$min.mean)

# HVG parameters
hvgpars=mat$hvgpars[,,1]
do_poissonvar = as.logical(hvgpars$do.poissonvar)
do_densityweights = as.logical(hvgpars$do.densityweights)
min_mean_hvg = as.numeric(hvgpars$min.mean.hvg)
var_thr = as.numeric(hvgpars$var.thr)
n_hvg = as.numeric(hvgpars$n.features)
if (n_hvg==0) { n_hvg = NULL}
fdr_thr = as.numeric(hvgpars$fdr.thr)
if (fdr_thr==1) { fdr_thr = NULL}

# fastMNN parameters
mnnpars=mat$mnnpars[,,1]
d <- as.numeric(mnnpars$d)
ndist <- as.numeric(mnnpars$ndist)
k <- as.numeric(mnnpars$k)
prop_k = as.numeric(mnnpars$prop.k)
if (prop_k==0) { prop_k = NULL}

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

#blocking on most granular grouping...
block <- sce@colData[,splits[1]]

print(Sys.time() - start_time)
print(dim(sce))

qclust <- sce$qclust
sfs <- sce$sfs

if (is.null(sfs) == F) {
  print('using precomputed size factors')
  sizeFactors(sce) <- sfs
  do_pooledsizefactors = F
  do_multibatch = F ###############################################
}

start_time = Sys.time()

if (do_pooledsizefactors==T) {
  print('computePooledFactors...')
  if (is.null(sce$qclust) == T) {
    qclust <- quickCluster(sce, block=block)
  } else {
    print('using precomputed clusters')
    }
  sce <- computePooledFactors(sce, cluster=qclust, min.mean=min_mean)
}

if (do_multibatch==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=block, min.mean=min_mean)
}

if (do_pooledsizefactors==F && do_multibatch==F) {
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}

print(Sys.time() - start_time)


print('Model gene variance and select HVGs...')
start_time = Sys.time()
if (do_poissonvar==T) {
  blk <- modelGeneVarByPoisson(sce, block=block, min.mean=min_mean_hvg)
} else
{
  blk <- modelGeneVar(sce, block=block, min.mean=min_mean_hvg, density.weights=do_densityweights)
}

chosen.hvgs <- getTopHVGs(blk, n=n_hvg, var.threshold = var_thr, fdr.threshold = fdr_thr)

print(Sys.time() - start_time)


#figure out the split/merge order
# TODO: just go to a list of lists type approach using one batch var

split1 <- sce@colData[,splits[1]]
# split1 <- cellinfo[,splits[1]]

merge.cats <- unique(split1)
if (length(splits)==2) {
  split2<-sce@colData[,splits[2]]
  # split2 <- cellinfo[,splits[2]]
  merge.cats <- split(split1,split2)
  merge.cats <- lapply(merge.cats, unique)
}
print(merge.cats)

mergeix<-seq(from=1, to=length(merge.cats))
if (length(mnnpars$merge.order)!=0) {
  mergeix <- mnnpars$merge.order
  mergeix <- lapply(mergeix, unlist, use.names=F)
}
print(mergeix)

# merge.order <- merge.cats[c(mergeix)]
merge.order <- lapply(mergeix, function(X){merge.cats[X]})

# merge.cats <- unique(block)
# mergeix <- seq(from=1, to=length(merge.cats))
# if (length(mnnpars$merge.order)!=0) {
#   mergeix <- mnnpars$merge.order
#   mergeix <- lapply(mergeix, unlist, recursive=T, use.names=F)
# }
# merge.order <- lapply(mergeix, function(X){merge.cats[X]})
# print(merge.order)


print('Performing fastMNN correction...')
start_time = Sys.time()

# TODO: use multiBatchPCA + reducedMNN and return the mbPCA

fp <-FastMnnParam(k=k, prop.k=prop_k, d=d, ndist=ndist, merge.order=merge.order, cos.norm = doCosNorm)

sce.mnn <- correctExperiments(sce,
                              batch = block,
                              subset.row=chosen.hvgs,
                              PARAM = fp)

print(Sys.time() - start_time)


print('Writing results to CSV files...')
start_time = Sys.time()

merge.info <- metadata(sce.mnn)$merge.info
# pairs <- merge.info$pairs #could try the mnnDeltas thing?
merge.info$pairs <- NULL
merge.info <- data.frame(merge.info)
# merge.info$num.pairs <- sapply(pairs, nrow)
merge.info <- apply(merge.info, 2, as.character)

# pca.info <- metadata(sce.mnn)$pca.info

mnn <- reducedDim(sce.mnn, "corrected")
rot <- as.matrix(rowData(sce.mnn))


infofile<-file.path(tmp_path,paste0(tmp_fileroot,"_minfo.csv"))
mnnfile<-file.path(tmp_path,paste0(tmp_fileroot,"_mnn.csv"))
rotfile<-file.path(tmp_path,paste0(tmp_fileroot,"_rot.csv"))
hvgfile<-file.path(tmp_path,paste0(tmp_fileroot,"_hvgs.txt"))
sfsfile<-file.path(tmp_path,paste0(tmp_fileroot,"_sfs.txt"))

write.csv(merge.info, file=infofile, row.names = F)
write.table(mnn, file=mnnfile, sep=',', row.names = F, col.names = F)
write.table(rot, file=rotfile, sep=',', row.names = F, col.names = F)
writeLines(chosen.hvgs, con = hvgfile)
write.table(sizeFactors(sce), file=sfsfile, sep=',', row.names = F, col.names = F)

clustfile<-file.path(tmp_path,paste0(tmp_fileroot,"_qclust.txt"))
if (do_pooledsizefactors==T) {
  write.table(qclust, file=clustfile, sep=',', row.names = F, col.names = F)
}

print(Sys.time() - start_time)
