# command line script to run Seurat scran normalization + fastMNN
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

args <- commandArgs(TRUE)

tmp_path <- "D:/tmp/tmp_scran_fastMNN/"
tmp_fileroot <- "tmp"
tmp_path <- args[1]
tmp_fileroot <- args[2]

parsmatfile<-file.path(tmp_path,paste0(tmp_fileroot,"_pars.mat"))

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
splits <- unlist(mat$splitby)

normpars=mat$normpars[,,1]
# print(normpars)

do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch = as.logical(normpars$do.multibatch)
min_mean = as.numeric(normpars$min.mean)

hvgpars=mat$hvgpars[,,1]
do_poissonvar = as.logical(hvgpars$do.poissonvar)
do_densityweights = as.logical(hvgpars$do.densityweights)
do_topn = as.logical(hvgpars$do.topn)
n_hvg = as.numeric(hvgpars$n.features)
var_thr = as.numeric(hvgpars$var.thr)
fdr_thr = as.numeric(hvgpars$fdr.thr)


mnnpars=mat$mnnpars[,,1]
# print(mnnpars)

k <- as.numeric(mnnpars$k)
d <- as.numeric(mnnpars$d)
ndist <- as.numeric(mnnpars$ndist)

prop_k <- NULL
if (as.logical(mnnpars$do.propk)) { prop_k = as.numeric(mnnpars$prop.k)}

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")

# restrict=cellsubset$restrict

data <- Seurat::Read10X_h5(datafile, use.names = F)
# data<-read10xCounts(datafile)
# data <- as.array(counts(data))
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

# print(dim(sce))

gene_sub = as.logical(normpars$gene.subset)
# print(sum(gene_sub==T))

if (length(gene_sub)==dim(sce)[1]) {
  sce <- sce[gene_sub==T,]
}

print(dim(sce))


#blocking on most granular grouping...
block <- sce@colData[,splits[1]]

#options? libsize
# - seems like a different option for partitioning might work... 

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

if (do_topn==T) {
  chosen.hvgs <- getTopHVGs(blk,  n=n_hvg, var.threshold = var_thr, fdr.threshold = fdr_thr)
} else
{
  chosen.hvgs <- getTopHVGs(blk, var.threshold = var_thr, fdr.threshold = fdr_thr)
}


# myUnlist <- function(X){
#   U <- sapply(mnnpars$merge.order, unlist, recursive=F)
# }


#figure out the split/merge order
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
  mergeix <- lapply(mergeix, unlist, use.names =F)
}
print(mergeix)

# make_order <- function(groups, or){
#   for ()
# }

# merge.order <- merge.cats[c(mergeix)]
merge.order <- lapply(mergeix, function(X){merge.cats[X]})


# merge.order = NULL
# if (length(mnnpars$merge.order)!=0) {
#   # merge.order <- mnnpars$merge.order
#   merge.order <- sapply(mnnpars$merge.order, unlist)
# }
print(merge.order)


print('Performing fastMNN correction...')
start_time = Sys.time()

fp <-FastMnnParam(k=k, prop.k=prop_k, d=d, ndist=ndist, merge.order=merge.order, get.variance = FALSE)

sce.mnn <- correctExperiments(sce,
                              batch = block,
                              subset.row=chosen.hvgs,
                              PARAM = fp)
end_time = Sys.time()
print(end_time - start_time)

merge.info <- metadata(sce.mnn)$merge.info
pairs <- merge.info$pairs
merge.info$pairs <- NULL
merge.info <- data.frame(merge.info)
merge.info$num.pairs <- sapply(pairs, nrow)
merge.info <- apply(merge.info, 2, as.character)

# pca.info <- metadata(sce.mnn)$pca.info

mnn <- reducedDim(sce.mnn, "corrected")
rot <- as.matrix(rowData(sce.mnn))

print('Writing results to CSV files...')
start_time = Sys.time()

# library(rhdf5)
# h5write(as.matrix(X),"output.h5","X")

infofile<-file.path(tmp_path,paste0(tmp_fileroot,"_minfo.csv"))
mnnfile<-file.path(tmp_path,paste0(tmp_fileroot,"_mnn.csv"))
rotfile<-file.path(tmp_path,paste0(tmp_fileroot,"_rot.csv"))
hvgfile<-file.path(tmp_path,paste0(tmp_fileroot,"_hvgs.txt"))
sfsfile<-file.path(tmp_path,paste0(tmp_fileroot,"_sfs.txt"))

write.csv(merge.info, file=infofile, row.names = F)
write.table(mnn, file=mnnfile, sep=',', row.names = F, col.names = F)
write.table(rot, file=rotfile, sep=',', row.names = F, col.names = F)
# write.table(chosen.hvgs, file=hvgfile, sep=',', row.names = F, col.names = F)
writeLines(chosen.hvgs, con = hvgfile)
write.table(sizeFactors(sce), file=sfsfile, sep=',', row.names = F, col.names = F)

clustfile<-file.path(tmp_path,paste0(tmp_fileroot,"_qclust.txt"))
if (do_pooledsizefactors==T) {
  write.table(clust.pin, file=clustfile, sep=',', row.names = F, col.names = F)
}


# writeMat(resultsmatfile, mnn=mnn, hvgs=chosen.hvgs, rot=rot, sizefactors=sizeFactors(sce))


end_time = Sys.time()
print(end_time - start_time)
