# command line script to run multibatchnorm and get batch-adjusted sizefactors
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(R.matlab))

tmp_path <- "D:/tmp/tmp_r_multibatchnorm/"
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
min_mean = as.numeric(normpars$min.mean) #support a list here

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
# block <- cellinfo$sample

print(Sys.time() - start_time)


print('Computing size factors...')
start_time = Sys.time()

SFS = matrix(data=NA_real_,nrow=length(block),ncol=length(min_mean))
for (i in seq(length(min_mean))) {
  mm = min_mean[[i]]
  sce.norm <- multiBatchNorm(sce, batch=block, min.mean=mm)
  SFS[,i] <- sizeFactors(sce.norm)
}

print(Sys.time() - start_time)


print('Writing results to CSV files...')
start_time = Sys.time()

sfsfile<-file.path(tmp_path,paste0(tmp_fileroot,"_sfs.txt"))
write.table(SFS, file=sfsfile, sep=',', row.names = F, col.names = F)

print(Sys.time() - start_time)
