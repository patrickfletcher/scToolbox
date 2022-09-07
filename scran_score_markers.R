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
suppressPackageStartupMessages(library(rmatio))

args <- commandArgs(TRUE)
parsmatfile <- args[1]
resultfile <- args[2]

mat <- readMat(parsmatfile)
# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile)

groupvar <- unlist(mat$groupvar)

batchvar <- NULL
par <- unlist(mat$batchvar)
if (par!="NULL") {
  batchvar<-par
}
print(batchvar)


normpars=mat$normpars[,,1]
# print(normpars)

min.mean = as.numeric(normpars$min.mean)
do_pooledsizefactors = as.logical(normpars$do.pooledsizefactors)
do_multibatch_norm = as.logical(normpars$do.multibatch.norm)


scorepars=mat$scorepars[,,1]
nPairGroups = as.numeric(scorepars$nPairGroups)

# print(scorepars)

pairings <- NULL
par <- unlist(scorepars$pairings)
if (length(par)>1 || (length(par)==1 & par!="NULL")) {
  pairings<-par
}

par2 <- unlist(scorepars$pairings2)
if (length(par2)>1|| (length(par2)==1 & par2!="NULL")) {
  pairings=list(pairings,par2)
}


print(pairings)

lfc <- as.numeric(scorepars$lfc)
full.stats = as.logical(scorepars$full.stats)
min.cells = as.numeric(scorepars$min.cells)

normpars=mat$normpars[,,1]

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
cellinfo <- as(cellinfo, "DataFrame")
# print(head(cellinfo))

#TODO: is the gene info same as matlab????

data <- Seurat::Read10X_h5(datafile)
# data<-read10xCounts(datafile)
# data <- as.array(counts(data))
sce <- SingleCellExperiment(assays=list(counts=data))
colData(sce) <- cellinfo
sce <- sce[,sce$keep==1]

# print(head(colData(sce)))

groups <- sce@colData[,groupvar]

block<-NULL
if (!is.null(batchvar)) {
  block<-sce@colData[,batchvar]
}

print(block)

if (do_pooledsizefactors==T) {
  print('computeSumFactors...')
  clust.pin <- quickCluster(sce, block=block)
  sce <- computeSumFactors(sce, cluster=clust.pin, min.mean=min.mean)
}

if (do_multibatch_norm==T) {
  print('multiBatchNorm...')
  sce <- multiBatchNorm(sce, batch=block)
} else
{
  print('logNormCounts...')
  sce <- logNormCounts(sce)
}


subset.row=rowSums(counts(sce))>=min.cells

# print(head(groups))
# print(head(block))
# print(lfc)
# print(full.stats)
# print(head(subset.row))
print(sum(subset.row==T))

scores <- scoreMarkers(sce, groups=groups, block=block, pairings=pairings, 
                       lfc=lfc, full.stats=full.stats, subset.row=subset.row)

# how to get table into matlab?
# Without full.stats:
# - simple list --> struct with cell type names, but missing column/row names
# - also pass those.
# Full stats needs extra. 
# - each list element has sub-dataframes...

#this handles all cases. Flattens the sub-dataframes into one big one.
# -- except writeMat gets stuck trying to write it...

ids=groups
if (!is.null(block)) {
  ids=DataFrame(group=groups,block=block)
}

summaries  <- summarizeAssayByGroup(logcounts(sce),
                                    ids, 
                                    statistics=c("mean", "prop.detected"),
                                    subset.row=subset.row)
colnames(summaries)

averages <- assay(summaries, "mean")
props <- assay(summaries, "prop.detected")
if (!is.null(batchvar)) {
averages <- correctGroupSummary(averages, group=summaries$group, block=summaries$block)
props <- correctGroupSummary(props, group=summaries$group, block=summaries$block, transform="logit")
}

S=lapply(scores, FUN=function(x) as.data.frame(x))

S$groupnames=c(names(scores))
S$genenames=row.names(scores[[1]])
S$averages= as.data.frame(averages)
S$props= as.data.frame(props)

write.mat(S, resultfile)

# writeMat(resultfile, )
# writeMat(resultfile, scores=S, groupnames=groupnames, genenames=genenames)
