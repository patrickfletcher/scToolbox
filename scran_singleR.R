#singleR interface for matlab
suppressPackageStartupMessages(library(R.matlab))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SingleR))

args <- commandArgs(TRUE)
parsmatfile <- args[1]
resultsmatfile <- args[2]

####################################
# Notes
####################################
# singleR wraps trainSingleR and classifySingleR, optionally first calling aggregateReference
# Depending on choice of marker genes, also calls pairwiseDE functions (e.g., pairwiseWilcox)

# SingleR(
#   test,
#   ref,
#   labels,
#   clusters = NULL,        <---if clusters provided, aggregates them before assigning labels

#   genes = "de",           <---classic mode: top genes with largest lfc
#   sd.thresh = 1,
#   de.method = "classic",  <---can set user genes, or use test (eg., pairwiseWilcox)
#   de.n = NULL,
#   de.args = list(),   

#   aggr.ref = FALSE,       <---aggregateReference
#   aggr.args = list(),  

#   recompute = TRUE,     *multireference stuff
#   restrict = NULL,

#   quantile = 0.8,         <---singleR parameters
#   fine.tune = TRUE,
#   tune.thresh = 0.05,
#   prune = TRUE,

#   assay.type.test = "logcounts",
#   assay.type.ref = "logcounts",
#   check.missing = TRUE,
#   BNPARAM = KmknnParam(),
#   BPPARAM = SerialParam()
# )
# 
# aggregateReference(
#   ref,
#   labels,
#   ncenters = NULL,
#   power = 0.5,
#   ntop = 1000,
#   assay.type = "logcounts",
#   rank = 20,
#   subset.row = NULL,
#   check.missing = TRUE,
#   BPPARAM = SerialParam(),
#   BSPARAM = bsparam()
# )
# 
# trainSingleR(
#   ref,
#   labels,
#   genes = "de",
#   sd.thresh = 1,
#   de.method = c("classic", "wilcox", "t"),
#   de.n = NULL,
#   de.args = list(),
#   aggr.ref = FALSE,
#   aggr.args = list(),
#   recompute = TRUE,
#   restrict = NULL,
#   assay.type = "logcounts",
#   check.missing = TRUE,
#   BNPARAM = KmknnParam(),
#   BPPARAM = SerialParam()
# )
# 
# classifySingleR(
#   test,
#   trained,
#   quantile = 0.8,
#   fine.tune = TRUE,
#   tune.thresh = 0.05,
#   sd.thresh = NULL,
#   prune = TRUE,
#   assay.type = "logcounts",
#   check.missing = TRUE,
#   BPPARAM = SerialParam()
# )
# 
# pairwiseWilcox(
#   x,
#   groups,
#   block = NULL,
#   restrict = NULL,
#   exclude = NULL,
#   direction = c("any", "up", "down"),
#   lfc = 0,
#   log.p = FALSE,
#   gene.names = NULL,
#   subset.row = NULL,
#   BPPARAM = SerialParam()
# )

####################################
# parse parsmatfile
mat <- readMat(parsmatfile)

# reffile <- "merge_mays.h5"
# datapath <- "D:/scRNAseq/rat/pineal/dev_all5"
# datafile <- "dev_all5.h5"

refinfofile <- unlist(mat$refinfofile) #make sure to add a "keep" column
testinfofile <- unlist(mat$testinfofile) #make sure to add a "keep" column
# reflevs.fine <- c("aPinealocyte", "bPinealocyte", 
#                   "aAstrocyte", "bAstrocyte", "gAstrocyte",
#                   "Endothelial", "VLMCs", 
#                   "aMicroglia", "bMicroglia")
# reflevs.main <- c("Pinealocyte", "Astrocyte", 
#                   "Endothelial", "VLMCs", "Microglia")


# uselabs <- unlist(mat$uselabs)

# load the reference
####################################
refcellinfo <-read.csv(refinfofile, row.names = 'id', header = T)

data <- Seurat::Read10X_h5(file.path(refpath, reffile))
sce.ref <- SingleCellExperiment(assays=list(counts=data), colData=as(refcellinfo, "DataFrame"))
sce.ref$labels.fine <- factor(sce.ref$labels.fine, levels = reflevs.fine)
sce.ref$labels.main <- factor(sce.ref$labels.main, levels = reflevs.main)
sce.ref <- sce.ref[,sce.ref$keep==1]
sce.ref <- sce.ref[,is.na(sce.ref$labels.fine)==F]
sce.ref <- logNormCounts(sce.ref)

####################################
# load the test data
testcellinfo <-read.csv(testinfofile, row.names = 'id')

data <- Seurat::Read10X_h5(file.path(datapath,datafile)) 
sce <- SingleCellExperiment(assays=list(counts=data), colData=as(testcellinfo, "DataFrame"))
sce <- sce[,sce$keep==1]

####################################
# singleR
labs <- sce@colData[,uselabs]


# marker genes
#options for which marker genes to use....
# provided list, classic mode SingleR, new using scran...
de.method <- 'classic'
use.genes <- 'de'


quantile
tune.thresh
aggr.refs

# run singleR
start_time = Sys.time()
predictions <- SingleR(test=sce, ref=sce.ref, labels=labs, de.method=de.method, genes=use.genes, 
                       quantile=quantile, tune.thresh=tune.thresh, aggr.refs = aggr.refs)
end_time = Sys.time()
end_time-start_time
####################################
# store results

# table(predictions$labels)
# sce$labels <- predictions$labels
# sce$pruned.labels <- predictions$pruned.labels

#predictions$scores also useful, and the fine-tuning step (though seems to use 'classic mode' de genes)
