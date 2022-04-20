# command line script to run Seurat integration, then SCTransform + PCA.
library(Seurat)
library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

# parsmatfile <- "C:/Users/fletcherpa/Box/scRNAseq/human_adenoma/matlab/tmp_sct_integrate_pca_pars.mat"
# datafile <- "D:/scRNAseq/human/h5/Aggr_all.h5"
# cellinfofile <- "C:/Users/fletcherpa/Box/scRNAseq/human_adenoma/matlab/tmp_sct_cellsub.csv"

mat <- readMat(parsmatfile)
# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellinfofile <- unlist(mat$cellinfofile) #Required colnames = id, keep
splitby <- unlist(mat$splitby)[1]

sctpars=mat$sctpars[,,1]
# print(sctpars)

intpars=mat$intpars[,,1]
# print(intpars)

outopts=mat$outopts[,,1]
# print(outopts)

# summary of expected parameters:
# vars.to.regress = unlist(sctpars$vars.to.regress)
# variable.features.n = sctpars$variable.features.n
# batch.var = unlist(sctpars$batch.var)
# n.genes = sctpars$n.genes  # Number of genes to use when estimating parameters
# n.cells = sctpars$n.cells  # Number of cells to use when estimating parameters
# bin.size = sctpars$bin.size
# min.cells = sctpars$min.cells
# vst.method = unlist(sctpars$vst.method)
# vst.flavor = unlist(sctpars$vst.flavor) 
# nfeatures = intpars$n.anchor.features
# n.pcs = intpars$n.pcs
# reduction = unlist(intpars$reduction)
# k.anchor = intpars$k.anchor
# k.filter = intpars$k.filter
# k.score = intpars$k.score
# max.features = intpars$max.features
# k.weight = intpars$k.weight
# refs = intpars$refs
# sample.tree = intpars$sample.tree
# out.pcs = outopts$out.pcs
# save.so = outopts$save.so
# so.tag = outopts$so.tag

##################
#sct pars
# SCTransform(
#   object,
#   assay = "RNA",
#   new.assay.name = "SCT",
#   reference.SCT.model = NULL,
#   do.correct.umi = TRUE,
 #   ncells = 5000,                   #Number of subsampling cells used to build NB regression
#   residual.features = NULL,
 #   variable.features.n = 3000,
#   variable.features.rv.th = 1.3,
 #   vars.to.regress = NULL,          #Variables to regress out in a second non-regularized linear regression
#   do.scale = FALSE,
#   do.center = TRUE,
#   clip.range = c(-sqrt(x = ncol(x = object[[assay]])/30), sqrt(x = ncol(x = object[[assay]])/30)),
#   conserve.memory = FALSE,
#   return.only.var.genes = TRUE,
#   seed.use = 1448145,
#   verbose = TRUE,
#   ... #vst pars
# )
vars.to.regress <- NULL
par <- unlist(sctpars$vars.to.regress)
if (par!="NULL") {
  vars.to.regress <- par
}
variable.features.n = as.numeric(sctpars$variable.features.n)

##############
#vst pars - *are overridden by SCTransform
# vst(
#   umi,
#   cell_attr = NULL,
#   latent_var = c("log_umi"),
 #   batch_var = NULL,
#   latent_var_nonreg = NULL,                                 *
 #   n_genes = 2000,
#   n_cells = NULL,                                           *
 #   method = "poisson",
#   do_regularize = TRUE,
#   theta_regularization = "od_factor",
#   res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),    *
 #   bin_size = 500,
 #   min_cells = 5,
#   residual_type = "pearson",
#   return_cell_attr = FALSE,
#   return_gene_attr = TRUE,
#   return_corrected_umi = FALSE,
#   min_variance = -Inf,
#   bw_adjust = 3,
#   gmean_eps = 1,
#   theta_estimation_fun = "theta.ml",
#   theta_given = NULL,
#   exclude_poisson = FALSE,
#   use_geometric_mean = TRUE,
#   use_geometric_mean_offset = FALSE,
#   fix_intercept = FALSE,
#   fix_slope = FALSE,
#   scale_factor = NA,
 #   vst.flavor = NULL,
#   verbosity = 2,
#   verbose = NULL,
#   show_progress = NULL
# )
batch.var <- NULL
par <- unlist(sctpars$batch.var)
if (par!="NULL") {
  batch.var <- par
}
n.genes = as.numeric(sctpars$n.genes)  # Number of genes to use when estimating parameters
n.cells = as.numeric(sctpars$n.cells)  # Number of cells to use when estimating parameters
bin.size = as.numeric(sctpars$bin.size)
min.cells = as.numeric(sctpars$min.cells)
vst.method = unlist(sctpars$vst.method)

vst.flavor <- NULL
par <- unlist(sctpars$vst.flavor)
if (par!="NULL") {
  vst.flavor<-par
}
#When vst.flavor is set to 'v2' sets method = glmGamPoi_offset, n_cells=2000, and exclude_poisson = TRUE


###############
#integration pars

# SelectIntegrationFeatures(
#   object.list,
 #   nfeatures = 2000,
#   assay = NULL,
#   verbose = TRUE,
#   fvf.nfeatures = 2000,
#   ...
# )
nfeatures = as.numeric(intpars$n.anchor.features)

#pca
n.pcs = as.numeric(intpars$n.pcs)

#FindIntegrationAnchors(
# object.list = NULL,
# assay = NULL,
 # reference = NULL,
 # anchor.features = 2000,
# scale = TRUE,
# normalization.method = c("LogNormalize", "SCT"),
# sct.clip.range = NULL,
 # reduction = c("cca", "rpca", "rlsi"),
# l2.norm = TRUE,
# dims = 1:30,
 # k.anchor = 5,
 # k.filter = 200,
 # k.score = 30,
 # max.features = 200,
# nn.method = "annoy",
# n.trees = 50,
# eps = 0,
# verbose = TRUE
# )
reduction = unlist(intpars$reduction)
k.anchor = as.numeric(intpars$k.anchor)
k.filter = as.numeric(intpars$k.filter)
k.score = as.numeric(intpars$k.score)
max.features = as.numeric(intpars$max.features)

#IntegrateData(
# anchorset,
# new.assay.name = "integrated",
# normalization.method = c("LogNormalize", "SCT"),
# features = NULL,
# features.to.integrate = NULL,
# dims = 1:30,
 # k.weight = 100,
# weight.reduction = NULL,
# sd.weight = 1,
 # sample.tree = NULL,
# preserve.order = FALSE,
# eps = 0,
# verbose = TRUE
# )
k.weight = as.numeric(intpars$k.weight)

# reference - A vector specifying the object/s to be used as a reference during integration
refs <- NULL
par <- unlist(intpars$refs)
if (par!="NULL") {
  refs <- par
}

#sample.tree
# Specify the order of integration. Order of integration should be encoded in a
# matrix, where each row represents one of the pairwise integration steps.
# Negative numbers specify a dataset, positive numbers specify the integration
# results from a given row (the format of the merge matrix included in the
# hclust function output). For example: matrix(c(-2, 1, -3, -1), ncol = 2)
# gives:
#
#      [,1]  [,2] 
# [1,]   -2   -3 
# [2,]    1   -1 
# Which would cause dataset 2 and 3 to be integrated first, then the resulting
# object integrated with dataset 1.
sample.tree <- NULL
par <- unlist(intpars$sample.tree)
if (par!="NULL") {
  sample.tree<-par
}

# output options
out.pcs = as.numeric(outopts$out.pcs)
save.sctlist = as.numeric(outopts$save.sctlist)
save.so = as.numeric(outopts$save.so)
so.tag = unlist(outopts$so.tag)

# print(vars.to.regress)
# print(batch.var)
# print(vst.method)
# print(vst.flavor)
# print(refs)
# print(sample.tree)

##################
# cell info table: should include vars used in vars.to.regress, if applicable
cellinfo <-read.csv(cellinfofile, row.names = 'id', header = T)

# load the data
so <- CreateSeuratObject(Read10X_h5(datafile))
so <- AddMetaData(object = so, metadata = cellinfo)
so <- subset(so, subset = keep==1)

so.list <- SplitObject(so, split.by = splitby)

remove(so)
gc()

#batch_var = batch.var,
# min_cells = min.cells,

so.list <- lapply(X = so.list, FUN = SCTransform,
                  ncells = n.cells,
                  variable.features.n = variable.features.n,
                  vars.to.regress = vars.to.regress,
                  n_genes = n.genes,
                  bin_size = bin.size,
                  method = vst.method,
                  vst.flavor = vst.flavor)

#if save.sctlist, save so.list as RDS for later loading/subsetting/new PCA.
if (save.sctlist==1){
  rdsfile<-paste0(paste("sct_split",splitby,so.tag,sep="_"),".rds")
  saveRDS(so.list, file = file.path(dirname(datafile),rdsfile))
}

features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)

so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

so.list <- lapply(X = so.list, FUN = RunPCA, npcs = n.pcs, verbose=FALSE)

anchors <- FindIntegrationAnchors(object.list = so.list,
                                  reduction = reduction,
                                  normalization.method = "SCT",
                                  reference = refs,
                                  dims = 1:n.pcs,
                                  anchor.features = features,
                                  max.features=max.features,
                                  k.anchor=k.anchor,
                                  k.filter=k.filter,
                                  k.score = k.score)

remove(so.list)
gc()

so <- IntegrateData(anchorset = anchors,
                    normalization.method = "SCT",
                    dims = 1:n.pcs,
                    k.weight=k.weight,
                    sample.tree = sample.tree)

remove(anchors)
gc()

DefaultAssay(so) <- "integrated"
#don't do SCTransform again!!!!!
so <- RunPCA(so, verbose=FALSE, npcs = out.pcs)

pca <- Embeddings(so, reduction = "pca")
coeff <- Loadings(so, reduction = "pca")
pc_std <- Stdev(so, reduction = "pca")

writeMat(resultsmatfile, pca=pca, coeff=coeff, pc_std=pc_std, features=features)

#if saveSO, save so as RDS for later loading/subsetting/new PCA.
if (save.so==1){
  rdsfile<-paste0(paste("sct_integrated",reduction,splitby,so.tag,sep="_"),".rds")
  saveRDS(so, file = file.path(dirname(datafile),rdsfile))
}