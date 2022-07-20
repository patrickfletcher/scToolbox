# command line script to run Seurat SCTransform normalization + fastMNN
library(Seurat)
library(SeuratWrappers)
library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
splits <- unlist(mat$splitby)


sctpars=mat$sctpars[,,1]
# print(sctpars)

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
# mnn pars
mnnpars=mat$mnnpars[,,1]
# print(mnnpars)

# SelectIntegrationFeatures(
#   object.list,
#   nfeatures = 2000,
#   assay = NULL,
#   verbose = TRUE,
#   fvf.nfeatures = 2000,
#   ...
# )
n.anchor.features = as.numeric(mnnpars$n.anchor.features)

k <- as.numeric(mnnpars$k)
d <- as.numeric(mnnpars$d)
ndist <- as.numeric(mnnpars$ndist)
#auto.merge

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)

# restrict=cellsubset$restrict

data <- Read10X_h5(datafile)
so <- CreateSeuratObject(data, meta.data = cellinfo)
so <- subset(so, keep==1)

#figure out the split/merge order
split1<-unlist(so[[splits[1]]])

merge.cats <- unique(split1)

if (length(splits)==2) {
  split2<-unlist(so[[splits[2]]])
  merge.cats <- split(split1,split2)
  merge.cats <- lapply(merge.cats, unique)
}
# print(merge.cats)

mergeix<-seq(from=1, to=length(merge.cats))
if (length(mnnpars$merge.order)!=0) {
  mergeix <- mnnpars$merge.order
}
# print(mergeix)

merge.order<- merge.cats[c(mergeix)]
print(merge.order)

#must split on most granular split... 
so.list <- SplitObject(so, split.by=splits[1])
remove('so')
gc()

so.list <- lapply(X = so.list, FUN = SCTransform,
                  ncells = n.cells,
                  variable.features.n = variable.features.n,
                  vars.to.regress = vars.to.regress,
                  n_genes = n.genes,
                  bin_size = bin.size,
                  method = vst.method,
                  vst.flavor = vst.flavor,
                  return.only.var.genes=F) #to make sure we can combine?


var.features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = n.anchor.features)

# The Seurat Wrapper version 
so <- RunFastMNN(so.list,
                 assay="SCT",
                 features=var.features,
                 k = k,
                 ndist = ndist,
                 d = d,
                 merge.order=merge.order)

remove('so.list')
gc()
mnn <- Embeddings(so, reduction = "mnn")

writeMat(resultsmatfile, mnn=mnn,  hvgs=var.features)




# #what about combining the SCT assays into one large matrix, then running fastMNN directly?
# sce.list <- lapply(X=so.list, FUN=function(x){as.SingleCellExperiment(x,assay="SCT")})
# 
# universe <- rownames(sce.list[[1]])
# for (i in seq(2,length(sce.list))){
#   universe <- intersect(universe, rownames(sce.list[[i]]))
# }
# for (i in seq(1,length(sce.list))){
#   this <- sce.list[[i]]
#   # this <- this[universe,]
#   this <- this[var.features,]
#   sce.list[[i]] <- this
# }
# sce <- do.call(cbind, sce.list)
# 
# fp <-FastMnnParam(k=k, d=d, ndist=ndist, merge.order=merge.order)
# sce.mnn <- correctExperiments(sce.list,
#                               subset.row=var.features,
#                               PARAM = fp)
# 
# mnn <- reducedDim(sce.mnn, "corrected")
