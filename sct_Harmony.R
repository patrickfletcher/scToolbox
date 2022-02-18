# command line script to run Seurat SCTransform normalization + fastMNN
library(Seurat)
library(harmony)
library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)

# print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
splitby <- unlist(mat$splitby)

sctpars=mat$sctpars[,,1]
# print(sctpars)

vars.to.regress = unlist(sctpars$vars.to.regress)
if (vars.to.regress=="none") {
  vars.to.regress<-NULL
}
print(vars.to.regress)

n.features = sctpars$n.features
n.anchor.features = sctpars$n.anchor.features
n.pcs = sctpars$n.pcs #default RunPCA: 50
#method
#...


harmonypars=mat$harmonypars[,,1]
# print(mnnpars)

maxit <- harmonypars$maxit
saveit <- harmonypars$saveit
#auto.merge

# print(saveit)

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)
# restrict=cellsubset$restrict

data <- Read10X_h5(datafile)
so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellinfo)
so <- subset(so, keep==1)

assay.use="SCT"

so.list <- SplitObject(so, split.by=splitby)
remove('so')
gc()

so.list <- lapply(X = so.list,
                  FUN = SCTransform,
                  vst.flavor='v2',
                  variable.features.n = n.features,
                  vars.to.regress = vars.to.regress)


var.features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = n.anchor.features)
so <- merge(x = so.list[[1]], y = so.list[2:length(so.list)], merge.data=TRUE)

remove('so.list')
gc()

VariableFeatures(so) <- var.features
DefaultAssay(so) <- "SCT"

so <- RunPCA(so, verbose=FALSE, npcs=n.pcs)

if ( saveit==0 )  {
  
  so <- RunHarmony(so, 
                   group.by.vars=splitby,
                   assay.use = assay.use,
                   theta = NULL,
                   lambda = NULL,
                   sigma = 0.1,
                   nclust = NULL,
                   tau = 0,
                   block.size = 0.05,
                   max.iter.harmony = maxit,
                   max.iter.cluster = 20,
                   epsilon.cluster = 1e-05,
                   epsilon.harmony = 1e-04,
                   plot_convergence = FALSE)
  
  harmony <- Embeddings(so, reduction = "harmony")
  
} else {
  
  harmony=list()
  harmony[[1]] <- Embeddings(so, reduction = "pca")
  
  for (i in seq(maxit)) {
    set.seed(42)
    so <- RunHarmony(so, 
                     group.by.vars=splitby,
                     assay.use = assay.use,
                     theta = NULL,
                     lambda = NULL,
                     sigma = 0.1,
                     nclust = NULL,
                     tau = 0,
                     block.size = 0.05,
                     max.iter.harmony = i,
                     max.iter.cluster = 20,
                     epsilon.cluster = 1e-05,
                     epsilon.harmony = 1e-04,
                     plot_convergence = FALSE)
    
    harmony[[i+1]] <- Embeddings(so, reduction = "harmony")
  }
  harmony=setNames(harmony,paste0("it",seq(from=0,to=maxit)))
}


writeMat(resultsmatfile, harmony=harmony)