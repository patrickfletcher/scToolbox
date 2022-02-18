# command line script to run Seurat SCTransform normalization + fastMNN
library(Seurat)
library(SeuratWrappers)
library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)

print(mat)

datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile) #two column file with full list of cells. Colnames = id, keep
splitby <- unlist(mat$splitby)


sctpars=mat$sctpars[,,1]
print(sctpars)

vars.to.regress = unlist(sctpars$vars.to.regress)
if (vars.to.regress=="none") {
  vars.to.regress<-NULL
}
print(vars.to.regress)

n.features = sctpars$n.features
n.anchor.features = sctpars$n.anchor.features
#method
#...


mnnpars=mat$mnnpars[,,1]
# print(mnnpars)

k <- mnnpars$k
d <- mnnpars$d
ndist <- mnnpars$ndist
#auto.merge

#if 0<k<1, use prop.k
prop.k<-NULL
if (k<1) {
  prop.k<-k
  k<-NULL
}

cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)

# restrict=cellsubset$restrict

data <- Read10X_h5(datafile)
so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellinfo)
so <- subset(so, keep==1)

so.list <- SplitObject(so, split.by=splitby)
remove('so')
gc()
so.list <- lapply(X = so.list,
                  FUN = SCTransform,
                  vst.flavor='v2',
                  variable.features.n = n.features,
                  vars.to.regress = vars.to.regress)


var.features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = n.anchor.features)

# restrict = restrict,
so <- RunFastMNN(so.list,
                 features=var.features,
                 k = k, #increase to make integration stronger?
                 prop.k = prop.k,  #use a proportion of cells instead of K
                 ndist = ndist,
                 d = d,
                 auto.merge = TRUE)

remove('so.list')
gc()
mnn <- Embeddings(so, reduction = "mnn")
writeMat(resultsmatfile, mnn=mnn)