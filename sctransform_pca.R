# command line script to run Seurat SCTransform normalization, then PCA. Must specify n_pcs...
library(Seurat)

args <- commandArgs(TRUE)

data_file <- args[1] #full path to h5 file
cellsubset_file <- args[2] #two column file with full list of cells. Colnames = id, keep
n_pcs <- as.numeric(args[3])
out_file <- args[4] #output file to write pearson residuals
use_so <- args[5] #if 1, interpret data_file as a seurat object saved in RDS file
genesubset_file <- args[6]

cellsubset <-read.csv(cellsubset_file, row.names = 'id', header = T)
if (genesubset_file!=1) {
  genesubset <-read.csv(genesubset_file, header = T)
}

#as of 8/3/21 - sctransform has a bug with batch_var
# batch<-ifelse("batch" %in% colnames(cellsubset), c("batch"), NULL)

#vars to regress are all other columns besides "keep"
# vars.to.regress=colnames(cellsubset)

if (use_so==1){
  so <- readRDS(data_file)
  so <- AddMetaData(object = so, metadata = cellsubset)
  so <- subset(so, keep==1)
  if (genesubset_file!=1) {
    genesubset<-genesubset[genesubset$name%in%VariableFeatures(so),]
    print(length(genesubset$name))
    so <- so[VariableFeatures(so)%in%genesubset$name,]
    print(so)
    so <- so[genesubset$keep==1,]
    print(so)
  }
} else {
  
  data <- Read10X_h5(data_file)
  data <- data[,cellsubset$keep==1]
  if (genesubset_file!=1) {
    data <- data[genesubset$keep==1,]
  }
  so <- CreateSeuratObject(data)
  so <- SCTransform(so, vst.flavor='v2') #sets default assay to SCT , batch_var=batch
}

so <- RunPCA(so, verbose=FALSE, npcs = n_pcs)

pcs <- Embeddings(so, reduction = "pca")
write.csv(pcs, out_file)

filename<-basename(out_file)
resultpath<-dirname(out_file)

loadings <- Loadings(so, reduction = "pca")
write.csv(loadings, file.path(resultpath,paste0("feature_loadings_",filename)), row.names = T)

std <- Stdev(so, reduction = "pca")
write.csv(std, file.path(resultpath,paste0("stdev_",filename)), row.names = T)
