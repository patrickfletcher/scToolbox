# command line script to run Seurat SCTransform normalization, then PCA. Must specify n_pcs...
library(Seurat)

args <- commandArgs(TRUE)

data_file <- args[1] #full path to h5 file
cellsubset_file <- args[2] #two column file with full list of cells. Colnames = id, keep
n_pcs <- as.numeric(args[3])
out_file <- args[4] #output file to write pearson residuals

cellsubset <-read.csv(cellsubset_file, row.names = 'id', header = T)

#as of 8/3/21 - sctransform has a bug with batch_var
# batch<-ifelse("batch" %in% colnames(cellsubset), c("batch"), NULL)

data<-Read10X_h5(data_file, use.names=F)
data <- data[,rownames(cellsubset)]

so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellsubset)

so <- SCTransform(so, method = "glmGamPoi") #sets default assay to SCT , batch_var=batch

so <- RunPCA(so, verbose=FALSE, npcs = n_pcs)

pcs <- Embeddings(so, reduction = "pca")
write.csv(pcs,out_file)


#trying to write counts/scaled.data is prohibitive: too big

# resid <- so[["SCT"]]@scale.data
# write.csv(resid, resid_file)

# library(DropletUtils)
# resid <- so[["SCT"]]@counts
# write10xCounts(path=resid_file, x=resid, overwrite = T,
#                version="3", gene.id = rownames(resid))