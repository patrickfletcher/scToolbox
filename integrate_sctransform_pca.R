# command line script to run Seurat integration, then SCTransform + PCA.
library(Seurat)

args <- commandArgs(TRUE)

data_file <- args[1] #full path to h5 file
cellsubset_file <- args[2] #two column file with full list of cells. Colnames = id, keep
n_pcs <- as.numeric(args[3])
out_file <- args[4] #output file to write pearson residuals
splitby <- args[5]

cellsubset <-read.csv(cellsubset_file, row.names = 'id', header = T)

#as of 8/3/21 - sctransform has a bug with batch_var
# batch<-ifelse("batch" %in% colnames(cellsubset), c("batch"), NULL)

print(data_file)
print(cellsubset_file)
print(n_pcs)
print(out_file)
print(splitby)

data <- Read10X_h5(data_file, use.names=F)

so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellsubset)
so <- subset(so, subset=keep==1)

so.list <- SplitObject(so, split.by=splitby)

remove(so)
gc()

nfeatures <- 3000
ncells <- 5000 #for vst
n_genes <- 2000 #for vst
min_cells <- 5 #for vst
so.list <- lapply(X = so.list, FUN = SCTransform, ncells=ncells, variable.features.n=nfeatures,
                   method = "glmGamPoi", n_genes=n_genes, min_cells=min_cells) #vst args

features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

so.list <- lapply(X = so.list, FUN = RunPCA, npcs = n_pcs, verbose=FALSE)

norm_method<-"SCT"
# int_method <- "cca"
int_method <- "rpca"
k.anchor <- 5
refs <- NULL
anchors <- FindIntegrationAnchors(
  object.list = so.list, reduction = int_method, dims = 1:n_pcs, k.anchor=k.anchor, 
  anchor.features = features, normalization.method = norm_method, reference = refs) 

remove(so.list)
gc()

k.weight <-100
sample.tree <- NULL
so <- IntegrateData(anchorset = anchors, normalization.method = norm_method,
                        dims = 1:n_pcs) #, k.weight=k.weight, sample.tree = sample.tree

remove(anchors)
gc()

DefaultAssay(so) <- "integrated"
so <- SCTransform(so, method = "glmGamPoi") #sets default assay to SCT , batch_var=batch
so <- RunPCA(so, verbose=FALSE, npcs = n_pcs)
pcs <- Embeddings(so, reduction = "pca")

write.csv(pcs,out_file)
