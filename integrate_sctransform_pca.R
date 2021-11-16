# command line script to run Seurat integration, then SCTransform + PCA.

library(Seurat)

args <- commandArgs(TRUE)

#TODO: move to R.matlab - pass struct of params???
# input_mat_file <- args[1]
# inputs <- readMat(input_mat_file)
# output_mat_file <- args[2]
# writeMat(output_moutput_mat_filet_file)

data_file <- args[1] #full path to h5 file
cellsubset_file <- args[2] #two column file with full list of cells. Colnames = id, keep
n_pcs <- as.numeric(args[3])
pc_file <- args[4] #output file to write PCs
splitby <- args[5]
saveSO <- args[6]
so_tag <- args[7]

cellsubset <-read.csv(cellsubset_file, row.names = 'id', header = T)

nfeatures <- 3000
ncells <- 5000 #for vst
n_genes <- 2000 #for vst
min_cells <- 5 #for vst
norm_method<-"SCT"
# int_method <- "cca"
int_method <- "rpca"
k.anchor <- 5
refs <- NULL
k.weight <-100
sample.tree <- NULL


data <- Read10X_h5(data_file, use.names=F)

so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellsubset)
so <- subset(so, subset=keep==1)

so.list <- SplitObject(so, split.by=splitby)

# remove(so)
# gc()

so.list <- lapply(X = so.list, FUN = SCTransform, ncells=ncells, variable.features.n=nfeatures,
                   method = "glmGamPoi", n_genes=n_genes, min_cells=min_cells) #vst args
features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
so.list <- lapply(X = so.list, FUN = RunPCA, npcs = n_pcs, verbose=FALSE)
anchors <- FindIntegrationAnchors(
  object.list = so.list, reduction = int_method, dims = 1:n_pcs, k.anchor=k.anchor, 
  anchor.features = features, normalization.method = norm_method, reference = refs) 

# remove(so.list)
# gc()

so <- IntegrateData(anchorset = anchors, normalization.method = norm_method,
                        dims = 1:n_pcs, k.weight=k.weight, sample.tree = sample.tree)

# remove(anchors)
# gc()

DefaultAssay(so) <- "integrated"
#don't do SCTransform again!!!!!
so <- RunPCA(so, verbose=FALSE, npcs = n_pcs)
pcs <- Embeddings(so, reduction = "pca")

write.csv(pcs,pc_file)

#if saveSO, save so as RDS for later loading/subsetting/new PCA.
if (saveSO==1){
  rdsfile<-paste0(paste("integrated",splitby,so_tag,int_method,sep="_"),".rds")
  saveRDS(so, file = file.path(dirname(datapath),rdsfile))
}