#script to run DoubletFinder
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(R.matlab))

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

mat <- readMat(parsmatfile)
datafile <- unlist(mat$datafile) #full path to h5 file
cellsubsetfile <- unlist(mat$cellsubsetfile)

dblpars=mat$dblpars[,,1]


cellinfo <-read.csv(cellsubsetfile, row.names = 'id', header = T)

data <- Read10X_h5(datafile)
so <- CreateSeuratObject(data)
so <- AddMetaData(object = so, metadata = cellinfo)
so <- subset(so, keep==1)



# ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
# seu_kidney <- CreateSeuratObject(kidney.data)
# seu_kidney <- NormalizeData(seu_kidney)
# seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
# seu_kidney <- ScaleData(seu_kidney)
# seu_kidney <- RunPCA(seu_kidney)
# seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)
# 
# ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
# seu_kidney <- CreateSeuratObject(kidney.data)
# seu_kidney <- SCTransform(seu_kidney)
# seu_kidney <- RunPCA(seu_kidney)
# seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)
# 
# ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)
# 
# ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)
# 
# ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
# nExp_poi <- round(0.075*nrow(seu_kidney@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
