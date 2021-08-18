#script to run SoupX from command line
# - one merged filtered file, raw file per sample
# - cluster ID, sample ID for filtered cells
# - genes specific to each cluster, per sample ("stop" genes for estimating contam)

# stop genes: high fold change relative to all other clusters, but high % in all clusters...

args <- commandArgs(TRUE)

countfile <- args[1] #raw data h5 file
sampleinfofile <- args[2] #csv indicating the sample folders
cellinfofile <- args[3] #csv table with vars: id, sample, cluster
adjustedfile <- args[]

samp.folders <- c("Anterior1","A_P_1","A_P_2","Posterior1")

samp.name <- "merge"
datapath <- paste0("D:/scRNAseq/rat/Feb2020/",samp.name,"/outs/")
data.filt <- Seurat::Read10X_h5(file.path(datapath,"filtered_feature_bc_matrix.h5")) 

clust <- cellinfo$cluster
samp <- cellinfo$sample

samp.names <- unique(cellinfo$sample_orig)

#loop over each separate sample to apply SoupX sample-wise, merge results for final matrix
adj <- list()
sc_samp <- list()
samp_est<-list()
for (i in seq(length(samp.names))) {
  thissub <- samp==samp.names[i]
  filt.sub<-data.filt[,thissub]
  
  #could run SCT/PCA/findClusters here.
  # cluster.id <- cellinfo$leiden_samp[thissub]
  cluster.id <- cellinfo$leiden_fine_samp[thissub]
  # cluster.id <- cellinfo$leiden_cleaned_samp[thissub]
  names(cluster.id) <- rownames(cellinfo[thissub,])
  
  #raw matrix for this sample
  datapath <- paste0("D:/scRNAseq/rat/Feb2020/",samp.folders[i],"/outs/")
  data.raw <- Read10X_h5(file.path(datapath,"raw_feature_bc_matrix.h5")) 
  rownames(data.raw)<-geneinfo$name
  
  sc <- SoupChannel(data.raw, filt.sub)
  sc <- setClusters(sc, cluster.id)
  
  sc = autoEstCont(sc)
  print(sc$fit$rhoFWHM[2])
  sc <- setContaminationFraction(sc, contFrac = sc$fit$rhoFWHM[2])
  thisadj <-adjustCounts(sc, roundToInt=TRUE)
  
  #could also make/save plots here
  
  sc_samp[[i]] <- sc
  samp_est[[i]]<-setNames(sc$soupProfile$est, sc$toc@Dimnames[[1]])
  adj[[i]] <- thisadj
}

data.adj <- do.call(cbind,adj)
write10xCounts(path=file.path(datapath,file), x=data.adj, overwrite = T,
               version="3", gene.id = rownames(geneinfo), gene.symbol = geneinfo$name)

df <- do.call(data.frame,samp_est)
colnames(df) <- samp.names
write.csv(df, 'SoupX_estimates.csv')

saveRDS(sc_samp,file.path(datapath,'SoupChannels_samplewise.rds'))