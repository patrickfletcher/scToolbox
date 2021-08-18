#script to run SoupX's quickMarkers from command line

#need to get count matrix, sparse
#subset to desired cells -> cluster as csv with cellIDs as one column
#options as args
#write resulting dataframe back to a csv

args <- commandArgs(TRUE)

countfile <- args[1]
clustfile <- args[2]
N <- args[3]
FDR <- args[4]
expressCut <- args[5]
resultfile <- args[6]

counts <- Seurat::Read10X_h5(countfile)
clusters <- read.csv(clustfile, header=T, row.names = 'id', colClasses = "character")

# print(class(counts))
# print(dim(counts))
# print(dim(clusters))

counts <- counts[,rownames(clusters)]

# print(dim(counts))
# print(head(clusters))
# print(class(clusters[1]))

clusters <- setNames(as.character(clusters$cluster),rownames(clusters))

result <- SoupX::quickMarkers(counts, clusters, N = N, FDR = FDR, expressCut = expressCut)

write.csv(result, file=resultfile)