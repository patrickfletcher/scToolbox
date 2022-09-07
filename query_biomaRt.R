library(biomaRt)

args <- commandArgs(TRUE)

cat("# args passed at the command line\n")
cat(length(args));cat("\n")

#inputs: 1) dataset, 2) attributes, 3) filters, 4) values
# 2, 3, 4 must be external files (pass file names). Finally optionally specify an output file

dataset <- args[1]
attributes_file <- args[2]
filter_name <- args[3]
vals_file <- args[4] #should contain a single column of values to filter on
outfile <- args[5]

df <- read.csv(attributes_file, header = F)
attributes <- as.character(df[,1]) #make a vector from the dataframe
print(length(attributes))

df <- read.csv(vals_file, header = F)
vals <- as.character(df[,1]) #make a vector from the dataframe
print(length(vals))

df=listEnsemblArchives()
version="104"
host_url=df$url[df$version==version]

ensembl = useMart("ensembl",dataset=dataset, host=host_url)

gene.data <- getBM(attributes=c('ensembl_gene_id','external_gene_name', attributes),
                   filters = filter_name, values = vals, mart = ensembl)


write.table(gene.data, file=outfile, col.names = TRUE, row.names = FALSE)