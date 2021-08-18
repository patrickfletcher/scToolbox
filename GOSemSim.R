library(GOSemSim)

args <- commandArgs(TRUE)

cat("# args passed at the command line\n")
cat(length(args));cat("\n")

file <- args[1]
ont <- args[2]
sim_method <- args[3]
outfile <- args[4]

go_sub <- read.table(file, sep=',', header = F)
go_terms <- setNames(go_sub,go_sub)

if (sim_method=="Wang"){
  doIC <- FALSE
} else {
  doIC <- TRUE
}

hsGO <- GOSemSim::godata('org.Rn.eg.db', ont=ont, computeIC=doIC)
simMatrix <- mgoSim(go_terms, go_terms, measure=sim_method, semData = hsGO, combine = NULL)

write.table(simMatrix, file=outfile,col.names = FALSE)