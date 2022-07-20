#eulerr command line script

#pass in a data table with columns:
# - intersection name
# - intersection count
# - intersection color

# alternative:
# pass in a CSV with named logical columns representing the presence or absence of an item in the set.

# how to handle toggling labels/quantities?

args <- commandArgs(TRUE)

setsfile <- args[1] 
shape <- args[2] 
alpha <- as.numeric(args[3]) 
width <- as.numeric(args[4])
height <- as.numeric(args[5]) 
fsz_counts <- as.numeric(args[6])
fsz_labels <- as.numeric(args[7])
rng_seed <- as.numeric(args[8])
pdf_file <- args[9]

sets_df <-read.csv(setsfile, header = T)

combos <- setNames(unlist(sets_df$counts), unlist(sets_df$names))


control <- list(extraopt=T,
                extraopt_threshold = 0,
                extraopt_control = list(seed = rng_seed, max.time=30,  
                                        verbose=T, smooth=F, temperature=1000))
# list(seed = rng_seed, smooth=F, 
#      maxit=1e100, temperature=10000,verbose=T,
#      nb.stop.improvement=1e3)
print(control)

fit1<-eulerr::euler(combos, shape=shape, control=control)

print(fit1)

fills <- list()
if ("cols" %in% colnames(sets_df) ) {
  fills$fill = unlist(sets_df$cols)
}
# print(fills)

fills$alpha = alpha

quantities <- FALSE
if (fsz_counts>0) {
  quantities <- list(fontsize=fsz_counts)
}

labels <- FALSE
if (fsz_labels>0) {
  labels <- list(fontsize=fsz_labels)
}

pdf_width <- NULL
if (width>0) {
  pdf_width <- width
}
pdf_height <- NULL
if (height>0) {
  pdf_height <- height
}


pdf(pdf_file, width = pdf_width, height = pdf_height)
plot(fit1, 
     quantities = quantities, 
     labels = labels,
     fills = fills)
dev.off()