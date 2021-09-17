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
alpha <- args[3] 
width <- args[4] 
height <- args[5] 
fsz_counts <- args[6]
fsz_labels <- args[7]
pdf_file <- args[8]

sets_df <-read.csv(setsfile, header = T)

combos <- setNames(unlist(sets_df$counts), unlist(sets_df$names))

fit1<-eulerr::euler(combos, shape=shape)

fills <- list()
if ("cols" %in% colnames(sets_df) ) {
  fills$fill = unlist(sets_df$cols)
}
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