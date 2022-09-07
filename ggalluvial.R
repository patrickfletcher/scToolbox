library(ggplot2)
library(ggalluvial)
#stratum: rectangles, x: xgroupings, 
#alluvia: splines connecting 

# args <- commandArgs(TRUE)
# parsmatfile <- args[1]

# mat <- readMat(parsmatfile)
# 
# datafile <- unlist(mat$datafile) 
# savefile <- unlist(mat$savefile) 
# 
# xvar <- unlist(mat$xvar) 
# yvar <- unlist(mat$yvar)
# strat <- unlist(mat$strat)
# alluv <- unlist(mat$alluv)
# fillvar <- unlist(mat$fillvar)    

setwd("C:/Users/fletcherpa/Box/scRNAseq/rat/pineal/matlab")
datafile <- "cleaned_pinast_age_by_clust.csv"


# datafile='D:/tmp/tmp_ggalluvial.csv';

data <-read.csv(datafile, header = T)

clustnames<-as.character(sort(unique(data$Cluster)))

colsfile <- "cleaned_pinast_age_by_clust_cols.csv"
cols <-read.csv(colsfile, header = F)
cols <-setNames(unlist(cols),clustnames)


clust <- factor(data$Cluster,levels = clustnames)
age <- factor(data$Age,levels = c("E19","E21","P5","P20","P40"))
tab <- table(Cluster=clust, Age=age)

df<-as.data.frame(tab)
df2<-as.data.frame(proportions(tab,2))
df$Proportion<-df2$Freq*100

x<-c(rep(-3,11),rep(-1,11),rep(5,11),rep(20,11),rep(40,11))

df$x=x
# alluvial plot
ggplot(df,
       aes(x=x, 
           y=Proportion,
           stratum = Cluster, 
           alluvium = Cluster,
           fill = Cluster, 
           label = Cluster)) +
  geom_flow(stat = "alluvium", 
            color = "gray32",
            alpha = .5) +
  geom_stratum(stat="stratum",
               alpha = 1) +
  scale_x_continuous(breaks=c(0,20,40),
                     expand = c(.05, .05)) + 
  scale_fill_manual(values=cols) & theme_minimal()

# lode.guidance = "frontback",
# 
# scale_x_discrete(expand = c(.05, .05)) + 
