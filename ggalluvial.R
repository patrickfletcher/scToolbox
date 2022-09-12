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

clustnames<-as.character(unique(data$Cluster))
# clustnames<-clustnames[seq(2,11)]
clustnames<-clustnames[clustnames!="<undefined>"]
K=length(clustnames)

colsfile <- "cleaned_pinast_age_by_clust_cols.csv"
coltab <-read.csv(colsfile, header = T)
clustnames<-coltab$clust
cols <- coltab$hexcol
cols <-setNames(unlist(cols),clustnames)


ages<-c(-3,-1,5,20,40)
agelabs<-c("E19","E21","P5","P20","P40")

clust <- factor(data$Cluster,levels = clustnames)
age <- factor(data$Age,levels = agelabs)
tab <- table(Cluster=clust, Age=age)

df<-as.data.frame(tab)
df2<-as.data.frame(proportions(tab,2))
df$Proportion<-df2$Freq*100

age<-c(rep(-3,K),rep(-1,K),rep(5,K),rep(20,K),rep(40,K))
# age<-c(rep(0,K),rep(5,K),rep(10,K),rep(15,K),rep(20,K))

df$age=age

# alluvial plot - equal spacing between ages
plotfile<-"./figs/Fig1_cleaned_clust_by_age.pdf"
pdf(plotfile, width=3, height=2, pointsize = 8)

plotfile<-"./figs/Fig1_cleaned_clust_by_age.png"
png(plotfile, units ="in", width=3, height=2,  pointsize = 8, res=300,
    type="cairo", antialias = "subpixel")
ggplot(df,
       aes(x=Age, 
           y=Proportion,
           stratum = Cluster, 
           alluvium = Cluster,
           fill = Cluster, 
           label = Cluster)) + 
  geom_flow(stat = "alluvium", color = "gray32", knot.pos = 1/10, size=0.2,
            alpha = .6, width=.1,  show.legend = F) +
  geom_stratum(stat="stratum", size=0.2,
               alpha = 1, width=.1, show.legend = F) + 
  scale_x_discrete(expand = c(.015, .015)) +
  scale_y_continuous(expand = c(.01, .01)) + ylab("% Cells Per Age") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_line(size = 0.2)) + 
  scale_fill_manual(values=cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()


# alluvial plot - actual spacing between ages
plotfile<-"./figs/Fig1_cleaned_clust_by_age2.pdf"
pdf(plotfile, width=3, height=2, pointsize = 8)

plotfile<-"./figs/Fig1_cleaned_clust_by_age2.png"
png(plotfile, units ="in", width=3, height=2,  pointsize = 8, res=300,
    type="cairo", antialias = "subpixel")
ggplot(df,
       aes(x=age, 
           y=Proportion,
           stratum = Cluster, 
           alluvium = Cluster,
           fill = Cluster, 
           label = Cluster)) + 
  geom_flow(stat = "alluvium", color = "gray32", knot.pos = 1/10, size=0.2,
            alpha = .6, width=0.75,  show.legend = F) +
  geom_stratum(stat="stratum", size=0.2,
               alpha = 1, width=0.75, show.legend = F) + 
  scale_x_continuous(expand = c(.005, .005), breaks = ages, labels = agelabs) + xlab("Age") + 
  scale_y_continuous(expand = c(.01, .01)) + ylab("% Cells Per Age") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_line(size = 0.2)) + 
  scale_fill_manual(values=cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

