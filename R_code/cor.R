#install.packages("corrplot")
#install.packages("sva")
#correlation
library(corrplot)
library(RColorBrewer)

#set wd
setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/R/")

#load data
#load data
all.data <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.data)){
  tmp1 <- grep("NaN",all.data[,i])
  a <- c(a,tmp1)
}
all <- all.data[-a,-(1:3)]
colnames(all) <- gsub(".mm10.bw","",colnames(all))
colnames(all) <- gsub("_treat_pileup.bw","",colnames(all))

EPS <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/eps_bin_100kb.tab",sep = "\t",header = T)
EPS <- EPS[,-(1:3)]
EPS <- EPS[which(rowSums(EPS) > 0),]
colnames(EPS) <- gsub("_treat_pileup.bw","",colnames(EPS))

all.promoter <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_promoter_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.promoter)){
  tmp1 <- grep("NaN",all.promoter[,i])
  a <- c(a,tmp1)
}
all.promoter <- all.promoter[-a,-(1:3)]
colnames(all.promoter) <- gsub(".mm10.bw","",colnames(all.promoter))
colnames(all.promoter) <- gsub("_treat_pileup.bw","",colnames(all.promoter))

###cor
cor_all <- cor (as.matrix(all))
pdf("cor_all.pdf")
corrplot.mixed(cor_all,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5
              )
dev.off()


cor_eps <- cor (as.matrix(EPS))
pdf("cor_EPS.pdf")
corrplot.mixed(cor_eps,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()

cor_all.promoter <- cor (as.matrix(all.promoter))
pdf("cor_all.promoter.pdf")
corrplot.mixed(cor_all.promoter,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()

#read data
all.repeat <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_repeat_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.repeat)){
  tmp1 <- grep("NaN",all.repeat[,i])
  a <- c(a,tmp1)
}
all.repeat <- all.repeat[-a,-(1:3)]
colnames(all.repeat) <- gsub(".mm10.bw","",colnames(all.repeat))
colnames(all.repeat) <- gsub("_treat_pileup.bw","",colnames(all.repeat))

cor_all.repeat <- cor (as.matrix(all.repeat))
pdf("cor_all.repeat.pdf")
corrplot.mixed(cor_all.repeat,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()

#read data
all.LINE <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_LINE_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.LINE)){
  tmp1 <- grep("NaN",all.LINE[,i])
  a <- c(a,tmp1)
}
all.LINE <- all.LINE[-a,-(1:3)]
colnames(all.LINE) <- gsub(".mm10.bw","",colnames(all.LINE))
colnames(all.LINE) <- gsub("_treat_pileup.bw","",colnames(all.LINE))

cor_all.LINE <- cor (as.matrix(all.LINE))
pdf("cor_all.LINE.pdf")
corrplot.mixed(cor_all.LINE,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()


#read data
all.SINE <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_SINE_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.SINE)){
  tmp1 <- grep("NaN",all.SINE[,i])
  a <- c(a,tmp1)
}
all.SINE <- all.SINE[-a,-(1:3)]
colnames(all.SINE) <- gsub(".mm10.bw","",colnames(all.SINE))
colnames(all.SINE) <- gsub("_treat_pileup.bw","",colnames(all.SINE))

cor_all.SINE <- cor (as.matrix(all.SINE))
pdf("cor_all.SINE.pdf")
corrplot.mixed(cor_all.SINE,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()

#read data
all.LTR <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_LTR_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.LTR)){
  tmp1 <- grep("NaN",all.LTR[,i])
  a <- c(a,tmp1)
}
all.LTR <- all.LTR[-a,-(1:3)]
colnames(all.LTR) <- gsub(".mm10.bw","",colnames(all.LTR))
colnames(all.LTR) <- gsub("_treat_pileup.bw","",colnames(all.LTR))

cor_all.LTR <- cor (as.matrix(all.LTR))
pdf("cor_all.LTR.pdf")
corrplot.mixed(cor_all.LTR,tl.pos="lt",cl.lim = c(0,1),order="hclust",
               tl.col = "black",
               number.cex=0.5)
dev.off()







colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
pheatmap(cor_all)

