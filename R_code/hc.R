#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("sva")
library(sva)
setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/R/")

#load data
all.data <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/bwfiles/SPMR/deeptools/all_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.data)){
  tmp1 <- grep("NaN",all.data[,i])
  a <- c(a,tmp1)
}
all <- all.data[-a,-(1:3)]
colnames(all) <- gsub(".mm10.bw","",colnames(all))
colnames(all) <- gsub("_treat_pileup.bw","",colnames(all))

EPS <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/bwfiles/SPMR/deeptools/eps_bin_100kb.tab",sep = "\t",header = T)
EPS <- EPS[,-(1:3)]
EPS <- EPS[which(rowSums(EPS) > 0),]
colnames(EPS) <- gsub("_treat_pileup.bw","",colnames(EPS))


#hc all
dist_mat <- dist(t(all)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(all)
pdf(file = 'hc_cluster_all.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
dev.off()


test <- all[,c(1:5,8:11,13:22)]
dist_mat <- dist(t(test)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(test)
pdf(file = 'hc_cluster_all.pdf',width = 6,height = 5)
plot(hc_cluster,cex=1,col="dark blue",labels = samplesname)
dev.off()

###combat
batch.info <- data.frame(1:19,row.names = colnames(test))
colnames(batch.info)  <- "samples"
batch.info$batch <- c(rep(1,3),rep(2,8),rep(1,8))
batch.info$type <- c(rep("embryo",3),rep("cellline",8),rep("embryo",8))
mod = model.matrix(~as.factor(type), data=batch.info)
combat_exp <- ComBat(dat = as.matrix(test), batch = batch.info$batch)
hc_cluster_combat <- hclust(dist(t(combat_exp)),method = "complete")
samplesname<- colnames(test)
pdf(file = 'hc_cluster_all_combat.pdf',width = 6,height = 5)
plot(hc_cluster_combat,cex=1,col="dark green",labels = samplesname,main="Cluster after combat")
dev.off()


dist_mat <- dist(t(EPS)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(EPS)
pdf(file = 'hc_cluster_eps.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
dev.off()

tmp <- EPS[,c(1,2,5:8,10,11)]
dist_mat <- dist(t(tmp))#聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(tmp)
pdf(file = 'hc_cluster_eps_3.pdf',width = 6,height = 6)
plot(hc_cluster,cex=1.2,col="dark blue",labels = samplesname)
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

dist_mat <- dist(t(all.LTR)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(all.LTR)
pdf(file = 'hc_cluster_all.LTR.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
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

dist_mat <- dist(t(all.LINE)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(all.LINE)
pdf(file = 'hc_cluster_all.LINE.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
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

dist_mat <- dist(t(all.SINE)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(all.SINE)
pdf(file = 'hc_cluster_all.SINE.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
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

dist_mat <- dist(t(all.repeat)) #聚类分析要先计算距离矩阵
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(all.repeat)
pdf(file = 'hc_cluster_all.repeat.pdf')
plot(hc_cluster,cex=0.8,col="dark blue",labels = samplesname)
dev.off()




