setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/R/")
#install.packages("devtools", repo="http://cran.us.r-project.org")
#library(devtools)
#install.packages("ggplot2")
#install_github("vqv/ggbiplot")
library(ggbiplot)

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
EPS <- EPS[,c(1,2,5:8,10,11)]
EPS <- EPS[which(rowSums(EPS) > 0),]

all.promoter <- read.table("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bwfiles/SPMR/all_promoter_bin_100kb.tab",sep = "\t",header = T)
a <- vector()
for (i in 1: ncol(all.promoter)){
  tmp1 <- grep("NaN",all.promoter[,i])
  a <- c(a,tmp1)
}
all.promoter <- all.promoter[-a,-(1:3)]
colnames(all.promoter) <- gsub(".mm10.bw","",colnames(all.promoter))
colnames(all.promoter) <- gsub("_treat_pileup.bw","",colnames(all.promoter))



#pca all
pca.all <- prcomp(t(all), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all$x)
pc_df$Condition <- rownames(pc_df)
pdf("pca_all.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
dev.off()

#pca eps
pca.EPS <- prcomp(t(EPS), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.EPS$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.EPS$x)
pc_df$Condition <- gsub("..$","",rownames(pc_df)) 

pdf("pca_EPS.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
dev.off()

#pca all promter
pca.all.promoter <- prcomp(t(all.promoter), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all.promoter$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all.promoter$x)
pc_df$Condition <- rownames(pc_df)

pdf("pca.all.promoter.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
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


#pca all repeat
pca.all.repeat <- prcomp(t(all.repeat), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all.repeat$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all.repeat$x)
pc_df$Condition <- rownames(pc_df)

pdf("pca.all.repeat.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
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


#pca all LINE
pca.all.LINE <- prcomp(t(all.LINE), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all.LINE$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all.LINE$x)
pc_df$Condition <- rownames(pc_df)

pdf("pca.all.LINE.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
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


#pca all SINE
pca.all.SINE <- prcomp(t(all.SINE), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all.SINE$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all.SINE$x)
pc_df$Condition <- rownames(pc_df)

pdf("pca.all.SINE.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
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


#pca all LTR
pca.all.LTR <- prcomp(t(all.LTR), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca.var <- pca.all.LTR$sdev^2 # 计算原始数据中的每个数据在每个PC上的比重
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #计算每个PC占所有PC的和的比列
pc_df <- as.data.frame(pca.all.LTR$x)
pc_df$Condition <- rownames(pc_df)

pdf("pca.all.LTR.pdf")
ggplot(pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+  
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = ""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = ""))
dev.off()
