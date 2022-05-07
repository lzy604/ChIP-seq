#read files
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/")
imprinting.gene <- read.csv("imprinting_gene.csv",sep = ",",header = T)
mm10_refGene.bed <-read.csv("mm10_refGene.bed",sep = "\t",header = T)
imprinting.gene.bed <- mm10_refGene.bed[mm10_refGene.bed$symbol %in% imprinting.gene$Gene,]
write.csv(imprinting.gene.bed,"imprinting.gene.bed",quote = F,sep = "\t",row.names = F)

