#BiocManager::install("ComplexHeatmap")
#install.packages("DescTools")
library("ComplexHeatmap")
library(DescTools)
setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/heatmap/results/")
H3K9me3<-read.table('H3K9me_cluster_signal.txt',sep = "\t",header = T,row.names = 1)
H3K9me3.filter <- data.frame()
for (i in 1:nrow(H3K9me3)) {
  tmp <- H3K9me3[i,]
  if (any(c(tmp[, grep(".[1-3]$",colnames(H3K9me3))]) >=5)) {
    H3K9me3.filter <- rbind(H3K9me3.filter,tmp)
  }
}

H3K9me3$C1.18.TdEPS.value <- ((H3K9me3$C1.18.TdEPS.1+H3K9me3$C1.18.TdEPS.2)/2)/(H3K9me3$C1.18.TdEPS.input+1)
H3K9me3$EPSC.2C.5.value <- ((H3K9me3$EPSC.2C.5.2+H3K9me3$EPSC.2C.5.3)/2)/(H3K9me3$EPSC.2C.5.input+1)
H3K9me3$EPSC.8C.4.value <- ((H3K9me3$EPSC.8C.4.1+H3K9me3$EPSC.8C.4.2)/2)/(H3K9me3$EPSC.8C.4.input+1)
H3K9me3$ESC.value <- ((H3K9me3$ESC.1+H3K9me3$ESC.2)/2)/(H3K9me3$ESC.input+1)

H3K9me3$H <- apply(H3K9me3[,grep("*.value",colnames(H3K9me3))],1,Entropy)

write.table(H3K9me3[,18:22],"EPS.1kb.k9.vsinput.csv",sep = "\t",row.names = T,col.names = T,quote = F)
set.seed(100)
pm1<-as.matrix(t(scale(t(log2(H3K9me3[,grep("*.value",colnames(H3K9me3))]+1)))))
write.table(pm1,"EPS.1kb.log2.csv",sep = "\t",row.names = T,col.names = T,quote = F)
colnames(pm1) <- gsub("*.value","",colnames(pm1))

pdf("eps.h3k9me3.k7.pdf", useDingbats=FALSE)
Heatmap(as.matrix(pm1[,c(4,1,3,2)]),km=7,name = "Normalized H3K9me3/input ratio",
              use_raster=FALSE,cluster_columns = FALSE,
              show_row_names = F) 
dev.off()

set.seed(100)
pdf("eps.h3k9me3.k5.pdf")
Heatmap(as.matrix(pm1[,c(4,1,3,2)]),km=5,name = "Normalized H3K9me3/input ratio",
              use_raster=FALSE,cluster_columns = FALSE,
              show_row_names = F) 
dev.off()

set.seed(100)
pm2<-as.matrix(t(scale(t(log2(H3K9me3[H3K9me3$H<1.5,grep("*.value",colnames(H3K9me3))]+1)))))
colnames(pm2) <- gsub("*.value","",colnames(pm2))
pdf("eps.h3k9me3.k5.h1.5.pdf")
ht1 = Heatmap(pm2[,c(4,1:3)],km=4,name = "Normalized H3K9me3/input ratio",
              use_raster=FALSE,cluster_columns = FALSE,
              show_row_names = F) 
ht1 = draw(ht1)
dev.off()

set.seed(100)
pm2<-as.matrix(t(scale(t(log2(H3K9me3[H3K9me3$H<1.7,grep("*.value",colnames(H3K9me3))]+1)))))
colnames(pm2) <- gsub("*.value","",colnames(pm2))
pdf("eps.h3k9me3.k5.h1.7.pdf")
ht1 = Heatmap(pm2[,c(4,1:3)],km=6,name = "Normalized H3K9me3/input ratio",
              use_raster=FALSE,cluster_columns = FALSE,
              show_row_names = F) 
ht1 = draw(ht1)
dev.off()

set.seed(100)
pm3<-as.matrix(t(scale(t(log2(H3K9me3[H3K9me3$H<1.95,grep("*.value",colnames(H3K9me3))]+1)))))
colnames(pm3) <- gsub("*.value","",colnames(pm3))
pdf("eps.h3k9me3.k5.h1.9.pdf")
ht1 = Heatmap(as.matrix(pm3[,c(4,1:3)]),km=5,name = "Normalized H3K9me3/input ratio",
              use_raster=FALSE,cluster_columns = FALSE,
              show_row_names = F) 
ht1 = draw(ht1)
dev.off()
##output cluster
r_o <- row_order(ht1)
pm3 <-as.data.frame(pm3)
pm3$cl=100
pm3[r_o$`1`,]$cl =1
pm3[r_o$`2`,]$cl =2 
pm3[r_o$`3`,]$cl =3
pm3[r_o$`4`,]$cl =4 
pm3[r_o$`5`,]$cl =5
# to bed
library(stringr)
bed.2c5 <- as.data.frame(rownames(pm3)[pm3$cl == "1"])
bed.2c5<- as.data.frame(str_split_fixed(bed.2c5$`rownames(pm3)[pm3$cl == "1"]`, "_", 2))
bed.2c5$V2 <-as.numeric(as.character(bed.2c5$V2))
bed.2c5$V3 <- (bed.2c5$V2)+1000
write.table(bed.2c5,"bed.2c5",sep = "\t",col.names = F,row.names = F,quote = F)
bed.8c4<- as.data.frame(str_split_fixed(rownames(pm3)[pm3$cl == "2"], "_", 2))
bed.8c4$V2 <-as.numeric(as.character(bed.8c4$V2))
bed.8c4$V3 <- (bed.8c4$V2)+1000
write.table(bed.8c4,"bed.8c4",sep = "\t",col.names = F,row.names = F,quote = F)

bed.tdeps<- as.data.frame(str_split_fixed(rownames(pm3)[pm3$cl == "3"], "_", 2))
bed.tdeps$V2 <-as.numeric(as.character(bed.tdeps$V2))
bed.tdeps$V3 <- (bed.tdeps$V2)+1000
write.table(bed.tdeps,"bed.tdeps",sep = "\t",col.names = F,row.names = F,quote = F)

bed.esc<- as.data.frame(str_split_fixed(rownames(pm3)[pm3$cl == "5"], "_", 2))
bed.esc$V2 <-as.numeric(as.character(bed.esc$V2))
bed.esc$V3 <- (bed.esc$V2)+1000
write.table(bed.esc,"bed.esc",sep = "\t",col.names = F,row.names = F,quote = F)


###embryo

embryo.H3K9me3 <- read.csv("embryo.H3k9.mm10.1kb.txt",sep = "\t",header = T)
rownames(embryo.H3K9me3) <- paste(embryo.H3K9me3$X..chr.,embryo.H3K9me3$X.start.,sep = "_")
embryo.H3K9me3$Oocyte.value <- embryo.H3K9me3$X.Oocyte_M2./(embryo.H3K9me3$X.Oocyte_M2_ctr.+1)
embryo.H3K9me3$Zygote.value <- embryo.H3K9me3$X.Zygote./(embryo.H3K9me3$X.Zygote_ctr.+1)
embryo.H3K9me3$E2cell.value <- embryo.H3K9me3$X.2cell./(embryo.H3K9me3$X.2cell_ctr.+1)
embryo.H3K9me3$E4cell.value <- embryo.H3K9me3$X.4cell./(embryo.H3K9me3$X.4cell_ctr.+1)
embryo.H3K9me3$E8cell.value <- embryo.H3K9me3$X.8cell./(embryo.H3K9me3$X.8cell_ctr.+1)
embryo.H3K9me3$Morula.value <- embryo.H3K9me3$X.Morula./(embryo.H3K9me3$X.Morula_ctr.+1)
embryo.H3K9me3$ICM.value <- embryo.H3K9me3$X.ICM./(embryo.H3K9me3$X.ICM_ctr.+1)


  
  

col_fpkm1 <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
col_fpkm2 <- colorRampPalette(c("#4292C6", "white", "#EF3B2C"))(50)
col_fpkm3 <- colorRampPalette(c("#8073ac", "white", "#e08214"))(50)
col_fpkm3 <- colorRampPalette(c("#1f78b4","#a6cee3","#6a3d9a","#cab2d6","#f5f5f5",
                                "#fdbf6f","#ff7f00","#fb9a99", "#e31a1c"))(100)


set.seed(100)
pdf("embryo.vs.EPS.K9.heatmap.pdf")
pm.eps <- as.matrix(t(scale(t(log2(H3K9me3[H3K9me3$H<1.9,grep("*.value",colnames(H3K9me3))]+1)))))
ht3 <- Heatmap(pm3,name="EPS",k=5,
               col = col_fpkm1,use_raster=FALSE,cluster_rows=T,
               cluster_columns = FALSE,show_row_names = F)
#ht1 = draw(ht3)
tmp <- embryo.H3K9me3[rownames(pm3),
                      grep(".value",colnames(embryo.H3K9me3))]

pm.embryo <-as.matrix(t(scale(t(log2(tmp+1)))))
ht4 <- Heatmap(pm.embryo,name="Embryo",
               col = col_fpkm2,use_raster=FALSE,
               cluster_columns = FALSE,show_row_names = F)
ht3+ht4
dev.off()
##eps8c.esc
set.seed(100)
pdf("embryo.vs.EPS8c.K9.heatmap.pdf")
eps.8c <- H3K9me3[,c(20,21)]
eps.8c <- eps.8c[-(which(rowSums(eps.8c) == 0)),]
eps.8c$fc <- log2((eps.8c$EPSC.8C.4.value+1)/(eps.8c$ESC.value+1))


tmp <- eps.8c[abs(eps.8c$fc)>1,]
pm.eps8c <- as.matrix(log2(tmp[,c(1,2)]+1))
pm.eps8c <- as.data.frame(pm.eps8c)
colnames(pm.eps8c) <- gsub(".value","",colnames(pm.eps8c))
library(reshape2)
tmp2 <- melt(pm.eps8c)
library(ggplot2)
ggplot(tmp2, aes(x=value, color=Var2)) +
  geom_density()
library(dplyr)

pdf("pm.eps8c.logfc1.pdf")
ht3 <- Heatmap(as.matrix(pm.eps8c),name="EPS8c",k=2,
               col = col_fpkm3,use_raster=FALSE,cluster_rows=T,
               cluster_columns = FALSE,show_row_names = F)
draw(ht3)
dev.off()

tmp <- eps.8c[abs(eps.8c$fc)>0.585,]
pm.eps8c <- as.matrix(log2(tmp[,c(1,2)]+1))
pdf("pm.eps8c.logfc0.585.pdf")
ht3 <- Heatmap(as.matrix(pm.eps8c),name="EPS8c",k=2,
               col = col_fpkm3,use_raster=FALSE,cluster_rows=T,
               cluster_columns = FALSE,show_row_names = F)
draw(ht3)
dev.off()


