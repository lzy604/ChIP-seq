setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/R/")
library("ComplexHeatmap")
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")

eps.peak <- read.table('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/bedfiles/eps.H3K9me3.peak.enrichment',
                           sep='\t',header=T,row.names=1)
embyro.peak <- read.table('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/bedfiles/embryo.H3K9me3.peak.enrichment',
                           sep='\t',header=T,row.names=1)
H3K9me3.peak <- rbind(eps.peak,embyro.peak)
H3K9me3.promoter <- H3K9me3.peak[,1]
H3K9me3.LINE <- H3K9me3.peak[,2]
H3K9me3.SINE <- H3K9me3.peak[,3]
H3K9me3.LTR <- H3K9me3.peak[,4]
H3K9me3 <- H3K9me3.peak[,5]

promoter.total<-81348868/2662783660
SINE.total<-206172774/2662783660
LINE.total<-548842036/2662783660
LTR.total<-320032086/2662783660


pdf('peakEnrichment.pdf',width=7.5,height=8)
par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(log((H3K9me3.promoter+1)/H3K9me3/promoter.total,2),
        names=rownames(H3K9me3.peak),
        ylim=c(-2,1),las=2,border=NA,ylab="Observed/Expected Ratio(log2)",
        main="Promoter region",beside=T,col=cccol[1])
barplot(log((H3K9me3.SINE+1)/H3K9me3/SINE.total,2),
        names=rownames(H3K9me3.peak),
        las=2,border=NA,ylab="Observed/Expected Ratio(log2)",
        main="SINE region",beside=T,col=cccol[2])
barplot(log((H3K9me3.LINE+1)/H3K9me3/LINE.total,2),
        names=rownames(H3K9me3.peak),
        las=2,ylim=c(0,1.0),border=NA,ylab="Observed/Expected Ratio(log2)",
        main="LINE region",beside=T,col=cccol[3])
barplot(log((H3K9me3.LTR+1)/H3K9me3/LTR.total,2),
        names=rownames(H3K9me3.peak),
        las=2,ylim=c(0,2),border=NA,ylab="Observed/Expected Ratio(log2)",
        main="LTR region",beside=T,col=cccol[4])
dev.off()

pdf("H3K9me3.region.pdf",width = 5,height = 4)
barplot(H3K9me3[c(1,2,5:8,10,11)]/2662783660*100,
        names=rownames(H3K9me3.peak[c(1,2,5:8,10,11),]),
        las=2,ylim=c(0,4),border=NA,ylab="H3K9me3-marked regions 
        compared to the genome (%)",
        main="",col=cccol[4])
dev.off()


H3K9me3.peak.test <- H3K9me3.peak[c(1,2,5:8,10,11),]
H3K9me3.peak.average <- data.frame()

vals = c(1,3,5,7)
for (i in vals ) {
        tmp <- H3K9me3.peak.test[i:(i+1),]
        tmp2 <- as.data.frame(t(colMeans((tmp))))
        H3K9me3.peak.average <- rbind(H3K9me3.peak.average,tmp2)
}

rownames(H3K9me3.peak.average) <- gsub("..$","",rownames(H3K9me3.peak)[c(1,4,7,10)])

H3K9me3.regio <- data.frame(value=H3K9me3.peak.average[,5]/2662783660*100)
H3K9me3.regio$sample <- row.names(H3K9me3.peak.average)
ggplot(data=H3K9me3.regio, aes(x=sample, y=value)) +
        geom_bar(stat="identity", fill="steelblue")+
        theme_classic()



pdf("H3K9me3.region.pdf",width = 4,height = 5)
par( mar=c(7, 5, 5.1, 5.1), xpd=TRUE)
barplot(H3K9me3.peak.average[,5]/2662783660*100,
        names=rownames(H3K9me3.peak.average),
        las=2,ylim=c(0,4),border=NA,ylab="H3K9me3-marked regions 
        compared to the genome (%)",
        main="",col=cccol[2],font.main = 4)
dev.off()
###number
H3K9me3.peak.number <- read.table('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/bedfiles/eps.H3K9me3.peak.number.enrichment',
                           sep='\t',header=T,row.names=1)

H3K9me3.peak.number.LTR <- H3K9me3.peak.number[,4]

###embryo


