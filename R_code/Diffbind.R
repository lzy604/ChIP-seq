# BiocManager::install("DiffBind")
library(DiffBind)
#read sample
sampleSheet <- read.csv("SampleSheet.csv")
dbObj <- dba(sampleSheet=sampleSheet)
pdf("all.cor.heatmap.pdf")
plot(dbObj)
dev.off()
pdf("all.plotPCA.pdf")
dba.plotPCA(dbObj,  attributes=DBA_TISSUE, label=DBA_ID)
dev.off()

#EPS2c_vs_ESC
EPS2c_vs_ESC <- sampleSheet[c(grep("EPSC-2C*|ESC*",sampleSheet$SampleID )),]
dbObj <- dba(sampleSheet=EPS2c_vs_ESC)
pdf("EPS2c_vs_ESC.cor.heatmap.pdf")
plot(dbObj)
dev.off()
#count
dbObj <- dba.count(dbObj)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,bRemoveDuplicates=TRUE)
#contrast and analyze
dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("EPS2c_vs_ESC.Venn.pdf")
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()
# EdgeR
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="EPS2c_vs_ESC_edgeR.txt", sep="\t", quote=F, col.names = NA)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="EPS2c_vs_ESC_edgeR_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="EPS2c_vs_ESC_deseq2.txt", sep="\t", quote=F, col.names = NA)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="EPS2c_vs_ESC_deseq2_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)

#EPS8c_vs_ESC
EPS8c_vs_ESC <- sampleSheet[c(grep("EPSC-8C*|ESC*",sampleSheet$SampleID )),]
dbObj <- dba(sampleSheet=EPS8c_vs_ESC)
pdf("EPS8c_vs_ESC.cor.heatmap.pdf")
plot(dbObj)
dev.off()
#count
dbObj <- dba.count(dbObj)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,bRemoveDuplicates=TRUE)
#contrast and analyze
dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("EPS8c_vs_ESC.Venn.pdf")
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()
# EdgeR
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="EPS8c_vs_ESC_edgeR.txt", sep="\t", quote=F, col.names = NA)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="EPS8c_vs_ESC_edgeR_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="EPS8c_vs_ESC_deseq2.txt", sep="\t", quote=F, col.names = NA)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="EPS8c_vs_ESC_deseq2_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)

#TdEPS_vs_ESC
TdEPS_vs_ESC <- sampleSheet[c(grep("TdEPS*|ESC*",sampleSheet$SampleID )),]
dbObj <- dba(sampleSheet=TdEPS_vs_ESC)
pdf("TdEPS_vs_ESC.cor.heatmap.pdf")
plot(dbObj)
dev.off()
#count
dbObj <- dba.count(dbObj)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,bRemoveDuplicates=TRUE)
#contrast and analyze
dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("TdEPS_vs_ESC.Venn.pdf")
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()
# EdgeR
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="TdEPS_vs_ESC_edgeR.txt", sep="\t", quote=F, col.names = NA)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="TdEPS_vs_ESC_edgeR_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="TdEPS_vs_ESC_deseq2.txt", sep="\t", quote=F, col.names = NA)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="TdEPS_vs_ESC_deseq2_sig0.05.bed", sep="\t", quote=F, row.names=F, col.names=F)
