library(ggplot2)
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/")

#read files
tmp <- read.csv("caculate.sh.out",header = F)
diff.result.eps2c.esc <- data.frame(row.names = c("Total","LTR","SINE","LINE","imprint_gene"))
rownames(diff.result.eps2c.esc) <-c("Total","LTR","SINE","LINE","imprint_gene")
for (i in seq(1, nrow(tmp), 5)) {
  tmp1 <- data.frame(tmp[(i:(i+4)),])
  colnames(tmp1) <- tmp[i,]
  diff.result.eps2c.esc <- cbind(diff.result.eps2c.esc,tmp1)
}

diff.result.eps2c.esc <-as.data.frame(t(diff.result.eps2c.esc))
diff.result.eps2c.esc$Total <- gsub(' .*','',row.names(diff.result.eps2c.esc))
diff.result.eps2c.esc <-as.data.frame(apply(diff.result.eps2c.esc, 2, as.numeric))
row.names(diff.result.eps2c.esc) <- c("EPS2c_vs_ESC_deseq2_sig0.05_nfc.bed",
                                      "EPS2c_vs_ESC_deseq2_sig0.05_pfc.bed",
                                      "EPS8c_vs_ESC_deseq2_sig0.05_nfc.bed",
                                      "EPS8c_vs_ESC_deseq2_sig0.05_pfc.bed",
                                      "TdEPS_vs_ESC_deseq2_sig0.05_nfc.bed",
                                      "TdEPS_vs_ESC_deseq2_sig0.05_pfc.bed")
diff.result.eps2c.esc.ratio <- diff.result.eps2c.esc[,2:5]/diff.result.eps2c.esc$Total
diff.result.pfc.ratio <- diff.result.eps2c.esc.ratio[grep('pfc*',row.names(diff.result.eps2c.esc)),]
diff.result.nfc.ratio <- diff.result.eps2c.esc.ratio[grep('nfc*',row.names(diff.result.eps2c.esc)),]

#plot
pfc.barplot <- data.frame(sample=rep(gsub("_deseq2_sig0.05_pfc.bed","",row.names(diff.result.pfc.ratio)), each=ncol(diff.result.pfc.ratio)),
                          region=rep(colnames(diff.result.pfc.ratio),nrow(diff.result.pfc.ratio)),
                          value =as.vector(t(diff.result.pfc.ratio)))
pdf("pfc.barplot.pdf")
ggplot(data=pfc.barplot, aes(x=sample, y=value,fill=region)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()
dev.off()

nfc.barplot <- data.frame(sample=rep(gsub("_deseq2_sig0.05_nfc.bed","",row.names(diff.result.nfc.ratio)), each=ncol(diff.result.nfc.ratio)),
                          region=rep(colnames(diff.result.nfc.ratio),nrow(diff.result.nfc.ratio)),
                          value =as.vector(t(diff.result.nfc.ratio)))
pdf("nfc.barplot.pdf")
ggplot(data=nfc.barplot, aes(x=sample, y=value,fill=region)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()
dev.off()
##Barplot with error

