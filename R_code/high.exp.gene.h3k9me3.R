#read file
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/count/")
eps2c.esc.diff.exp <- read.csv("EPS2C.vs.esc.results.csv")
eps8c.esc.diff.exp <- read.csv("eps8c.vs.esc.results.csv")
TdEPS.esc.diff.exp <- read.csv("TdEPS.vs.esc.results.csv")
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/bedfiles/")
EPS2c.h3k9me3.promoter <- read.table("EPSC-2C-5_overlaps.promoter.bed",sep = "\t")
ESC.h3k9me3.promoter <- read.table("ESC_overlaps.promoter.bed",sep = "\t")
EPS8c.h3k9me3.promoter <- read.table("EPSC-8C-4_overlaps.promoter.bed",sep = "\t")
TdEPS.h3k9me3.promoter <- read.table("TdEPS_overlaps.promoter.bed",sep = "\t")
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/bedfiles/")
EPS2c.h3k9me3.gene.promoter <- read.table("EPSC-2C-5_overlaps.gene.promoter.bed",sep = "\t")
ESC.h3k9me3.gene.promoter <- read.table("ESC_overlaps.gene.promoter.bed",sep = "\t")
EPS8c.h3k9me3.gene.promoter <- read.table("EPSC-8C-4_overlaps.gene.promoter.bed",sep = "\t")
TdEPS.h3k9me3.gene.promoter <- read.table("TdEPS_overlaps.gene.promoter.bed",sep = "\t")
####EPS2C highly expressed gene
eps2c.esc.diff.exp$EPSC.2C.mean <- (eps2c.esc.diff.exp$EPSC.2C.129+eps2c.esc.diff.exp$EPSC.2C.5)/2
EPS2c.high.no.h3k9.gene <- setdiff(eps2c.esc.diff.exp$Row.names[eps2c.esc.diff.exp$log2FoldChange>1&eps2c.esc.diff.exp$EPSC.2C.mean>1],
                                   unique(EPS2c.h3k9me3.gene.promoter$V14))
EPS2c.high.ESC.h3k9.gene<- intersect(EPS2c.high.no.h3k9.gene,ESC.h3k9me3.gene.promoter$V14)
EPS2c.high.ESC.h3k9.gene.5fc <- intersect(EPS2c.high.ESC.h3k9.gene,ESC.h3k9me3.gene.promoter$V14[ESC.h3k9me3.gene.promoter$V7>= 5])

####EPS8c highly expressed gene
EPS8c.high.no.h3k9.gene <- setdiff(eps8c.esc.diff.exp$Row.names[eps8c.esc.diff.exp$log2FoldChange>1&eps8c.esc.diff.exp$EPSC.8C.4>1],
                                   unique(EPS8c.h3k9me3.gene.promoter$V14))
EPS8c.high.ESC.h3k9.gene<- intersect(EPS8c.high.no.h3k9.gene,ESC.h3k9me3.gene.promoter$V14)
EPS8c.high.ESC.h3k9.gene.5fc <- intersect(EPS8c.high.ESC.h3k9.gene,ESC.h3k9me3.gene.promoter$V14[ESC.h3k9me3.gene.promoter$V7>= 5])

####TdEPS highly expressed gene
TdEPS.high.no.h3k9.gene <- setdiff(TdEPS.esc.diff.exp$Row.names[TdEPS.esc.diff.exp$log2FoldChange>1&TdEPS.esc.diff.exp$C1.18.TdEPS>1],
                                   unique(TdEPS.h3k9me3.gene.promoter$V14))
TdEPS.high.ESC.h3k9.gene<- intersect(TdEPS.high.no.h3k9.gene,ESC.h3k9me3.gene.promoter$V14)
TdEPS.high.ESC.h3k9.gene.5fc <- intersect(TdEPS.high.ESC.h3k9.gene,ESC.h3k9me3.gene.promoter$V14[ESC.h3k9me3.gene.promoter$V7>= 5])

####EPS2C highly expressed gene
eps2c.esc.diff.exp$EPSC.2C.mean <- (eps2c.esc.diff.exp$EPSC.2C.129+eps2c.esc.diff.exp$EPSC.2C.5)/2
EPS2c.high.no.h3k9.gene <- setdiff(eps2c.esc.diff.exp$Row.names[eps2c.esc.diff.exp$log2FoldChange>1&eps2c.esc.diff.exp$EPSC.2C.mean>1],
                                   unique(EPS2c.h3k9me3.promoter$V14))
EPS2c.high.ESC.h3k9.gene<- intersect(EPS2c.high.no.h3k9.gene,ESC.h3k9me3.promoter$V14)
EPS2c.high.ESC.h3k9.gene.5fc <- intersect(EPS2c.high.ESC.h3k9.gene,ESC.h3k9me3.promoter$V14[ESC.h3k9me3.promoter$V7>= 5])

####EPS8c highly expressed gene
EPS8c.high.no.h3k9.gene <- setdiff(eps8c.esc.diff.exp$Row.names[eps8c.esc.diff.exp$log2FoldChange>1&eps8c.esc.diff.exp$EPSC.8C.4>1],
                                   unique(EPS8c.h3k9me3.promoter$V14))
EPS8c.high.ESC.h3k9.gene<- intersect(EPS8c.high.no.h3k9.gene,ESC.h3k9me3.promoter$V14)
EPS8c.high.ESC.h3k9.gene.5fc <- intersect(EPS8c.high.ESC.h3k9.gene,ESC.h3k9me3.promoter$V14[ESC.h3k9me3.promoter$V7>= 5])

####TdEPS highly expressed gene
TdEPS.high.no.h3k9.gene <- setdiff(TdEPS.esc.diff.exp$Row.names[TdEPS.esc.diff.exp$log2FoldChange>1&TdEPS.esc.diff.exp$C1.18.TdEPS>1],
                                   unique(TdEPS.h3k9me3.promoter$V14))
TdEPS.high.ESC.h3k9.gene<- intersect(TdEPS.high.no.h3k9.gene,ESC.h3k9me3.promoter$V14)
TdEPS.high.ESC.h3k9.gene.5fc <- intersect(TdEPS.high.ESC.h3k9.gene,ESC.h3k9me3.promoter$V14[ESC.h3k9me3.promoter$V7>= 5])

