library(reshape)
library(ggplot2)
setwd("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/R/")
mapping <- read.csv("/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/align/mapping.csv")
mapping$ratio <- mapping$UniquelyMapped/mapping$Total
mapping$Unmapped <- mapping$Total-mapping$UniquelyMapped
mapping<- mapping[c(1,2,4,6:8,9,10,12:15),c(1:2,4:6)]

data.frame=melt(mapping[,c(1,3,5)],id.vars = 'Sample')
data.frame$variable=factor(data.frame$variable,levels = c('Unmapped','UniquelyMapped')) 
data.frame=merge(data.frame,mapping,by = 'Sample')
data.frame$ratio <- round(data.frame$ratio,2)
p <- ggplot(data=data.frame, aes(x=Sample, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=UniquelyMapped, label=ratio), vjust=1.6, 
            color="white", size=6)+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,hjust=1),axis.text = element_text(size = 15))+
  guides(fill = guide_legend(title = NULL))+xlab('Sample')+ylab('Reads')+ 
  theme(legend.position = 'top')
p
ggsave("maping2.pdf",plot = p,width = 8,height = 8)
