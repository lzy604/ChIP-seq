#read files
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/")
#read file
mapping <- read.csv("mapping.csv",row.names = 1)
peakStats <- read.csv("peakStats.csv",row.names = 1)
meta.data <- read.csv("metasheet.csv",row.names = 1)
#sample with peak and reads
meta.data <- meta.data[-c(1:7),]
rownames(peakStats) <- gsub(".rep1","",rownames(peakStats))
peakStats <- merge(peakStats,meta.data,by = 0)
rownames(peakStats) <- peakStats$X
sample.peak <- data.frame(peak=peakStats[,2])
rownames(sample.peak) <- rownames(peakStats)
sample.mapping <- merge(sample.peak,mapping,by=0)
row.names(sample.mapping) <- sample.mapping$Row.names
sample.mapping$peak.read<- sample.mapping$peak/sample.mapping$UniquelyMapped*1e6

##Barplot with error bars
barplot.error.bar <- data.frame(peak.read=sample.mapping$peak.read,row.names = row.names(sample.mapping))
barplot.error.bar$sample <-sub('..$','',row.names(sample.mapping))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
library(plyr)
error.bar <- data_summary(barplot.error.bar, varname="peak.read", 
                          groupnames="sample")
library(ggplot2)
pdf("peak.per.reads.pdf")
ggplot(error.bar, aes(x=sample, y=peak.read)) + 
  geom_bar(stat="identity", fill="steelblue",position=position_dodge()) +
  geom_errorbar(aes(ymin=peak.read-sd, ymax=peak.read+sd), width=.2,position=position_dodge(.9))+
  theme_minimal()
dev.off()
##