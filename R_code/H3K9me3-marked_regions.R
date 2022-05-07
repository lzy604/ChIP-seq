setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/")
library(ggplot2)
library(plyr)
#read file
tmp <- read.csv("H3K9me3-marked_regions.out",header = F)
H3K9me3_marked_regions <- as.data.frame(tmp[seq(2, nrow(tmp), 2),])
colnames(H3K9me3_marked_regions) <- "length"
row.names(H3K9me3_marked_regions) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/','', tmp[seq(1, nrow(tmp), 2),])
row.names(H3K9me3_marked_regions) =sub('/.*','', row.names(H3K9me3_marked_regions))
H3K9me3_marked.region.ratio <-as.data.frame(apply(H3K9me3_marked_regions, 2, as.numeric))
H3K9me3_marked.region.ratio$ratio <- H3K9me3_marked.region.ratio$length/(H3K9me3_marked.region.ratio$length[nrow(H3K9me3_marked.region.ratio)])
row.names(H3K9me3_marked.region.ratio) <- row.names(H3K9me3_marked_regions)

tmp2 <- read.csv("H3K9me3-marked_regions_overlap_bed.out",header = F)
H3K9me3_marked_regions2 <- as.data.frame(tmp2[seq(2, nrow(tmp2), 2),])
colnames(H3K9me3_marked_regions2) <- "length"
row.names(H3K9me3_marked_regions2) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bedfiles/','', tmp2[seq(1, nrow(tmp2), 2),])

H3K9me3_marked.region2.ratio <-as.data.frame(apply(H3K9me3_marked_regions2, 2, as.numeric))
H3K9me3_marked.region2.ratio$ratio <- H3K9me3_marked.region2.ratio$length/(H3K9me3_marked.region.ratio$length[nrow(H3K9me3_marked.region.ratio)])
row.names(H3K9me3_marked.region2.ratio) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/bedfiles/','', tmp2[seq(1, nrow(tmp2), 2),])


#plot
H3K9me3_marked.region.ratio.barplot <- data.frame(sample=row.names(H3K9me3_marked_regions)[-nrow(H3K9me3_marked.region.ratio)],
                                                  value =H3K9me3_marked.region.ratio$ratio[-nrow(H3K9me3_marked.region.ratio)])
ggplot(data=H3K9me3_marked.region.ratio.barplot, aes(x=sample, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()

H3K9me3_marked.region.ratio.barplot2 <- data.frame(sample=row.names(H3K9me3_marked.region2.ratio),
                                                  value =H3K9me3_marked.region2.ratio$ratio)
ggplot(data=H3K9me3_marked.region.ratio.barplot2, aes(x=sample, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()

##Barplot with error bars
barplot.error.bar <- H3K9me3_marked.region.ratio.barplot
barplot.error.bar$value <- barplot.error.bar$value*100
barplot.error.bar$sample <-sub('..$','',barplot.error.bar$sample)
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
error.bar <- data_summary(barplot.error.bar, varname="value", 
                          groupnames="sample")

pdf("H3K9me3_marked.region.ratio.barplot.pdf")
ggplot(error.bar, aes(x=sample, y=value)) + 
  geom_bar(stat="identity", fill="steelblue", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()+
  ylab("H3K9me3-marked regions 
       compared to the genome (%)")+
  theme_classic()
dev.off()
####LINE
tmp <- read.csv("LINE.region.out",header = F)
LINE_marked_regions <- as.data.frame(tmp[seq(2, nrow(tmp), 2),])
colnames(LINE_marked_regions) <- "length"
row.names(LINE_marked_regions) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/','', tmp[seq(1, nrow(tmp), 2),])
row.names(LINE_marked_regions) =sub('/.*','', row.names(LINE_marked_regions))
LINE_marked.region.ratio <-as.data.frame(apply(LINE_marked_regions, 2, as.numeric))
LINE_marked.region.ratio$ratio <- LINE_marked.region.ratio$length/(548842036)
#plot
LINE_marked.region.ratio.barplot <- data.frame(sample=row.names(LINE_marked_regions)[-nrow(LINE_marked.region.ratio)],
                                                  value =LINE_marked.region.ratio$ratio[-nrow(LINE_marked.region.ratio)])
ggplot(data=LINE_marked.region.ratio.barplot, aes(x=sample, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
##Barplot with error bars
barplot.error.bar <- LINE_marked.region.ratio.barplot
barplot.error.bar$sample <-sub('..$','',barplot.error.bar$sample)
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
error.bar <- data_summary(barplot.error.bar, varname="value", 
                          groupnames="sample")
pdf("LINE_marked.region.ratio.barplot.pdf")
ggplot(error.bar, aes(x=sample, y=value)) + 
  geom_bar(stat="identity", fill="steelblue", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()
dev.off()


#####SINE
tmp <- read.csv("SINE.region.out",header = F)
SINE_marked_regions <- as.data.frame(tmp[seq(2, nrow(tmp), 2),])
colnames(SINE_marked_regions) <- "length"
row.names(SINE_marked_regions) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/','', tmp[seq(1, nrow(tmp), 2),])
row.names(SINE_marked_regions) =sub('/.*','', row.names(SINE_marked_regions))
SINE_marked.region.ratio <-as.data.frame(apply(SINE_marked_regions, 2, as.numeric))
SINE_marked.region.ratio$ratio <- SINE_marked.region.ratio$length/(206172774)
#plot
SINE_marked.region.ratio.barplot <- data.frame(sample=row.names(SINE_marked_regions)[-nrow(SINE_marked.region.ratio)],
                                               value =SINE_marked.region.ratio$ratio[-nrow(SINE_marked.region.ratio)])
ggplot(data=SINE_marked.region.ratio.barplot, aes(x=sample, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
##Barplot with error bars
barplot.error.bar <- SINE_marked.region.ratio.barplot
barplot.error.bar$sample <-sub('..$','',barplot.error.bar$sample)
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
error.bar <- data_summary(barplot.error.bar, varname="value", 
                          groupnames="sample")
pdf("SINE_marked.region.ratio.barplot.pdf")
ggplot(error.bar, aes(x=sample, y=value)) + 
  geom_bar(stat="identity", fill="steelblue", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()
dev.off()


####LTR

tmp <- read.csv("LTR.region.out",header = F)
LTR_marked_regions <- as.data.frame(tmp[seq(2, nrow(tmp), 2),])
colnames(LTR_marked_regions) <- "length"
row.names(LTR_marked_regions) =sub('/home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/','', tmp[seq(1, nrow(tmp), 2),])
row.names(LTR_marked_regions) =sub('/.*','', row.names(LTR_marked_regions))
LTR_marked.region.ratio <-as.data.frame(apply(LTR_marked_regions, 2, as.numeric))
LTR_marked.region.ratio$ratio <- LTR_marked.region.ratio$length/(320032086)
#plot
LTR_marked.region.ratio.barplot <- data.frame(sample=row.names(LTR_marked_regions)[-nrow(LTR_marked.region.ratio)],
                                               value =LTR_marked.region.ratio$ratio[-nrow(LTR_marked.region.ratio)])
ggplot(data=LTR_marked.region.ratio.barplot, aes(x=sample, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
##Barplot with error bars
barplot.error.bar <- LTR_marked.region.ratio.barplot
barplot.error.bar$sample <-sub('..$','',barplot.error.bar$sample)
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
error.bar <- data_summary(barplot.error.bar, varname="value", 
                          groupnames="sample")
pdf("LTR_marked.region.ratio.barplot.pdf")
ggplot(error.bar, aes(x=sample, y=value)) + 
  geom_bar(stat="identity", fill="steelblue", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()
dev.off()
