setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/Chip-seq/Sun_H3K9me3_ChIP/")
library(ggplot2)
library(plyr)
#read file
ceas.peak <- read.csv("meta.csv",row.names = 1)
metasheet <- read.csv("metasheet.csv")
metasheet <- metasheet[8:18,]
rownames(metasheet) <- metasheet$X.metasheet.csv..This.sheet.groups.your.analysis.into.runs.aka.experiments..
row.names(ceas.peak) <- metasheet$X[8:18]

#peak ratio
ceas.peak.ratio <- ceas.peak[,2:5]/ceas.peak$Total
#barplot(as.matrix(t(ceas.peak.ratio)),legend.text = TRUE,
#        axisnames = T,las=1, angle = 45,
#        axis.lty=0)
ceas.peak.ratio.barplot <- data.frame(sample=rep(row.names(ceas.peak.ratio), each=4),
                                      region=rep(colnames(ceas.peak.ratio),11),
                                      value =as.vector(t(ceas.peak.ratio)) )
pdf("ceas.peak.ratio.barplot.pdf",width = 15,height = 6)
ggplot(data=ceas.peak.ratio.barplot, aes(x=sample, y=value, fill=region)) +
  geom_bar(stat="identity")                                    
dev.off()
#####Barplot mean
peak.ratio.mean <- t(data.frame(apply(ceas.peak.ratio[1:3,], 2, mean),
                              apply(ceas.peak.ratio[4:6,], 2, mean),
                              apply(ceas.peak.ratio[7:9,], 2, mean),
                              apply(ceas.peak.ratio[10:11,], 2, mean)))
rownames(peak.ratio.mean) <- sub('..$','',row.names(ceas.peak.ratio)[c(1,4,7,10)])
peak.ratio.mean.plot <- data.frame(sample=rep(sub('..$','',row.names(ceas.peak.ratio)[1]), each=4),
                              region=colnames(ceas.peak.ratio),
                              value =apply(ceas.peak.ratio[1:3,], 2, mean))
tmp <- data.frame(sample=rep(sub('..$','',row.names(ceas.peak.ratio)[4]), each=4),
                 region=colnames(ceas.peak.ratio),
                 value =apply(ceas.peak.ratio[4:6,], 2, mean))
peak.ratio.mean.plot <- rbind(peak.ratio.mean.plot,tmp) 
tmp <- data.frame(sample=rep(sub('..$','',row.names(ceas.peak.ratio)[7]), each=4),
                  region=colnames(ceas.peak.ratio),
                  value =apply(ceas.peak.ratio[7:9,], 2, mean)) 
peak.ratio.mean.plot <- rbind(peak.ratio.mean.plot,tmp)   
tmp <- data.frame(sample=rep(sub('..$','',row.names(ceas.peak.ratio)[10]), each=4),
                  region=colnames(ceas.peak.ratio),
                  value =apply(ceas.peak.ratio[10:11,], 2, mean)) 
peak.ratio.mean.plot <- rbind(peak.ratio.mean.plot,tmp)  
pdf("ceas.peak.ratio.mean.barplot.pdf")
ggplot(data=peak.ratio.mean.plot, aes(x=sample, y=value, fill=region)) +
  geom_bar(stat="identity")        
dev.off()
##Barplot with error bars
barplot.error.bar <- ceas.peak.ratio.barplot
barplot.error.bar$sample <-sub('..$','',ceas.peak.ratio.barplot$sample)
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
                    groupnames=c("region", "sample"))
ggplot(error.bar, aes(x=sample, y=value, fill=region)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()

###repeat
library(stringr)
repeat.region <- as.data.frame(t(read.csv("repeat.csv",sep = ",",header = F)))
colnames(repeat.region) <- c("LTR","SINE","LINE","TOTAL")
row.names(repeat.region) <- str_extract(repeat.region$TOTAL,'Run.*.rep')
row.names(repeat.region) <- gsub(".rep1/R.*","",row.names(repeat.region))
repeat.region$TOTAL <-gsub("/home.*","",repeat.region$TOTAL)
repeat.region <- merge(repeat.region,metasheet,by=0)
row.names(repeat.region) <- repeat.region$X
repeat.region <- repeat.region[,2:5]
repeat.region <- repeat.region[ order(row.names(repeat.region)), ]
write.csv(repeat.region,"repeat.region.csv")

barplot.error.bar <- repeat.region
barplot.error.bar <-as.data.frame(apply(barplot.error.bar, 2, as.numeric))
barplot.error.bar$sample <-sub('..$','',row.names(repeat.region))
barplot.error.bar <- barplot.error.bar[,-(2:4)] 
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

error.bar <- data_summary(barplot.error.bar, varname="LTR", 
                          groupnames="sample")

barplot.error.bar <-as.data.frame(apply(barplot.error.bar, 2, as.numeric))

ggplot(error.bar, aes(x=sample, y=LTR)) + 
  geom_bar(stat="identity", fill="steelblue", position=position_dodge()) +
  geom_errorbar(aes(ymin=LTR-sd, ymax=LTR+sd), width=.2,position=position_dodge(.9))+ 
  scale_fill_brewer(palette="Paired") + theme_minimal()

repeat.region.ratio <-as.data.frame(apply(repeat.region, 2, as.numeric))
repeat.region.ratio <- repeat.region.ratio[,1:3]/(repeat.region.ratio$TOTAL)
rownames(repeat.region.ratio) <- row.names(repeat.region)

repeat.region.ratio.barplot <- data.frame(sample=rep(row.names(repeat.region.ratio), each=3),
                                      region=rep(colnames(repeat.region.ratio),11),
                                      value =as.vector(t(repeat.region.ratio)) )
pdf("repeat.region.ratio.barplot.pdf",width = 15,height = 6)
ggplot(data=repeat.region.ratio.barplot, aes(x=sample, y=value, fill=region)) +
  geom_bar(stat="identity")                                    
dev.off()

repeat.ratio.mean.plot <- repeat.region.ratio.barplot
repeat.ratio.mean.plot$sample <- sub('..$','',repeat.ratio.mean.plot$sample)
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
df3 <- data_summary(repeat.ratio.mean.plot, varname="value", 
                    groupnames=c("sample", "region"))
pdf("ceas.peak.ratio.mean.barplot.pdf")
ggplot(data=df3, aes(x=sample, y=value, fill=region)) +
  geom_bar(stat="identity",position=position_dodge())+ 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9))        
dev.off()
pdf("ceas.peak.ratio.mean.barplot1.pdf")
ggplot(data=repeat.ratio.mean.plot, aes(x=region, y=value, fill=sample)) +
  geom_bar(stat="identity")     
dev.off()


ggbarplot(repeat.ratio.mean.plot[repeat.ratio.mean.plot$region =="LINE",], x="sample", y="value", add = "mean_se")+
  stat_compare_means() +                                         # Global p-value
  stat_compare_means(ref.group = "ESC", label = "p.signif",label.y = c(0.5, 0.6))


ggbarplot(repeat.ratio.mean.plot[repeat.ratio.mean.plot$region =="SINE",], x="sample", y="value", add = "mean_se")+
  stat_compare_means() +                                         # Global p-value
  stat_compare_means(ref.group = "ESC", label = "p.signif",label.y = c(0.15, 0.2))

ggbarplot(repeat.ratio.mean.plot[repeat.ratio.mean.plot$region =="LTR",], x="sample", y="value", add = "mean_se")+
  stat_compare_means() +                                         # Global p-value
  stat_compare_means(ref.group = "ESC", label = "p.signif",label.y = c(0.4, 0.5))




                                   