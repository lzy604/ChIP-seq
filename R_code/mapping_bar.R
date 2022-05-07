library(reshape)
library(ggplot2)
dat_dir='./'
output='./mapping_ratio.pdf'

files=list.files(dat_dir)
sum_df=data.frame()
for (i in files){
  if (substr(i,start = nchar(i)-8,stop = nchar(i))=='final.out'){
    f=read.table(paste(dat_dir,i,sep = '/'),fill = T,sep = '\t')
    sam_name=unlist(strsplit(i,split = '.',fixed = T))[1]
    tmp_df=data.frame(sam_name = sam_name,Uniquely_mapped=as.numeric(as.character(f[8,2])),
                      total=as.numeric(as.character(f[5,2])),
                      map_ratio=(as.character(f[9,2])))
    sum_df=rbind(sum_df,tmp_df)
  }
}
sum_df$Unmapped=sum_df$total-sum_df$Uniquely_mapped # calculate unmapped reads
sum_df$map_ratio=gsub('\\%','',sum_df$map_ratio) # remove %

dd=melt(sum_df[,c(1,2,5)],id.vars = 'sam_name')
dd$variable=factor(dd$variable,levels = c('Unmapped','Uniquely_mapped')) 
dd=merge(dd,sum_df,by = 'sam_name')
p=ggplot(data=dd, aes(x=sam_name, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Uniquely_mapped, label=map_ratio), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,hjust=1))+
  guides(fill = guide_legend(title = NULL))+xlab('Sample')+ylab('Reads')
p
ggsave(output,plot = p,width = 8,height = 6)






