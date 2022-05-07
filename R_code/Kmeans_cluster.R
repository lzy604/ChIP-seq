

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (input file and output dir).n", call.=FALSE)
} else{
  input_file=args[1]
  out_dir=args[2]
}

# input_file='H3K9me3_1kb_peaks.txt'
# out_dir='./'

### filter peaks commonly identified across all conditions 
H3K9me3<-read.table(input_file,header=T)
H3K9me3$score=apply(H3K9me3,1,function(x)sum(as.numeric(x[3:ncol(H3K9me3)])))
H3K9me3_diff=H3K9me3[which(H3K9me3$score<(ncol(H3K9me3)-3)),]
H3K9me3_diff=H3K9me3_diff[,-ncol(H3K9me3_diff)]
### k-mean cluster
set.seed(100)
km<-kmeans(H3K9me3_diff[,3:ncol(H3K9me3_diff)],4)
H3K9me3_diff<-cbind(H3K9me3_diff,km$cluster);colnames(H3K9me3_diff)[ncol(H3K9me3_diff)]='km_cluster'
H3K9me3_diff<-H3K9me3_diff[order(H3K9me3_diff[,ncol(H3K9me3_diff)]),]
write.table(H3K9me3_diff,file.path(out_dir,'H3K9me3_peaks_kmeans4.txt'),sep='\t',quote=F,col.names=T,row.names = F)
