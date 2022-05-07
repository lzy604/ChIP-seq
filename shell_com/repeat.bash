filenames=$(ls /home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/*/*rep1_peaks.bed)
for i in ${filenames}
do
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.promoter.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.HCP.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.ICP.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.LCP.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.exon.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.intron.bed -u | wc -l
 awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/repeat/LTR.bed -u -u | wc -l
 awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/repeat/SINE.bed -u | wc -l
 awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/repeat/LINE.bed -u | wc -l
 # awk '{print $1"\t"$2"\t"$2+1}' $i | intersectBed -wa -a - -b $HOME/annotations/hg38.telomere.bed -u | wc -l
wc -l $i
done