filenames=$(ls /home1/liziyi/EPS/ChIP-seq/Sun_H3K9me3_ChIP-seq/trim_analysis/peaks/*/*rep1_peaks.bed)
for i in ${filenames}
do
ls $i
awk '{print ($3-$2)}' $i | awk '{sum += $1};END {print sum}'
done

ls /home1/liziyi/anno/WGBS/mm10_genome_length.bed
awk '{print $3}' /home1/liziyi/anno/WGBS/mm10_genome_length.bed | awk '{sum += $1};END {print sum}'