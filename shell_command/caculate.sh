file=$(ls *deseq2_sig0.05.bed)
for i in ${file}
do
sample=$(basename $i .bed )
awk '{if($5 >= 0){print}}' $i > ${sample}_pfc.bed
awk '{if($5 <= 0){print}}' $i > ${sample}_nfc.bed

wc -l ${sample}_nfc.bed
awk '{print $1"\t"$2"\t"$3}' ${sample}_nfc.bed | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/LTR.bed -u  | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_nfc.bed| intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/LINE.bed -u | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_nfc.bed | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/SINE.bed -u  | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_nfc.bed| intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/imprinting.gene.bed -u | wc -l

wc -l ${sample}_pfc.bed
awk '{print $1"\t"$2"\t"$3}' ${sample}_pfc.bed | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/LTR.bed -u | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_pfc.bed| intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/LINE.bed -u | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_pfc.bed | intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/repeat/SINE.bed -u | wc -l

awk '{print $1"\t"$2"\t"$3}' ${sample}_pfc.bed| intersectBed -wa -a - -b /home1/liziyi/anno/mm10.ucsc/bed/imprinting.gene.bed -u | wc -l
done
