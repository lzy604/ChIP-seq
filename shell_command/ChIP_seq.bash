#Use the following commands to install Minicoda3：
 
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh


### Create an isolated environment for ChIP-seq
conda create -n chip-seq
conda activate chip-seq
 

### Install tools
Tools needed for this analysis are: R, bedtools, FastQC, Trim Galore, MACS2, BWA, Lisa. 
 
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c r r 
conda install -c bioconda bedtools 
conda install -c bioconda fastqc
conda install trim-galore
conda install -c bioconda macs2
conda install -c bioconda bwa 
conda install -c liulab-dfci lisa2 
conda install -c bioconda homer 


#### Genome files



## Quality control on FastQ files 
fastqc seqfile1 seqfile2 .. seqfileN


## Adapter Trim[OPTIONAL]

trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2  \
            --gzip -o $input_data


## Alignment
#set global variables
GENOME=/anno/human/genome/Homo_sapiens.GRCh38.dna.primary_assemblyexport
BASE_PATH=/home/user/human/ChIP-seq/
SUFFIX=_trim_paired.fastq.gz
cd $BASE_PATH/GENE_ID/replicate1

# replace _read_ID_ with your IDs
R1=QC/_read_ID_$SUFFIX
R2=QC/_read_ID_$SUFFIX

# do alignment
bowtie2 -p 4 -x $GENOME -1 $R1 -2 $R2 -S alignment.sam 2> log.txt

# if single-read, use
READ=QC/_read_ID_$SUFFIX
bowtie2 -p 4 -x $GENOME -U $READ -S alignment.sam 2> log.txt

#remove unmapped reads and duplicated reads (268= Read unmapped (4) or  Mate unmapped (8) or  Not primary alignment (256))
samtools view -h -F 268 -q 5 -bS alignment.sam > unique_alignment.bam
samtools sort unique_alignment.bam -o unique_alignment_sorted.bam
samtools rmdup unique_alignment_sorted.bam unique_alignment_sorted_rd.bam
rm alignment.sam unique_alignment.bam unique_alignment_sorted.bam


## Peak calling
Peak calling was performed using MACS2 (200bp each, ±100bp from the peak summit), followed by removal of peak regions overlapping any blacklisted regions.

## for narrow peaks:
macs2 callpeak -t IP.bam -c Input.bam -n test -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs

## for borad regions:
macs2 callpeak -t IP.bam -c Input.bam -n test --broad -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs


## Motif discovery
findMotifs.pl <inputfile.txt> <promoter set> <output directory> [options]

##findMotifsGenome.pl H3K4Me3.bed hg19 H3K4Me3_motif -len 8,10,12 
##This will search for motifs of length 8，10，and 12 from -400 to +100 relative to the TSS, using 4 threads (i.e. 4 CPUs)

##findMotifs.pl will produce a number of output files in the "output directory".  The primary output will be in HTML files that should be opened with you favorite web browser.



### Differential Peak calling
DiffBind is an R Bioconductor package that is used for identifying sites that are differentially enriched between two or more sample groups. 

Rscript ~rcode/DiffBind.R
