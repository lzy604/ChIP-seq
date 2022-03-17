## ChIP-seq PIPELINE
This pipeline performs the following tasks:
 * Intallation and Reference Genomes
 * Raw read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
 * Adapter trimming ([Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
 * Alignment ([BWA](https://sourceforge.net/projects/bio-bwa/files/) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
 * Peak calling([MACS2](https://github.com/taoliu/MACS))
 * Differential Peak calling([DiffBind](https://bioconductor.riken.jp/packages/3.9/bioc/html/DiffBind.html))
 * Motif identification([HOMER](http://homer.ucsd.edu/homer/download.html))
 * Inferring transcriptional regulators([LISA](https://github.com/qinqian/lisa)）
 * Peaks data visualisation ([IGV](https://software.broadinstitute.org/software/igv/))
 * Visualization(R)
 * Running the pepline

## ChIP-seq Data Standards 
1. Experiments should have two or more biological replicates
2. Each ChIP-seq experiment should have a corresponding input control experiment with matching run type, read length, and replicate structure.
3. Library complexity is measured using the Non-Redundant Fraction (NRF) and PCR Bottlenecking Coefficients 1 and 2, or PBC1 and PBC2. Preferred values are as follows: NRF>0.9, PBC1>0.9, and PBC2>10.
4. For narrow-peak histone experiments, each replicate should have 20 million usable fragments.（H3F3A, H3K27me3, H3K36me3, H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, H4K20me1）
5. For broad-peak histone experiments, each replicate should have 45 million usable fragments.(H2AFZ, H3ac, H3K27ac, H3K4me2, H3K4me3, H3K9ac)
6. H3K9me3 is an exception as it is enriched in repetitive regions of the genome. Compared to other broad marks, there are few H3K9me3 peaks in non-repetitive regions of the genome in tissues and primary cells. This results in many ChIP-seq reads that map to a non-unique position in the genome. Tissues and primary cells should have 45 million total mapped reads per replicate.
7. For transcription factor experiments, each replicate should have 10 million usable fragments.

## System requirements
- Linux/Unix
- Python
- R 

## Installation
We uses the Miniconda3 package management system to harmonize all of the software packages. 
Use the following commands to install Minicoda3：
``` bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Create an isolated environment for ChIP-seq
``` bash
conda create -n chip-seq
conda activate chip-seq
``` 

### Install tools
Tools needed for this analysis are: R, bedtools, FastQC, Trim Galore, MACS2, BWA, Lisa. 
``` bash
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
```

#### Genome files
Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC. In this example analysis we will use the mouse mm10 version of the genome from UCSC.
```bash

```

## Quality control on FastQ files 
FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. 

run FastQC interactively or using ht CLI, which offers the following options:
```bash
fastqc seqfile1 seqfile2 .. seqfileN
```

## Adapter Trim[OPTIONAL]
Use trim_glore to trim sequence adapter from the read FASTQ files.
```bash
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2  \
            --gzip -o $input_data
```

## Alignment
Reads were preprocessed and aligned reads to the appropriate reference genome (e.g., GRCh38 for human) using bowtie2, as shown below.

```bash
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
```

## Peak calling
Peak calling was performed using MACS2 (200bp each, ±100bp from the peak summit), followed by removal of peak regions overlapping any blacklisted regions.
```bash
## for narrow peaks:
macs2 callpeak -t IP.bam -c Input.bam -n test -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs

## for borad regions:
macs2 callpeak -t IP.bam -c Input.bam -n test --broad -p 0.01 --nomodel --extsize fragment_length --keep-dup all -g hs
```

## Motif discovery
HOMER contains two tools, findMotifs.pl and findMotifsGenome.pl, that manage all the steps for discovering motifs in promoter and genomic regions, respectively. 


```Shell
findMotifs.pl <inputfile.txt> <promoter set> <output directory> [options]

##findMotifsGenome.pl H3K4Me3.bed hg19 H3K4Me3_motif -len 8,10,12 
##This will search for motifs of length 8，10，and 12 from -400 to +100 relative to the TSS, using 4 threads (i.e. 4 CPUs)

##findMotifs.pl will produce a number of output files in the "output directory".  The primary output will be in HTML files that should be opened with you favorite web browser.

```

### Differential Peak calling
DiffBind is an R Bioconductor package that is used for identifying sites that are differentially enriched between two or more sample groups. 
```bash
Rscript ~rcode/DiffBind.R
```

## Running the pepline using [CHIPS(https://github.com/liulab-dfci/CHIPS)]

### Resources for ChIP-seq 
1. [ENCODE: Encyclopedia of DNA Elements](https://www.encodeproject.org/)  [ENCODExplorer](https://www.bioconductor.org/packages/release/bioc/html/ENCODExplorer.html): A compilation of metadata from ENCODE. A bioc package to access the meta data of ENCODE and download the raw files.
2. [ENCODE Factorbook](https://www.encodeproject.org/)  
3. [ChromNet ChIP-seq interactions](http://chromnet.cs.washington.edu/#/?search=&threshold=0.5)  
    paper: [Learning the human chromatin network using all ENCODE ChIP-seq datasets](http://biorxiv.org/content/early/2015/08/04/023911)  
4. [The International Human Epigenome Consortium (IHEC) epigenome data portal](http://epigenomesportal.ca/ihec/index.html?as=1)
5. [GEO](http://www.ncbi.nlm.nih.gov/gds/?term=). Sequences are in .sra format, need to use sratools to dump into fastq.
6. [European Nucleotide Archive](http://www.ebi.ac.uk/ena). Sequences are available in fastq format.
7. [Data bases and software from Sheirly Liu's lab at Harvard](http://liulab.dfci.harvard.edu/WEBSITE/software.htm)
8. [Blueprint epigenome](http://dcc.blueprint-epigenome.eu/#/home)
9. [A collection of tools and papers for nucelosome positioning and TF ChIP-seq](http://generegulation.info/)
10. [review paper:Deciphering ENCODE](http://www.cell.com/trends/genetics/fulltext/S0168-9525(16)00017-2)
11. [EpiFactors](http://epifactors.autosome.ru/) is a database for epigenetic factors, corresponding genes and products.
12. [biostar handbook](https://read.biostarhandbook.com/). My [ChIP-seq chapter](https://read.biostarhandbook.com/chip-seq/chip-seq-analysis.html) is out April 2017!
13. [ReMap 2018](http://tagc.univ-mrs.fr/remap/) An integrative ChIP-seq analysis of regulatory regions. The ReMap atlas consits of 80 million peaks from 485 transcription factors (TFs), transcription coactivators (TCAs) and chromatin-remodeling factors (CRFs) from public data sets. The atlas is available to browse or download either for a given TF or cell line, or for the entire dataset. 
