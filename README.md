## ChIP-seq PIPELINE
This pipeline performs the following tasks:
- Create an isolated environment for RNA-seq analysis
- Intallation and Reference Genomes
- Quality control on FastQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Adapter Trim([Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) 
- Alignment([STAR](https://github.com/alexdobin/STAR))
- QC reports for RNA-seq([RSeQC](http://rseqc.sourceforge.net))
- Quantifying gene expression([HTSeq-count](https://github.com/htseq/htseq))
- RPKMs and TPM (using edgeR or [StringTie](https://ccb.jhu.edu/software/stringtie/))
- Differential expression analysis for standard designs ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
- Visualization(R)
- Running the pepline 

## System requirements
- Linux/Unix
- Python
- R 


## Installation
We uses the Miniconda3 package management system to harmonize all of the software packages. 
Use the following commands to install Minicoda3：
``` bash
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
### create an isolated environment for RNA-seq
``` bash
$ conda create -n rna-seq
$ conda activate rna-seq
``` 

### install tools
Tools needed for this analysis are: R, samtools, FastQC, Trim Galore, STAR, RSeQC, stringtie, gffcompare, htseq-count. 
``` bash
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda install -c r r 
$ conda install -c bioconda samtools
$ conda install -c bioconda fastqc
$ conda install trim-galore
$ conda install STAR
$ conda install -c bioconda rseqc 
$ conda install -c bioconda htseq
$ conda install -c bioconda bioconductor-deseq2
$ conda install -c bioconda stringtie 
```

#### Genome files
Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC. In this example analysis we will use the mouse mm10 version of the genome from UCSC.
```bash
$ mkdir anno
$ cd anno
$ mkdir mm10
$ cd mm10
$ wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

```

Generate genome indexes files for STAR mapping
```bash
$ tar zxvf Mus_musculus_UCSC_mm10.tar.gz
$ STAR --runThreadN 30 --runMode genomeGenerate --genomeDir star_index_mm10 --genomeFastaFiles /Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /Annotation/Archives/archive-current/Genes/genes.gtf 
```
#### Download the public data

```
for ((i = 12;i<=15;i++)); #
do
fastq-dump --split-3 -O data/ SRR0020$i.sra.sra 
done
```

## Quality control on FastQ files 
FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. 

run FastQC interactively or using ht CLI, which offers the following options:
```bash
$ fastqc seqfile1 seqfile2 .. seqfileN
```

## Adapter Trim[OPTIONAL]
Use trim_glore to trim sequence adapter from the read FASTQ files.
```bash
$ trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2  \
            --gzip -o $input_data
```

## Alignment
Perform alignments with STAR to the genome and transcriptome.

```bash
$ STAR --runThreadN 10 --genomeDir ~/anno/mm10/ --readFilesCommand zcat --readFilesIn R1.fastq.gz R2.fastq.gz --outFileNamePrefix samplename   --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5
```
Visualization the mapping ratio by R
```bash
$ Rscript ~rcode/mapping.R
```

## QC reports for RNA-seq
RSeQC package comprehensively evaluate different aspects of RNA-seq experiments, such as sequence quality, GC bias, polymerase chain reaction bias, nucleotide composition bias, sequencing depth, strand specificity, coverage uniformity and read distribution over the genome structure. 

‘geneBody_coverage.py’ scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position. Finally, it generates a plot illustrating the coverage profile along the gene body.
```bash
$ samtools index input.sorted.bam 
$ geneBody_coverage.py -r hg19.housekeeping.bed -i test1.bam,test2.bam,test3.bam  -o output
```

'clipping_profile.py' calculate the distributions of clipped nucleotides across reads.This program is used to estimate clipping profile of RNA-seq reads from BAM or SAM file. Note that to use this funciton, CIGAR strings within SAM/BAM file should have ‘S’ operation (This means your reads aligner should support clipped mapping).
```bash
$ clipping_profile.py -i test1.bam -s "PE" -o out
```

## Quantifying gene expression
Counting reads in features with htseq-count
```bash
$ htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty yourfile_name.bam ~/reference/hisat2_reference/Homo_sapiens.GRCh38.86.chr_patch_hapl_scaff.gtf > counts.txt
```

Calculate RPKM and TPM by R
```bash
$ Rscript ~rcode/counts2rpkm_tpm.R
```

## Differential expression analysis
DE analysis is done using the DESeq2 Bioconductor package. It takes the merged raw read counts (from HTseq-count) as an input:
```bash
$ Rscript ~rcode/deseq2.R
```

## Running the pepline using the bash file
One command for running
```bash
$ bash rna_seq.bash
```
