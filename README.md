## ChIP-seq PIPELINE
This pipeline performs the following tasks:
  * [Data download](#1-data-download)
  * [Control selection](#2-control-selection)
  * [Read quality control](#3-read-quality-control)
  * [Read mapping](#4-read-mapping)
  * [Peak calling](#5-peak-calling)
  * [Motif discovery](#6-motif-discovery)
  * [Replicate quality control and merging](#7-replicate-quality-control)
  * [Representative motif selection](#8-representative-motif-selection)
  * [Resources for ChIP-seq](#) 

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
#### create an isolated environment for RNA-seq
``` bash
conda create -n rna-seq
conda activate rna-seq
``` 

#### install tools
Tools needed for this analysis are: R, samtools, FastQC, Trim Galore, STAR, RSeQC, stringtie, gffcompare, htseq-count. 
``` bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c r r 
conda install -c bioconda samtools
conda install -c bioconda fastqc
conda install trim-galore
conda install STAR
conda install -c bioconda rseqc 
conda install -c bioconda htseq
conda install -c bioconda bioconductor-deseq2
conda install -c bioconda stringtie 
```

### (1) Data download

Raw transcription factor ChIP-seq FASTQ reads were downloaded directly from ENCODE. For example, for *GATA1*, to obtain the SE36nt reads associated with Accession ID ENCFF000YND (experiment ENCSR000EFT, library ENCLB209AJT), the following URL was used:

`https://www.encodeproject.org/files/ENCFF000YND/@@download/ENCFF000YND.fastq.gz`

### (2) Control selection

Controls were selected as described in the manuscript, with preference for type input DNA, as shown in the figure below.

![Choosing controls](images/choosing_controls.png)

### (3) Read quality control

Read quality control (QC) was performed using Trimmomatic, with command-line argument values as decribed in the manuscript. For paired-end (PE) reads, the following commands were used:

```Shell
#set global variables
BASE=/home/your/working/directory
Trimmomatic_Path=/where/is/Trimmomatic
Trimmomatic='java -jar /path/Trimmomatic-0.36/trimmomatic-0.36.jar'
FASTQC=/where/is/fastqc

# replace LGJ20-XQ41_R1_001 and LGJ20-XQ41_R2_001 with your read IDs
cd $BASE/Sample_one/
R1=LGJ20-XQ41_R1_001.fastq.gz
R2=LGJ20-XQ41_R2_001.fastq.gz
mkdir QC fastqc

# QC
$Trimmomatic PE -threads 1 $R1 $R2 QC/${R1%.fastq.gz}_trim_paired.fastq.gz QC/${R1%.fastq.gz}_trim_unpaired.fastq.gz QC/${R2%.fastq.gz}_trim_paired.fastq.gz QC/${R2%.fastq.gz}_trim_unpaired.fastq.gz ILLUMINACLIP:$Trimmomatic_Path/adapters/TruSeq3-PE-2.fa:2:40:12:8:true LEADING:10 SLIDINGWINDOW:4:15 MINLEN:50 2> read_processing.log

# check read quality
$FASTQC QC/*trim_paired.fastq.gz -o fastqc
```

For single-end (SE) reads, the Trimmomatic call was replaced with the following:

```Shell
$Trimmomatic SE -threads 1 $READ QC/${READ%.fastq.gz}_trim_paired.fastq.gz ILLUMINACLIP:/home/cpyu/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:40:12 LEADING:10 SLIDINGWINDOW:4:15 MINLEN:30 2> read_processing.log
```

### (4) Read mapping

Reads were preprocessed and aligned reads to the appropriate reference genome (e.g., GRCh38 for human) using bowtie2, as shown below.

```Shell
#set global variables
GENOME=/home/user/human/genome/Homo_sapiens.GRCh38.dna.primary_assemblyexport
BASE_PATH=/home/user/human/ChIP-seq/
SUFFIX=_trim_paired.fastq.gz
cd $BASE_PATH/GENE_ID/replicate1

# replace _read_ID_ with your IDs
R1=QC/_read_ID_$SUFFIX
R2=QC/_read_ID_$SUFFIX

# do alignment
time bowtie2 -p 4 -x $GENOME -1 $R1 -2 $R2 -S alignment.sam 2> log.txt

# if single-read, use
READ=QC/_read_ID_$SUFFIX
time bowtie2 -p 4 -x $GENOME -U $READ -S alignment.sam 2> log.txt

#remove unmapped reads and duplicated reads (268= Read unmapped (4) or  Mate unmapped (8) or  Not primary alignment (256))
samtools view -h -F 268 -q 5 -bS alignment.sam > unique_alignment.bam
samtools sort unique_alignment.bam -o unique_alignment_sorted.bam
samtools rmdup unique_alignment_sorted.bam unique_alignment_sorted_rd.bam
rm alignment.sam unique_alignment.bam unique_alignment_sorted.bam
```

### (5) Peak calling

Peak calling was performed using MACS2 (200bp each, ±100bp from the peak summit), followed by removal of peak regions overlapping any blacklisted regions.

![Motif discovery blacklist](images/motif_discovery_blacklist.png)

```Shell
#set global variables
CURR=/home/user/human/ChIP-seq/
cd $CURR/GENE_ID/replicate1
CONTROL="control1.bam control2.bam" # control string may include multiple files

# individual replicate, paired-end (PE) reads
macs2 callpeak -t unique_alignment_sorted_rd.bam -c $CONTROL -f BAMPE --gsize hs --outdir macs2

# individual replicate, single-read (SE) reads
macs2 callpeak -t unique_alignment_sorted_rd.bam -c $CONTROL -f BAM --gsize hs --outdir macs2

## merging two replicates, PE reads
macs2 callpeak -t replicate1.bam replicate2.bam -c $CONTROL $ CONTROL -f BAMPE --gsize hs --outdir macs2

# mergin two replicates, SE reads 
macs2 callpeak -t replicate1.bam replicate2.bam -c $CONTROL $ CONTROL -f BAM --gsize hs --outdir macs2
```

### (6) Motif discovery

Motif discovery was performed using MEME-chip as described in the manuscript, using the top 500 peaks to determine five motifs per analysis.

```Shell
#set global variables
CURR=/home/user/human/ChIP-seq/
GENOME=/home/user/human/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
BlackList=/home/human/blacklist/ENCFF023CZC_sorted.bed
NoTopPeak=500 # number of top peaks to analyze

# navigate to replicate directory, extract peaks
cd $CURR/GENE_ID/replicate1/
awk 'BEGIN {WIDTH=100} {if($2<=WIDTH) print $1 "\t1\t" $2+WIDTH "\t" $4 "\t" $5; else print $1 "\t" $2-WIDTH "\t" $2+WIDTH "\t" $4 "\t" $5}' macs2/NA_summits.bed > extended_peaks.bed

# remove blacklist regions from peaks
bedtools subtract -a extended_peaks.bed -b $BlackList -A > extended_bk_removal_peaks.bed
sort -r -k5 -n extended_bk_removal_peaks.bed| head -n $NoTopPeak > top_peaks.bed
bedtools getfasta -bed top_peaks.bed -fi $GENOME > top_peaks.fa

# motif discovery
meme-chip -meme-nmotifs 5 top_peaks.fa
```

### (7) Replicate quality control and merging

Replicates were merged as described in the manuscript. Briefly, replicates were required to yield:

1. at least one position weight matrix (PWM) supported by >100 peaks and an E-value of <0.0001
2. and IDR score (-log[IDR]) ≥1.5 when compared to the other replicate of the same experiment

Experiments were required to have at least 2 passing replicates, and excluded otherwise.

``` Shell
# normalizing scores for each peak
python score_quantile.py -s sample1_bk_removal.bed -g sample1 -o sample1_score.bed # repeat for all samples

# concatenating all scoring peaks and grouping peaks into clusters
cat sample1_score.bed > total.bed
cat sample2_score.bed >> total.bed #... for all samples

# sort and cluster
bedtools sort -i total.bed > temp.bed
bedtools cluster -i temp.bed > total_cluter.bed
rm *_score.bed total.bed temp.bed

# select top 500 peaks
python top_peak_sel.py -i total_cluter.bed -o top_combined.bed

# calculate motif position weight matrix (PWM) correlations
python correlation.py -m1 combined.meme -o motif_pcc.txt

# <<<<<<< Updated upstream
python motif_cluster.py motif_pcc.txt occurrences.txt
```

### (8) Representative motif selection

**Single experiment**.
For transcription factors represented by **one passing experiment** in the ENCODE data, replicates were merged and step 5 ([peak calling](#peak-calling)) was repeated using the merged replicates. Final representative motifs were called by repeating and 6 ([motif discovery](#motif-discovery)) using these merged-replicate peaks.

**!!TODO** *CP or CH please check that the BELOW section is correct*

**Multiple experiments**.
For transcription factors represented by **more than one passing experiment** from **more than one biosample (*i.e.*, cell line or tissue type)** in the ENCODE data, one experiment was selected to represent each biosample. To do this, we first computed the similarities (KFV cos; Xu and Su 2010) between all pairs of motifs (*i.e.*, position weight matrices; PWMs) for all pairs of experiments utilizing different biosamples. For each biosample, the experiment was selected which shows the highest KFV cos value with another PWM from another experiment using a different biosample.

For TFs with two biosamples, if neither biosample’s representative experiment had at least one PWM with cos>=0.80, we selected the experiment whose top PWM had the highest number of occurrences in peak regions. For TFs with >2 biosamples, experiments lacking at least one motif with KFV cos>=0.80 were excluded. Finally, using the representative experiments for all available biosamples, the ranking method described in the manuscript was used to select the top 500 peaks across all experiments. Final representative motifs were called by repeating step 6 ([motif discovery](#motif-discovery)) using these top peaks.

For transcription factors represented by **more than one passing experiment** from **only one biosample** in the ENCODE data, different experiments were not compared. Instead, we directly used the ranking method described in the manuscript to select the top 500 peaks across all experiments. Final representative motifs were called by repeating step 6 ([motif discovery](#motif-discovery)) using these top peaks.

## Software and Data Versions Used

* Genome
  * Human: GRCh38
  * Mouse: GRCm38
* Blacklists
  * Amemiya et al. 2019
* FASTQC v0.11.8
* TRIMMOMATIC v0.39
  * Bolger et al. 2014
* BOWTIE2 v2.3.5 (64-bit)
  * Langmead et al. 2012
* MACS2 v2.1.2
  * <https://github.com/taoliu/MACS>
  * Zhang et al. 2008
* MEME_CHIP v5.0.5
  * Bailey et al. 2009
* SAMTOOLS v1.9
  * using htslib 1.9
  * Li et al. 2009
* BEDTOOLS v2.28.0
  * Quinlan et al. 2010
* IDR
  * <https://github.com/nboley/idr>
  * Li et al. 2011

## Analysis Scripts

In addition to published tools, our analyses utilized the following custom scripts:
(Yu: I'll add descriptions, Chen-Hao: Chen-Hao can help!)

1. `pwm.py`. **!!TODO: please check this is correct. Can the input have only one, or multiple PWMs? CH: Multiple PWMs in single MEME file**

    * **Decription**: Trims PWMs from their ends until the remaining termini have an information content greater than or equal to some threshold value.
    * **Input**: (1) `-i`, PWMs in MEME format; (2) `--trim`, threshold information content; and (3) `-o` path/name of output file.
    * **Output**: A MEME file containing trimmed PWMs.
    * **Example**:
      * `python3.6 pwm.py -i myPWM.meme --trim 0.3 -o myPWM_trimmed.meme`

2. `motif_cluster.py`. **!!TODO: please check this is correct.**

    * **Description**: Determines groups (clusters) of similar PWMs based on their correlations.
    * **Input**: (1) the motif_pcc.txt files produced by `correlation.py` (first argument, unnamed); and (2) **!!TODO??**
    * **Output**: Groups (clusters) of PWMs... **!!TODO??**. Printed to STOUT, here redirected to the file `motif_cluster.txt`.
      * `python3.6 motif_cluster.py motif_pcc.txt occurrences.txt > motif_cluster.txt`

3. `consensus_pwm.py`. **!!TODO** please add

4. `correlation.py`. **!!TODO: please check this is correct.**

    * **Description**: Calculates similarity between PWMs in MEME format with user-defined correlation methods described in KFV.
    * **Input**: (1) `--method`, method for calculating correlation (cos, pcc, eucl, and kl); (2) `-m1`, PWMs in MEME format (if one PWM file, calculates pairwise correlations within the file; if two PWM files, calculates pairwise correlation between the two files); and (3) `-o` path/name of output file.
    * **Output**: a plain text file giving the similarity between each pair of motifs.
    * **Example**:
      * `python3.6 correlation.py --method pcc -m1 PWM_trimmed.meme -o motif_pcc.txt`

5. `peak_motif_ranges.R`. **!!TODO** please add

6. `score_quantile.py`. **!!TODO** please add

7. `top_peak_sel.py`. **!!TODO** please add

## Supplementary Perl Scripts

Exact commands for running the pipeline for a single transcription factor (TF) are provided in the following Perl scripts:

1. `chip_seq_download.pl`. This script downloads a user-provided list of ENCODE experiments, *i.e.*, the FASTQ data from ENCODE TF ChIP-seq assays. It should be called from the directory you wish to populate with directories and subdirectories containing the data. Input should be provided as a TAB-delimited input file with the following six (named) columns:

    * `Target name`. The name of the transcription factor (*e.g.*, **adb-A**).
    * `directory_name`. The name of the directory to create in which to store the transcription factor's data (*e.g.*, **adb_A**)
    * `AC`. Accession ID of the relevant experiment (*e.g.*, **ENCSR609JDR**). Note that one experiment will be associated with more than one replicate and control, i.e., each Accession ID will have more than one row in the file. A new subdirectory within `directory_name` will be created for each unique `AC`.
    * `type`. Must contain the value `replicate` or `control`.
    * `Library`. Library (assay) accession ID for the specific replicate or control (e.g., **ENCLB367MZH**). Individual assays that are used as controls will have `_control` appended to the name of their `Library` subdirectory.
    * `url`. The ENCODE url for direct download of the FASTQ data associated with this `AC`/`Library` (*e.g.*, **<https://www.encodeproject.org/files/ENCFF036RZV/@@download/ENCFF036RZV.fastq.gz>**). This can be constructed by appending the FASTQ file ID to the end of `https://www.encodeproject.org/files/ENCFF036RZV/@@download/`.

    For example, a short file for downloading some *Drosophila* FASTQ data is shown below, followed by a figure explaining the subdirectory structure that will be created in the process. This structure will be expected by the subsequent scripts (*e.g.*, `chip_seq_pipeline.pl`). Note that replicate and control subdirectories will contain one FASTQ file for single-end (SE) reads and two FASTQ files for paired-end (PE) reads; this is how the pipeline differentiates between the two.

    ``` CSV
    Target name directory_name AC type Library url
    abd-A abd_A ENCSR609JDR replicate ENCLB367MZH https://www.encodeproject.org/files/ENCFF036RZV/@@download/ENCFF036RZV.fastq.gz
    abd-A abd_A ENCSR609JDR replicate ENCLB506QBD https://www.encodeproject.org/files/ENCFF999PWE/@@download/ENCFF999PWE.fastq.gz
    abd-A abd_A ENCSR609JDR replicate ENCLB589GMA https://www.encodeproject.org/files/ENCFF646QDZ/@@download/ENCFF646QDZ.fastq.gz
    abd-A abd_A ENCSR609JDR control ENCLB428GIT_control https://www.encodeproject.org/files/ENCFF120UZS/@@download/ENCFF120UZS.fastq.gz
    Abd-B Abd_B ENCSR465HPZ replicate ENCLB009GTL https://www.encodeproject.org/files/ENCFF469WCT/@@download/ENCFF469WCT.fastq.gz
    Abd-B Abd_B ENCSR465HPZ replicate ENCLB205LSV https://www.encodeproject.org/files/ENCFF969SEG/@@download/ENCFF969SEG.fastq.gz
    Abd-B Abd_B ENCSR465HPZ control ENCLB730TLW_control https://www.encodeproject.org/files/ENCFF730IEL/@@download/ENCFF730IEL.fastq.gz
    Abd-B Abd_B ENCSR692UBK replicate ENCLB526KET https://www.encodeproject.org/files/ENCFF820QJV/@@download/ENCFF820QJV.fastq.gz
    Abd-B Abd_B ENCSR692UBK replicate ENCLB588NTL https://www.encodeproject.org/files/ENCFF252IVG/@@download/ENCFF252IVG.fastq.gz
    Abd-B Abd_B ENCSR692UBK control ENCLB728VIS_control https://www.encodeproject.org/files/ENCFF354CHN/@@download/ENCFF354CHN.fastq.gz
    achi achi ENCSR959SWC replicate ENCLB061SFK https://www.encodeproject.org/files/ENCFF979YUJ/@@download/ENCFF979YUJ.fastq.gz
    achi achi ENCSR959SWC replicate ENCLB064AEO https://www.encodeproject.org/files/ENCFF881ITO/@@download/ENCFF881ITO.fastq.gz
    achi achi ENCSR959SWC replicate ENCLB722DLY https://www.encodeproject.org/files/ENCFF721FQK/@@download/ENCFF721FQK.fastq.gz
    achi achi ENCSR959SWC control ENCLB240LXI_control https://www.encodeproject.org/files/ENCFF548BRW/@@download/ENCFF548BRW.fastq.gz
    ```

    ![Directory structure](images/directory_structure.png)

    ![Example](images/example.png)

2. `chip_seq_pipeline.pl`. This script can be used after the experimental data have been downloaded using `chip_seq_download.pl`. Alternatively, the user may download their own data, provided they have placed them in directories precisely matching the structure described above.

    This script **must be called from within the directory corresponding to a single Target (TF)**. This means the user must first navigate to one of the `directory_name` values specified above, e.g., abd_A (input file example above) or ATF3 (input figure example above). This means that multiple TFs can be run in parallel, but that a single TF cannot. The reason for this choice is that several steps in the pipeline (e.g., BOWTIE2) allow the user to specify multiple CPUs for further parallelism.

    This script navigates the directory structure starting at this point, carrying out the first steps of the pipeline, proceeding from quality control (Trimmomatic; FASTQC) to read mapping (BOWTIE2) to peak calling (MACS2) to motif discovery (MEME-chip). **Of utmost importance, the user must alter the code block at the top of the script, `### MANUALLY SET GLOBAL VARIABLES ###`, to set paths to each tool installed on their system.** Specifically, the user must provide paths or values from the following:

    * `FASTQC` (v0.11.8): path to software
    * `TRIMMOMATIC` (v0.39): path to software
    * `BOWTIE2` (v2.3.5 / 64-bit): path to software
    * `SAMTOOLS` (v1.9 using htslib 1.9): path to software
    * `BEDTOOLS` (v2.28.0): path to software
    * `MACS2` (v2.1.2): path to software
    * `MACS2_gsize`: a value, e.g. `hs` for *Homo sapeins* or `dm` for *Drosophila melanogaster*
    * `PEAK_RADIUS`: a value, the radius of read calling peaks, *e.g.*, `100` to extend peaks 100 bp in either direction from the peak summit.
    * `MIN_PEAK_SCORE`: a value, the minimum peak score required, *e.g.*, `13`.
    * `MEME_CHIP` (v5.0.5): path to software
    * `CCUT`: a value, the size (bp) to which peak regions should be trimmed for MEME-chip, e.g., `100`. This allows MEME-chip to examine the central region of the peaks for motifs while comparing to the flanking regions as a control for local sequence content. **!!TODO** please check this
    * `adaptor_SE`: path to file (FASTA format) containing sequencing adaptors for single-end (SE) experiments, *e.g.*, `TruSeq3-SE.fa`
    * `adaptor_PE`: path to file (FASTA format) containing sequencing adaptors for paired-end experiments, *e.g.*, `TruSeq3-PE-2.fa`
    * `NUM_TOP_PEAKS`: a value, the number of top peaks to consider, *e.g.*, 500
    * `MFOLD_MIN`: a value, the minimum fold depth enrichment required to call a peak for MACS2, *e.g.*, `5`. Not employed in our final analyses.
    * `MFOLD_MAX`: a value, the maximum fold depth enrichment allowed to call a peak for MACS2, *e.g.*, `50`. Not employed in our final analyses.
    * `blacklist`: path to file (BED format) containing a blacklist (excluded genome regions, including repetitive and low-complexity regions), *e.g.*, `ENCFF023CZC_sorted.bed`. See Amemiya *et al.* (2008).
    * `GENOME_FASTA`: path to file (FASTA format) containing the primary genome assembly for the organism of interest, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly.fa`.
    * `GENOME_IDX_PREFIX`: path to file (genome index) created using BOWTIE2, *e.g.*, `Homo_sapiens.GRCh38.dna.primary_assembly`. This can be accomplished using the following command: `bowtie2-build -f Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly`.

    The pipeline is shown below. Note that these tools use slightly different commands for single-end (SE) and paired-end (PE) reads (see script source code).

    ![Perl pipeline](images/perl_pipeline.png)

3. `PWM_pipeline.pl`. IDR and PCC. **!!TODO** Chase will add description

## Acknowledgments

We thank Frank Hsu, Federico Giorgi, Naoki Osato, Jeff Vierstra, and Meng Wang for helpful comments. This study was supported by Academia Sinica (AS-SUMMIT-109) and by Ministry of Science and Technology Taiwan (MOST 107-2311-B-001-016-MY3).

## Citation

When using this software, please refer to and cite:

> Yu CP, Kuo CH, Nelson CW, *et al.* 2021. <a target="_blank" href="https://www.pnas.org/content/118/20/e2026754118">Discovering unknown human and mouse transcription factor binding sites and their characteristics from ChIP-seq data</a>. *Proceedings of the National Academy of Sciences* **118**(20).

and this page:

> <https://github.com/chpngyu/chip-seq-pipeline/>

## References

* Amemiya HM, Kundaje A, Boyle AP. 2019. <a target="_blank" href="https://www.nature.com/articles/s41598-019-45839-z">The ENCODE Blacklist: Identification of Problematic Regions of the Genome</a>. *Scientific Reports* **9**: 9354.
* Bailey TL, Boden M, Buske FA, *et al.* 2009. <a target="_blank" href="https://academic.oup.com/nar/article/37/suppl_2/W202/1135092">MEME SUITE: tools for motif discovery and searching</a>. *Nucleic Acids Research* **37**: W202-W208.
* Bolger AM, Lohse M, Usadel B. 2014. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/30/15/2114/2390096">Trimmomatic: a flexible trimmer for Illumina sequence data</a>. *Bioinformatics* **30**(15): 2114-2120.
* Langmead B, Salzberg S. Langmead B, Salzberg SL. 2012. <a target="_blank" href="https://www.nature.com/articles/nmeth.1923">Fast gapped-read alignment with Bowtie 2</a>. *Nature Methods* **9**(4): 357-359.
* Li H, Handsaker B, Wysoker A, *et al.* 2009. <a target="_blank" href="https://pubmed.ncbi.nlm.nih.gov/19505943/">The Sequence Alignment/Map format and SAMtools</a>. *Bioinformatics* **25**(16): 2078-2079.
* Li Q, Brown JB, Huang H, Bickel PJ. 2011. <a target="_blank" href="https://projecteuclid.org/euclid.aoas/1318514284">Measuring reproducibility of high-throughput experiments</a>. *Annals of Applied Statistics* **5**(3): 1752--1779.
* Quinlan AR, Hall IM. 2010. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/26/6/841/244688">BEDTools: a flexible suite of utilities for comparing genomic features</a>. *Bioinformatics* **26**(6): 841-842.
* Xu M, Su Z. 2010. <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808352/">A Novel Alignment-Free Method for Comparing Transcription Factor Binding Site Motifs</a>. *PLoS ONE* **5**:e8797.
* Zhang Y, Liu T, Meyer CA, *et al.* 2008. <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137">Model-based analysis of ChIP-Seq (MACS)</a>. *Genome Biology* **9**(9): R137.
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
