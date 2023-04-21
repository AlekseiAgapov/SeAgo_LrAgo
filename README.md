# SeAgo & LrAgo

## What is this project for?
This is the code that was used for NGS analysis published in paper "Bacterial Argonaute proteins aid cell division in the presence of topoisomerase inhibitors in *Escherichia coli*" 
We sequenced small DNA molecules copurified with prokaryotic Argonaute proteins from *E. coli* cells.
This repository contains scripts that align NGS reads to the reference DNA and prepare data for visualization.

## How to use it?


Script reads_preprocessing.sh starts reads quality control with FastQC, then removes adaptors with cutadapt and makes quality control once again. For input one should specify (with `-d` argument) a path to a directory that contains single fq.gz file with reads. Output - directories with FastQC output and trimmed.fastq.gz file with the processed reads. **It is important to check the quality of the processed reads and manually adjust cutadapt arguments if needed.**

Bash script pipeline.sh in the <pipeline> directory launches the code. The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should include trimmed.fastq.gz file (result of reads_preproessing.sh script) and a directory that contains three fasta files: genome.fa (the header should be >genome), plasmid.fa (the header should be >plasmid), and genome.gff3 (note that the name of the chromosome should be "genome").
 
The script will calculate alignments and create a set of directories in the working directory, where R files for plot drawing will be copied. Manual adjusting of parameters (such as axis limits and labels) in R scripts is required.
The result will be the following plots:
- genome coverage with small DNA reads
- metaplot of small DNA reads coverage of regions around Chi-sites
- aligned reads logo
- GC-content along the guide length and in surrounding sequences of chromosomal DNA
- length distribution of the aligned reads
- coverage of genes and intergenic regions with smDNA


Directory <revision> contains scripts written by Dmitry Sutormin at the revision stage. These scripts are used to produce the following plots:
- metaplot of small DNA reads coverage of regions around Chi-sites
- metaplot of small DNA reads coverage of regions around gyrase cleavage sites
- metaplots of small DNA reads coverage of genes and intergenic regions


## Requirements
This script utilizes some commonly used programs for data analysis and NGS analysis:
- pandas library for Python https://github.com/pandas-dev/pandas
- FastQC https://github.com/s-andrews/FastQC
- cutadapt https://github.com/marcelm/cutadapt
- bowtie https://github.com/BenLangmead/bowtie
- samtools https://github.com/samtools/samtools
- bedtools https://github.com/arq5x/bedtools2
- ggplot2 library for R https://github.com/tidyverse/ggplot2
- ggseqlogo library for R https://github.com/omarwagih/ggseqlogo
- blast+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

 ## Where to find the data that was processed with this code?
 The raw sequencing reads are deposited in SRA in the BioProject PRJNA878808.
