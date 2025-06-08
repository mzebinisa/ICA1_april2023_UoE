
Bacterial RNA-seq Analysis Pipeline

This is a Bash script I wrote for processing bacterial RNA-seq data—from raw FASTQ files all the way to group-wise expression summaries. It’s meant to make the workflow easier and more reproducible by automating each major step: quality control, alignment, counting, and basic stats. 
This was part of an assignment in BPSM(Bioinformatics programming and systems management module) at UOE, and I wanted a clear, reproducible way to go from raw sequencing files to basic gene-level summaries. I’ve added some prompts and messages  so that the user potentially can follow what’s going on. 

What it does

Prompts the user for input/output directories, filenames, and references interactively
Copies only the required FASTQ files to the current working directory
Runs `FastQC` on raw reads and summarizes results with `MultiQC`
Builds a Bowtie2 index from a reference genome and aligns all paired-end reads
Converts SAM to sorted and indexed BAM files using `samtools`
Uses `bedtools` to:
 
  Identify overlaps between reads and annotation features
  Generate counts per region (via `multicov`)

Organizes counts by experimental groups (based on a `Tco.fqfiles` metadata file)
Aggregates and computes group-wise mean expression levels

Which Tools were used

fastqc, multiqc

bowtie2

samtools

bedtools

Bash (was mandatory at this point)


Input formats

Paired-end .fq.gz files with names like Tco-XXXX_1.fq.gz and Tco-XXXX_2.fq.gz
Reference genome in .fasta.gz format
BED file for annotated genomic regions
A Tco.fqfiles metadata file containing sample info (SampleName, SampleType, Replicate, etc.)
