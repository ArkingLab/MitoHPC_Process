#!/bin/bash

################### This script is for pre_processing fastq files from twist

#### Steps include trimming adapter, aligning with bwa-mem, and sub_sampling reads
#### The BAM file can then be used as input to MitoHPC
#### 05/14/2025

#### Create a list of all fastq files

## Total file counts
ls -1 *.fastq.gz | wc -l

## Get all unique ids (adjust field based on your file name)
ls -1 *.fastq.gz | awk -F '_' '{printf "%s\n",$2}' | uniq > id.txt


#### FASTQ QC

## Check read length (150 or 250) for fastq files
for file in *.fastq.gz; do zcat "$file" | head -n 2 | tail -n 1 | wc -c; done | sort | uniq -c

## Check all read length for all file (not really necessary as this take longer times)
for file in *.fastq.gz; do zcat $file | awk 'NR % 4 == 2' | awk '{print length}' | sort | uniq -c ; done
zcat $file | awk 'NR % 4 == 2' | awk '{print length}' | sort | uniq -c

## Check file sizes before trimming adapter (files with sequencing error or just water/control samples)
ls -la | grep .fastq.gz$ | awk ' $5 < 2000000 {print $0}'
# remove water and samller ones manually
# update id


#### TRIM adapters (adjust file names accordingly)

mkdir trimmed
## Adjust the adapter sequence in the trim_adapter.sh script accordingly
cat id.txt | perl -ane 'print "sbatch --mem=3G ./trim_adapter.sh Arking_$F[0]\_L002_R1_001.fastq.gz Arking_$F[0]\_L002_R2_001.fastq.gz\n";' | bash

## size check
ls -la | grep .fq.gz$ | awk ' $5 < 2000000 {print $0}' | wc -l

## Check length distribution of a single fastq file
zcat sample.fq.gz | awk 'NR % 4 == 2' | awk '{print length}' | sort | uniq -c


#### Alignment using BWA-MEM
mkdir mapped
cat id.txt | perl -ane 'print "sbatch --mem=20G ./bwa_mem.sh trimmed/Arking_$F[0]\_L002_R1_001_val_1.fq.gz trimmed/Arking_$F[0]\_L002_R2_001_val_2.fq.gz $F[0]\n";' | bash

## QC steps (remove failed jobs)
ls | grep tmp | perl -ane 'print "rm $F[0]\n";' | bash


#### BAM QC

## Original bam start position check
samtools view Sample-1_S12.bam | awk '$3 ~ /chrM/ {print $0}' | cut -f 4 | sort | uniq -c

## Original bam check file mapped read length (see if broken down after trimming adapters)
samtools view Sample-1_S12.bam | cut -f 10 | awk '{print length}' | sort | uniq -c

## check file size & delete if necessary
ls -la | grep .bam$ | awk ' $5 < 1000000 {print $0}' | wc -l
#delete command
ls -la | grep .bam$ | awk ' $5 < 1000000 {print $9}' | awk -F . '{printf "rm %s*\n", $1}' | bash


#### Sub-sample for even cvg across each file/sample (for samples with much higher coverage ~1G in file size)

## The downsample script can be found in the scripts/ directory in the MitoHPC repository
mkdir mapped_subs
find mapped/ -name "*.bam" | perl -ane '/.+\/(.+).bam/; \
print "sbatch --mem=5G MitoHPC/scripts/downsampleSam.sh mapped/$1 mapped_subs/$1 20 chrM\n";' | bash


#### Check sub_sample bam files size

ls -lah | grep .bam$ | awk '{print $5}' | sort | uniq -c
ls -la | grep .bam$ | awk ' $5 < 3000000 {print $0}' | wc -l
# get the number of bams less then 3 MB (14, here all 14 are water and can be deleted)
ls -la | grep .bam$ | awk ' $5 < 3000000 {print $9}' | awk -F . '{printf "rm %s*\n", $1}' | bash


#### Run MitoHPC (more details on GitHub page)


