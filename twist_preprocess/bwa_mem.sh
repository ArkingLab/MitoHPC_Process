#!/bin/bash
#SBATCH --job-name=bwamap
#SBATCH --cpus-per-task=1
#SBATCH --time 120:00:00

#--nodelist=compute-09*,compute-1*

#module load trimgalore
#module load bwa
module load samtools
#module load python

IN1=$1
IN2=$2
OUT=$3

echo $OUT

/dcs05/legacy-dcl01-arking/data/active/newseq/MitoHPC/bin/bwa mem /dcs04/legacy-dcs01-arking/arkinglab/software/shared/genomes/hg38.bwa/GRCh38 \
$IN1 $IN2 | \
samtools sort -o ./mapped/$OUT.bam -
