#!/bin/bash
#SBATCH --job-name=trimadp
#SBATCH --cpus-per-task=1

module load trimgalore
#module load bwa
module load samtools
#module load python

IN1=$1
IN2=$2

echo $IN1

## modify the Adapter sequences accordingly

trim_galore \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ./trimmed/ \
--no_report_file \
--suppress_warn \
--paired $IN1 $IN2
