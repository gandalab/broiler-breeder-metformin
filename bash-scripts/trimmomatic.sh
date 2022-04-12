#!/bin/bash

#PBS -N trim
#PBS -A evk5387_d_g_hc_default
#PBS -l nodes=1:ppn=1
#PBS -l pmem=4gb
#PBS -j oe
#PBS -l walltime=2:00:00

## run trimmomatic for length & quality
conda activate bioinfo
set -uex

# working directory
WORKDIR=/gpfs/group/evk5387/default/Novogene/metformin

# input directory
INDIR=/gpfs/group/evk5387/default/Novogene/adapter-removed/

# output directory
OUTDIR=/gpfs/group/evk5387/default/Novogene/metformin/trimmed/
mkdir -p $OUTDIR

# run trimmomatic
cat $WORKDIR/ids.txt | parallel /storage/work/epb5360/miniconda3/envs/bioinfo/bin/trimmomatic PE $INDIR/{}.trimmed_1.fastq $INDIR/{}.trimmed_2.fastq -baseout $OUTDIR/{}.trim.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100

# print finished message
echo "trimming complete"
