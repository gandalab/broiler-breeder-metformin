#!/bin/bash

#PBS -N cutad
#PBS -A open
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=10
#PBS -l pmem=8gb
#PBS -j oe

## basic bash file for cutadapt on all samples

conda activate bioinfo
set -uex

# input directory
INDIR=/gpfs/group/evk5387/default/Novogene/all-samples/

# working directory for sample names
WORKDIR=/gpfs/group/evk5387/default/Novogene/
mkdir -p $WORKDIR

# output directory
OUTDIR=/gpfs/group/evk5387/default/Novogene/adapter-removed/
mkdir -p $OUTDIR

# fastqc directory
QUALDIR=/gpfs/group/evk5387/default/Novogene/quality
mkdir -p $QUALDIR

# adapter sequences for cutadapt
AD1=GTGCCAGCMGCCGCGGTAA
AD2=GGACTACHVGGGTWTCTAAT

# get list of sample names
ls $INDIR | awk -F '_[[:digit:]]' '{print $1}' | grep -v raw | uniq  > $WORKDIR/names.txt

## ---- cutadapt ----

# run in parallel
#cat $WORKDIR/names.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/cutadapt -j 0 -a $AD1 -A $AD2 --max-n=0 -o $OUTDIR/{}.trimmed_1.fastq -p $OUTDIR/{}.trimmed_2.fastq $INDIR/{}_1.fq.gz $INDIR/{}_2.fq.gz

# print finished message
echo "cutadapt done"

## ---- fastqc on trimmed reads ----

# run in parallel - forward reads
cat $WORKDIR/names.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastqc $OUTDIR/{}.trimmed_1.fastq  -o $QUALDIR/

# run in parallel - reverse reads
cat $WORKDIR/names.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastqc $OUTDIR/{}.trimmed_2.fastq  -o $QUALDIR/

# get multiqc report for all fastqc
multiqc $QUALDIR/* -o $WORKDIR -n adapter-removed.multiqc

# print finished message
echo "quality check done"
