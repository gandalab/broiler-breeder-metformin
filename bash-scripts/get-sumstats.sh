#!/bin/bash

#PBS -N sumst
#PBS -l walltime=3:00:00
#PBS -A evk5387_d_g_hc_default
#PBS -l nodes=1:ppn=1
#PBS -l pmem=8gb
#PBS -j oe

# get summary stats on trimmed reads
conda activate bioinfo
set -uex

# set wd
WORKDIR=/gpfs/group/evk5387/default/Novogene/metformin

# TRIMMED READS
seqkit stat $WORKDIR/trimmed/*_1P.fastq -j 4 -T -o $WORKDIR/trimstatf.txt
# reverse
seqkit stat $WORKDIR/trimmed/*_2P.fastq -j 4 -T -o $WORKDIR/trimstatr.txt
echo "trimmed done"
