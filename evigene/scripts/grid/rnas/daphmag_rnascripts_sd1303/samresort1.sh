#!/bin/bash
## env bam=xxx datad=`pwd` qsub -q shared samsort1.sh
#PBS -N sams
#PBS -A ind114
### dang takes longer than 6hr;
#PBS -l nodes=1:ppn=2,walltime=15:55:00
#PBS -o sams.$$.out
#PBS -e sams.$$.err
#PBS -V

ncpu=2
# maxmem=16,000,000,000
maxmem=1500000000
## OPT: sort -n : namesort

samd=$HOME/bio/bin
export PATH=$PATH:$samd

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/fasta"; exit -1; fi
if [ "X" = "X$bam" ]; then echo "ERR: missing env bam=mydata.bam"; exit -1; fi

cd $datad
inam=`echo $bam | sed 's/\.bam//; s/\.merge//; s/\.sort//;'`

$samd/samtools sort -n -m $maxmem $bam $inam.nsort
# Usage: samtools sort [-on] [-m <maxMem>] <in.bam> <out.prefix>

