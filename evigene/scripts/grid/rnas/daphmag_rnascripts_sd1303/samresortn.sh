#!/bin/bash
## env bamset=xxx*.bam datad=`pwd` qsub -q normal samresortn.sh
#PBS -N sams
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=3:55:00
#PBS -o sams.$$.out
#PBS -e sams.$$.err
#PBS -V

ncpu=32
# samsort def 500M maxmem=500000000
maxmem=1500000000

samd=$HOME/bio/bin
# export PATH=$PATH:$samd

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/sam"; exit -1; fi
if [ "X" = "X$bamset" ]; then echo "ERR: missing env bamset=all\*.bam"; exit -1; fi

cd $datad
i=0;
for bam in $bamset; do {
  inam=`echo $bam | sed 's/\.bam//; s/\.sort//;'`
  if ! test -f $bam ; then  echo "missing $bam"; continue;  fi
  if test -f $inam.nsort.bam ; then  echo "done $inam.nsort.bam"; continue;  fi

  $samd/samtools sort -n -m $maxmem $bam $inam.nsort &
  # Usage: samtools sort [-on] [-m <maxMem>] <in.bam> <out.prefix>

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
  
} done
wait

