#!/bin/bash
## env datad=`pwd` qsub -q normal samsortn.sh
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
export PATH=$PATH:$samd
# trestles
# datad=$HOME/scratcht/chrs/daphmag/rnas/

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/sam"; exit -1; fi
# if [ "X" = "X$isamu" ]; then echo "ERR: missing env isamu=mydata.sam"; exit -1; fi
#?? is trsize needed as .sam header has tr sizes??  replace -t trsize w/ -S
# if [ "X" = "X$trsize" ]; then echo "ERR: missing env trsize=mydata.tr.count"; exit -1; fi

cd $datad
samset=`ls *.sam`

i=0;
for isamu in $samset; do {
  inam=`echo $isamu | sed 's/.sam//'`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  if test -f $inam.sort.bam ; then  echo "done $inam.sort.bam"; continue;  fi

  ( samtools view -u -S $isamu | samtools sort -m $maxmem - $inam.sort; \
    samtools index $inam.sort.bam; samtools idxstats $inam.sort.bam > $inam.count; ) &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
  
} done
wait

##?need trsize?:samtools view -u -t $trsize $isamu | samtools sort - $inam.sort;

