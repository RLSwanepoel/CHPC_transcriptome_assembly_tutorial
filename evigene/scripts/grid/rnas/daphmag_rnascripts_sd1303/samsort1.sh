#!/bin/bash
## env isamu=xxx datad=`pwd` qsub -q shared samsort1.sh
#PBS -N sams
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=15:55:00
#PBS -o sams.$$.out
#PBS -e sams.$$.err
#PBS -V

ncpu=2
namesort=1
maxmem=10000000000
## OPT: sort -n : namesort

samd=$HOME/bio/bin
export PATH=$PATH:$samd
# trestles
# datad=$HOME/scratcht/chrs/daphmag/rnas/

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/fasta"; exit -1; fi
if [ "X" = "X$isamu" ]; then echo "ERR: missing env isamu=mydata.sam"; exit -1; fi
## NO need ?? is trsize needed as .sam header has tr sizes??  replace -t trsize w/ -S
# if [ "X" = "X$trsize" ]; then echo "ERR: missing env trsize=mydata.tr.count"; exit -1; fi

cd $datad
inam=`echo $isamu | sed 's/.sam//'`

if [ $namesort = 1 ]; then
  samtools view -u -S $isamu | samtools sort -n -m $maxmem - $inam.nsort; 

else
  samtools view -u -S $isamu | samtools sort - $inam.sort; \
  samtools index $inam.sort.bam; samtools idxstats $inam.sort.bam > $inam.count; 
fi

##NO ?need trsize?:samtools view -u -t $trsize $isamu | samtools sort - $inam.sort;

