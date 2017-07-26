#!/bin/bash
## env sort=name bamset=xxx outna=xxx datad=`pwd` qsub -q shared samsort1.sh
#PBS -N sams
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=15:55:00
#PBS -o sams.$$.out
#PBS -e sams.$$.err
#PBS -V

ncpu=2
# maxmem=10000000000
# namesort=1; fixme: merge -n for namesort
opt=""
if [ "X" = "X$sort" ]; then 
  opt=""; 
else
  if [ $sort = "name" ]; then opt="-n"; fi
  if [ $sort = "n" ]; then opt="-n"; fi
fi

samd=$HOME/bio/bin
# export PATH=$PATH:$samd

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/fasta"; exit -1; fi
if [ "X" = "X$bamset" ]; then echo "ERR: missing env bamset=mydata.bam"; exit -1; fi
if [ "X" = "X$outna" ]; then echo "ERR: missing env outna=myoutput "; exit -1; fi

cd $datad

#  merge -1 = gzip --fast
$samd/samtools merge -1 $opt $outna.bam $bamset 
$samd/samtools flagstat $outna.bam > $outna.fstat

exit
