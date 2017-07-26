#!/bin/bash
## env gsnapdir="gsodm13xpa gsodm13xtr" datad=`pwd` qsub -q normal samindxn.sh
#PBS -N sams
#PBS -A ind114
#PBS -l nodes=2:ppn=16,walltime=29:55:00
#PBS -o sams.$$.out
#PBS -e sams.$$.err
#PBS -V

## dang outatime at 10hr for concordant_mults.bam sorts .. big files

ncpu=32
maxm=2000000000
samd=$HOME/bio/bin
export PATH=$PATH:$samd

if [ "X" = "X$gsnapdir" ]; then echo "ERR: miss gsnapdir=gsnap dir list "; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapdirs"; exit -1; fi

## gsnapdir="gsodm13xpa gsodm13xtr" ..
## for concordant_mult.bam only big one, split loop in gdir to fork 1 for that, others in 2nd fork
bigone=concordant_mult
keepset="concordant_uniq halfmapping_uniq paired_mult paired_uniq_long"
poorset="halfmapping_mult paired_uniq_inv paired_uniq_scr unpaired_mult unpaired_uniq nomapping"
smallkeep="$keepset $poorset"
smallone=unpaired_uniq

cd $datad/

i=0; 
for gdir in $gsnapdir ; do {
  echo PROCESS in $datad/$gdir;
  cd $datad/$gdir;
  blist=`ls *$bigone.bam`;
  for ibam in $blist; do {
    obam=`echo $ibam | sed "s/.bam/s/;"`; obb=$obam.bam;
    nam=`echo $ibam | sed "s/$bigone.bam//;"`
    if test -f $obam.count ; then  echo "HAVE $obam.count"; continue;  fi

    ( samtools sort -m $maxm $ibam $obam; samtools index $obb; samtools idxstats $obb > $obam.count; ) &
    i=$(( $i + 1 )); 

    bsmall=`ls $nam*.bam | grep -v $bigone`
    obam="$nam${smallone}s"
    if test -f $obam.count ; then  
      echo "HAVE $obam.count"; 
    else 
      ( for sbam in $bsmall; do { obam=`echo $sbam | sed 's/.bam/s/'`; obb=$obam.bam; \
       samtools sort -m $maxm $sbam $obam; samtools index $obb; samtools idxstats $obb > $obam.count; } done ) &
      i=$(( $i + 1 )); 
    fi

    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done
} done

wait

#.......
# i=0;
# for gdir in $gsnapdir ; do { 
#   echo PROCESS in $datad/$gdir; 
#   cd $datad/$gdir; 
#   blist=`ls *.bam`;
#   for ibam in $blist; do { 
#     if test -f $ibam.count ; then  echo "done $ibam.count"; continue;  fi
#     ( samtools index $ibam; samtools idxstats $ibam > $ibam.count; ) &
#     i=$(( $i + 1 )); if [ $i -ge $ncpu ]; then wait; i=0; fi
#   } done
# } done
# 
# wait
#..............

