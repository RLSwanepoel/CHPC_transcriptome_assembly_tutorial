#!/bin/bash
## env gsnapdir="gsodm13xpa gsodm13xtr" datad=`pwd` qsub -q normal samindxn.sh
#... -N sams
#... -A ind114
#... -l nodes=2:ppn=32,walltime=3:55:00
#... -o sams.$$.out
#... -e sams.$$.err
#... -V

## need to resort bams to index/count; now are all -n namesorted; 
## redo, problems .. takes >4hr, some fails mem?
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
    nam=`echo $ibam | sed "s/$bigone.bam//;"`
    obam=`echo $ibam | sed "s/.bam/s/;"`
    if test -f $obam.count ; then  echo "HAVE $obam.count"; continue;  fi

    (echo fork1 samtools sort -m $maxm $ibam $obam..samtools index $obam..samtools idxstats $obam T $obam.count; ) &
    i=$(( $i + 1 )); 

    bsmall=`ls $nam*.bam | grep -v $bigone`
    obam="$nam${smallone}s"
    if test -f $obam.count ; then  
      echo "HAVE $obam.count"; 
    else 
      ( for sbam in $bsmall; do { obam=`echo $sbam | sed 's/.bam/s/'`; \
 echo fork2 samtools sort -m $maxm $sbam $obam..samtools index $obam..samtools idxstats $obam T $obam.count; } done ) &
      i=$(( $i + 1 )); 
    fi

    if [ $i -ge $ncpu ]; then wait; i=0; fi
    echo icpu $i
  } done
} done

echo waiting
wait
echo done wait
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

