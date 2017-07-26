#!/bin/bash
##  env fqdir=xxx qsub -q normal velrun1.sh
#PBS -N makefa
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=4:55:00
#PBS -o makefa.$$.out
#PBS -e makefa.$$.err
#PBS -V

ncpu=16

datad=$HOME/scratchg/
# /home/ux455375/scratchg/chrs/daphmag/rnas
workd=$datad/chrs/daphmag/rnas/
evs=$HOME/bio/evigene/scripts

cd $workd

fqs=`ls $fqdir/*_1.fastq.gz`
# outdir=`echo $fqdir | sed 's/fastq/fasta/;'`
outdir=$fqdir

i=0; for fqin in $fqs ; do 
{
  echo $evs/rnaseq/fa2pairfa.pl -out $outdir -fq $fqin 
  $evs/rnaseq/fa2pairfa.pl -fq $fqin &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done
wait

fa2=`ls $outdir/*.fa2`
i=0; for fain in $fa2; do 
{
  gzip --fast $fain &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done
wait


