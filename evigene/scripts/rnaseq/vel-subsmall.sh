#!/bin/bash
##  qsub -q batch velsubt.sh
#PBS -N velsmal
#PBS -l mem=256gb,nodes=1:ppn=32,walltime=23:55:00
#PBS -m abe
#PBS -j oe
#PBS -k o
#PBS -M gilbertd@indiana.edu

## note velv OPENMP threads trades w/ memory, more mem w/ more cpu; also uses +1 threads
## note2: oases doesnt use threads; not much gain, use 1 core/subset or 2? : 10 subs/part
## factor in memory total = 10x subset; how many cores/node??
ncpu=3
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

workd=$HOME/scratch/chrs/nasv1

cd $workd/rnas/velsmal/

echo "#.. start subsetvel : `date`"

export subset=subset.small0 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small1 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small2 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small3 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small4 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small5 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small6 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small7 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small8 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.small9 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &

wait
echo "#.. end subsetvel : `date`"

