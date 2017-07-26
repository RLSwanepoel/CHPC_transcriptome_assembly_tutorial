#! /bin/bash 
### env protin=prot.aa qsub -q batch cdhitprot.sh
#PBS -N cdhitp
#PBS -l nodes=1:ppn=16,walltime=12:55:00
#PBS -o cdhit.$$.out
#PBS -e cdhit.$$.err
#PBS -V

ncpu=16
export LANG=C
bindir=/N/soft/mason/cd-hit/cd-hit-v4.5.6-gcc-4.4.5-openMP

## note2: protin must be cleaned of '\*$' stopcodons
## note: -c 0.50 requires -n 3  wordlen
## 0.75: cut 650k to 500k clusters, many tinyaa inside big at 75% id
# cdopt="-c 0.75 -T $ncpu -M 32000 -d 0" ;  bltag="cd75"
# cdopt="-c 0.51 -n 3 -T $ncpu -M 32000 -d 0";  bltag="cd51" # takes >6hr
cdopt="-c 0.60 -n 4 -T $ncpu -M 32000 -d 0";  bltag="cd60"

workd=$HOME/scratch/chrs/aabugs4
cd $workd/prot/

echo cd-hit $cdopt -i $protin -o $protin.$bltag 
$bindir/cd-hit $cdopt -i $protin -o $protin.$bltag 

