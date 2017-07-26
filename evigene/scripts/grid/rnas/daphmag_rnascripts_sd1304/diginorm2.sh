#! /bin/bash
### env fadir=fa4x fain1=xxx.fa2.gz fain2=yyy.fa2.gz qsub -q normal dignorm.sh
#PBS -N dignorm
#PBS -A ind114
#... -l mem=40gb,nodes=1:ppn=2,walltime=39:55:00
#PBS -l nodes=1:ppn=32,walltime=29:55:00
#PBS -o dignorm.$$.out
#PBS -e dignorm.$$.err
#PBS -V

##use env : fain=fasta5/all_nodam.fa2; 
## FIXME: normapp doesnt preserve fain path, puts keep in cur wd
## diginorm2: run this on *each* fastadir/Dman*.fa2.gz ; all_*.fa2 getting huge (>1TB)
##  .. run as normapp filelist, loop here, or qsub per file?
##  .. run 2 in 64gb node ?

kopt="-p";

# mem:  8e9 x 4 = 32Gb ; 10e9 x 4 = 40Gb; 64gb = 16e9
#62gb# xhash=15.6e9
#40gb# xhash=9.8e9
xhash=7.9e9
kmer=30
rkeep=8
ncpu=2

#** requires python2.6+ not on gordon.sdsc, on trestles module add python
module() { eval `/opt/modules/Modules/3.2.5/bin/modulecmd bash $*`; }
module add python

export PYTHONPATH=$HOME/bio/khmer/python
normapp=$HOME/bio/khmer/scripts/normalize-by-median.py
datad=$HOME/scratchn
workd=$datad/chrs/daphmag/rnas

cd $workd
# fadir=`dirname $fain`
# fain=`basename $fain`
cd $fadir

if [ -f $fain1 ]; then
 $normapp $kopt -C $rkeep -x $xhash  -k $kmer  $fain1 &
fi
if [ -f $fain2 ]; then
 $normapp $kopt -C $rkeep -x $xhash  -k $kmer  $fain2 &
fi

wait

# fao=`echo $fain | sed 's/\.fasta//; s/\.fa//;'`
#  mv $fain.keep $fao.dnorm$kmer.fa

