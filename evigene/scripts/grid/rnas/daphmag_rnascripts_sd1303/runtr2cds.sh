#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=29:55:00
#PBS -o tr2cds.$$.out
#PBS -e tr2cds.$$.err
#PBS -V

ncpu=32; maxmem=54000
## bump up walltime: 18>24hr; daphmag5iall9 blastn takes 12hr; big trset
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

## rerun time ~1hr; blastn -ungapped reduced main set 70k to 62k; keep pi=98 w/ ungapped..
## MINCDS=90 CDSBLAST_IDENT=98; CDSBLAST_EVALUE=1e-19; # is this ok?
opts="-MINCDS=120 -CDSBLAST_IDENT=98"

evigene=$HOME/bio/evigene/scripts
## NEW cd-hit version cures 32kseq-bug; cdhit461
export PATH=$HOME/bio/cdhit461/bin:$PATH
# fastanrdb
export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
# blastn:
export PATH=$HOME/bio/ncbi2227/bin:$PATH

if [ "X" = "X$trset" ]; then
  echo "missing env trset=xxxx.tr"; exit -1
fi
if [ "X" = "X$datad" ]; then
  echo "missing env datad=path/to/data"; exit -1
fi

cd $datad/
echo $evigene/prot/tr2aacds.pl $opts -debug -NCPU $ncpu -MAXMEM $maxmem -tidy -log -cdna $trset
$evigene/prot/tr2aacds.pl $opts -debug -NCPU $ncpu -MAXMEM $maxmem -tidy -log -cdna $trset

