#!/bin/bash
##  env qsub -q normal velrun1.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion #below: dv=2x

ncpu=14
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

## trestles
datad=$HOME/scratcht
## bin1 for simple, kmer<60 .. bin2/ for kmer<106
velbin=$HOME/bio/velvet127s/bin2
workd=$datad/chrs/daphmag/rnas

## velrun5xun.sh == reads unmapped to 1st trasm set; groups/clone combined
# dv=4xun1; inpe="$workd/fq4xun/*.unc.[12].fq"
dv=5nun1; inpe="$workd/fq5nun/*.unc.[12].fq"
## 5nun1 13G fail outamem at k49
## inner insize for fq5n: 150 ; outer insize = 200+150
## sam10 n=15035374; ains=164; mdins=152; sam4 n=15296969; ains=174; mdins=159

rund=velv$dv

## maybe should boost minpair for kmer < 45 ?
vopts="-ins_length 350 "
ooptklo="-min_pair_count 4 -scaffolding yes -min_trans_lgth 180 -ins_length 350"
oopts="-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 -ins_length 350"
# kset0="95 91 85 81 75 65 55 45 35 29 25"
kset="81 75 65 59 55 49 45 35 31"

echo "START `date` " 
cd $workd
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
   -shortPaired -fastq  -separate $inpe
fi

#  -fasta.gz -shortPaired -interleaved $inpe 
#  -shortPaired -fastq.gz -separate $velfad/SRR346404_[12].fastq.gz 

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  ooptk=$oopts
  if [ $k -lt 50 ]; then
   ooptk=$ooptklo
  fi
  $velbin/velveth $ksubdir $k  -reuse_Sequences
  $velbin/velvetg $ksubdir $vopts -read_trkg yes 
  $velbin/oases   $ksubdir $ooptk
  
  /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
  echo "#.. end velrun $ksubdir : `date`"
}
done

echo "DONE `date` " 

