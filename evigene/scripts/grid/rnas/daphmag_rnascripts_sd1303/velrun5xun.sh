#!/bin/bash
##  env qsub -q normal velrun1.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion #below: dv=2x
ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

## bin1 for simple, kmer<60 .. bin2/ for kmer<106
velbin=$HOME/bio/velvet127s/bin2
velbin1=$HOME/bio/velvet127s/bin
## trestles
datad=$HOME/scratcht
workd=$datad/chrs/daphmag/rnas

## velrun5nun.sh == reads unmapped to 1st trasm set; groups/clone combined
# dv=5nun1; inpe="$workd/fq5nun/*.unc.[12].fq"
## 5nun1 13G fail outamem at k49
## inner insize for fq5n: 150 ; outer insize = 200+150
## sam10 n=15035374; ains=164; mdins=152; sam4 n=15296969; ains=174; mdins=159
##....
##.. obs ins sizes hiskx: 200..250; add 200 read pair = 425
# n=1519766; ains=288; mdins=265; n=1482344; ains=283; mdins=261; n=1456546; ains=224; mdins=203
# dv=5xun1; inpe="$workd/soap5xun1/hiskx2l47_unc.[12].fq"
dv=5xun2; inpe="$workd/soap5xun2/hiskx2l36_unc.[12].fq"

rund=velv$dv

## boost minpair for kmer < 45 ?
vopts="-ins_length 425 "
ooptklo="-min_pair_count 4 -scaffolding yes -min_trans_lgth 180 $vopts"
oopts="-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 $vopts"
kset="87 81 75 65 59 55 49 35 31 "

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
  vbin=$velbin
  ooptk=$oopts
  if [ $k -lt 57 ]; then
   vbin=$velbin1
   ooptk=$ooptklo
  fi
  $vbin/velveth $ksubdir $k  -reuse_Sequences
  $vbin/velvetg $ksubdir $vopts -read_trkg yes 
  $vbin/oases   $ksubdir $ooptk
  
  /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
  echo "#.. end velrun $ksubdir : `date`"
}
done

echo "DONE `date` " 

