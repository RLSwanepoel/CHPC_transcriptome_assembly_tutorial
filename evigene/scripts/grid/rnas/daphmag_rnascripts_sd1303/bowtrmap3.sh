#! /bin/bash
### env trdb=xxx readset=fqdir/xxx*_1.fastq.gz qsub -q normal bowtrmap2.sh
#PBS -N bowtrmap
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=11:55:00
#PBS -o bowtrmap.$$.out
#PBS -e bowtrmap.$$.err
#PBS -V

# ncpu=32; 
#badder#ncpuPerMap=$ncpu
## try 1 readset per  : runs in 47Gb for 1 smallest readset, 12 cpu
ncpu=12
ncpuPerMap=12

## .. also try replace --all w/ -k 99 or 199 ..
## .. try bow2build --offrate 3 : split diff for dbmem use. and use --mm memmap
## .. still mempig: mem = 63706816kb vmem = 68316776kb  after 30min
## .. cut down to 1 readset/16cpu, or 2 readset/8cpu each ?; test on gordon 16cpu?
## .. is --all the culprit? or --offrate < 5 ?
## ** FAIL outamem 32cpu/map for --all;  increase --offrate ??
## ** Too much mem use: vmem = 69837796kb  for --all 
## .. ? test allcpuPerMap > BADDER, test --mm  memmapdb ?
## .. mem hogging could also be perl wrapper + input fastq.gz .. try unzip fastq
## 8cpu x 21h not enough .. 16cpu per file : good, 8hr for pair.
## .. switch to -x bowxevg/transcripts.tr --offrate 1 ; skip --all for now
## ** need --all to get dup mappings, very common over alt tr..

bowd=$HOME/bio/bowtie2
samd=$HOME/bio/bin
export PATH=$PATH:$bowd:$samd
# trestles
datad=$HOME/scratcht
workd=$datad/chrs/daphmag/rnas/

##* change to env readset=fq4x2/Dman_7[45]_*_1.fastq.gz == file list
## readset=(fq4x2/Dman_7[45]_*_1.fastq.gz)  ; rd1=${readset[0]}
## or readset="fq4x2/Dman_79_*_1.fastq.gz fq4x2/Dman_80_*_1.fastq.gz"
tnam=`basename $trdb .tr`
reada=($readset)
rnam=`dirname ${reada[0]}`
odir=bowta$rnam
outna=$odir/$tnam-$rnam

#last# opts="-q --threads $ncpuPerMap --all -X 800 --end-to-end --fast "
#try# 
opts="-q --threads $ncpuPerMap -k 199 -X 800 --end-to-end --fast "

## NO HELP: --mm : memmapdb save mem for 2+? ; --offrate 2 save mem?
## new: dont use --offrate here, uses bow.index offrate #ERR* bowtie2-align: unrecognized option `--offrate'
## -q == fastq : default; --end-to-end : default;  -M 5 : default
## --un-conc <path>      write pairs that didn't align concordantly

cd $workd/
mkdir $odir
echo "START bowtie2 : `date`"

i=$ncpuPerMap;  # start here so next i > ncpu waits
for lreads in $readset; do {
  rreads=`echo $lreads | sed 's/_1\./_2./; '`
  rdnam=`basename $lreads .fastq.gz | sed 's/\.fastq.*//; s/\.fq.*//; s/_1\..*//; '`
  if [ ! -f $rreads ]; then 
    echo "ERR: missing right _2.fq of $lreads"; continue;
  fi 

  echo bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog
  bowtie2 $opts -x $trdb -1 $lreads -2 $rreads \
     --un-conc $odir/$rdnam.unc.fq -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog  &

  i=$(( $i + $ncpuPerMap ))
  if [ $i -gt $ncpu ]; then wait; i=$ncpuPerMap; fi

} done
wait

echo "DONE bowtie2 : `date`"

