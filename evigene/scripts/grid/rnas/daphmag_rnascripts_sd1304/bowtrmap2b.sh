#! /bin/bash
### env trdb=xxx readset=fqdir/xxx*_1.fastq.gz qsub -q normal bowtrmap2.sh
#PBS -N bowtrmap
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=18:55:00
#PBS -o bowtrmap.$$.out
#PBS -e bowtrmap.$$.err
#PBS -V

ncpu=32; ncpuPerMap=16
## 8cpu x 21h not enough .. 16cpu per file : good, 8hr for pair.
## .. switch to -x bowxevg/transcripts.tr --offrate 1 ; skip --all for now

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
# rnam=`basename $readset _1.fastq.gz`
reada=($readset)
rnam=`dirname ${reada[0]}`
odir=bowtr$rnam
outna=$odir/$tnam-$rnam
#noneed# trsize=`echo $trdb | sed 's/$/.count/'`

opts="-q --threads $ncpuPerMap -X 800 --end-to-end --fast "
## new: dont use --offrate here, uses bow.index offrate #ERR* bowtie2-align: unrecognized option `--offrate'
## -q == fastq : default; --end-to-end : default;  -M 5 : default
## --un-conc <path>      write pairs that didn't align concordantly

cd $workd/
mkdir $odir

#OLD# readset=`ls $reads/*_1.fastq.gz`

i=0;
for lreads in $readset; do {
  rreads=`echo $lreads | sed 's/_1\./_2./; '`
  rdnam=`basename $lreads .fastq.gz | sed 's/_1\..*//; '`
  if [ ! -f $rreads ]; then 
    echo "ERR: missing right _2.fq of $lreads"; continue;
  fi 

  echo bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog
  bowtie2 $opts -x $trdb -1 $lreads -2 $rreads \
     --un-conc $odir/$rdnam.unc.fq -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog  &

  i=$(( $i + $ncpuPerMap ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait
# reads loop

exit
#..............................................
##...... sam sort/merge........
#. samset=`ls $odir/$tnam-*.sam`
#. 
#. i=0;
#. for isamu in $samset; do {
#.   inam=`echo $isamu | sed 's/.sam//'`
#.   if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
#. 
#.   ( samtools view -u -S $isamu | samtools sort - $inam.sort; \
#.     samtools index $inam.sort.bam; samtools idxstats $inam.sort.bam > $inam.count; ) &
#. 
#.   i=$(( $i + 1 ))
#.   if [ $i -ge $ncpu ]; then wait; i=0; fi
#.   
#. } done
#. wait
# 
#....... update bowtie2 methods, see rna-express example, using -a all-hits slow but needed?
# bowtie-build --offrate 1 transcripts.fasta transcripts
##  build --offrate 1 for faster multi-map, bigger index
# bowtie -aS -X 800 --offrate 1 transcripts -1 rds_1.fastq -2 rds_2.fastq | samtools view -Sb - > hits.bam
##  -a = --all hits;
#  express transcripts.fasta hits.bam

