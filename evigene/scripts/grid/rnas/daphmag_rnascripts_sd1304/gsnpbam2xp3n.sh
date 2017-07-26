#! /bin/bash
### env gsnapdir="gsodm13xpa gsodm13xtr" datad=`pwd` qsub -q normal gsnpbam2xp3n.sh
#PBS -A ind114
#PBS -N rnax4
#PBS -l nodes=1:ppn=16,walltime=19:55:00
#PBS -o rnax4.$$.out
#PBS -e rnax4.$$.err
#PBS -V

ncpu=16
ncpu1=4
# gsnpbam2xp3n.sh : rnaexpress runs #  gsnapdir="gsodm13xpa gsodm13xtr" ..
## NO good: rnaexpress ... in1.bam in2.bam : only handles 1 dammit.. << NO, uses "in1.bam,in2.bam,.."
## .. still no good: fails w/ other error for 2+ bams ;; 
# problem may be paired_ = ERROR: Input BAM file contains no valid alignments.

xprset="concordant_mult concordant_uniq paired_mult paired_uniq_long paired_uniq_inv"
xpradd="concordant_uniq"
#NOxpradd="concordant_uniq paired_mult paired_uniq_long paired_uniq_inv"
# bigone="concordant_mult"
# keepset="concordant_uniq halfmapping_uniq paired_mult paired_uniq_long"
# poorset="halfmapping_mult paired_uniq_inv paired_uniq_scr unpaired_mult unpaired_uniq nomapping"
# skipset="concordant_transloc halfmapping_transloc unpaired_transloc"  
# allkeep="$bigone $keepset $poorset"

if [ "X" = "X$gsnapdir" ]; then echo "ERR: miss gsnapdir=gsnap dir list "; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapdirs"; exit -1; fi
if [ "X" = "X$trfasta" ]; then echo "ERR: missing env trfasta=mydata.tr"; exit -1; fi
biobin=$HOME/bio/bin

cd $datad/

tnam=`basename $trfasta .tr | sed 's/\.fasta//; s/\.fa$//; s/_ok.*//;'`
# datad=xxx/gsodm13sd  gdir=gsodm13cxpa_dn1/Dman_80_GCCAAT_L006_1-gsnap-concordant_mult.bam

echo "START rnax: `date`"
icpu=0;
for gdir in $gsnapdir ; do { 
  echo cd $datad;
  cd $datad

  glist=`ls $gdir/*concordant_mult.bam`; 
  for bam1 in $glist; do { 

    bnam=`echo $bam1| sed 's/.concordant_mult.bam//;'`
    bams=$bam1; for bsuf in $xpradd; do { bams="$bams,$bnam-$bsuf.bam"; } done
    bnamo=`basename $bnam | sed 's/.gsnap//; s/_1$//;'`
    outrx=rx$tnam-$bnamo

    ## no.need.app.does# echo mkdir $outrx

    ## which way? samerge takes a while at 1 cpu... other way fails? or can do only bam1
    # bamerge=$bnamo.merge.bam
    # bams=`echo $bams | sed 's/,/ /g;'`;
    # ( echo samtools merge -1 -n $bamerge $bams; \
    #   echo $biobin/rnaexpress -p $ncpu1 -o $outrx --no-update-check $trfasta $bamerge ) &

    echo $biobin/rnaexpress -p $ncpu1 -o $outrx --no-update-check $trfasta $bams 
    $biobin/rnaexpress -p $ncpu1 -o $outrx --no-update-check $trfasta $bams  &
    # echo forked;
   
    ## innerloop per readfile
    icpu=$(( $icpu + $ncpu1 )); if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi
    # echo atcpu  $icpu

   } done
   # glist

  cd $datad
} done
# gsnapdir

# echo waiting
wait
# echo done wait
echo "DONE rnax: `date`"

