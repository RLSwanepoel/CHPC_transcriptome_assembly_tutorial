#! /bin/bash
### env gsnapdir="gsodm13xpa gsodm13xtr" datad=`pwd` qsub -q normal gsnp2bam3n.sh
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=1:ppn=16,walltime=7:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

ncpu=16

# gsnp2bam3n.sh : rev to handle gsnapdir=list with 1,2+ data groups/dir
#  gsnapdir="gsodm13xpa gsodm13xtr" ..

## change to shared? NOT ON gordon; unavail que; only conmult is big, takes long; revise to run 2+ gsnapset
## 01:32 hr for Dman_80; 30Gb conmult.bam, 2.3Gb conuniq.bam, others <1Gb
## samtools merge for gsnap2genes.sh new split-output
## .. merge subparts to each split-part.bam; skip empties (transloc); 11 valid sparts

bigone="concordant_mult"
keepset="concordant_uniq halfmapping_uniq paired_mult paired_uniq_long"
poorset="halfmapping_mult paired_uniq_inv paired_uniq_scr unpaired_mult unpaired_uniq nomapping"
# skipset="concordant_transloc halfmapping_transloc unpaired_transloc"  
smallkeep="$keepset $poorset"
# allkeep="$bigone $keepset $poorset"

# if [ "X" = "X$gsnapset" ]; then echo "ERR: miss gsnapset=gsnapdataprefix"; exit -1; fi
if [ "X" = "X$gsnapdir" ]; then echo "ERR: miss gsnapdir=gsnap dir list "; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapdirs"; exit -1; fi
bindir=$HOME/bio/bin

cd $datad/

# gsnp2bam3n.sh : rev to handle gsnapdir=list with 1,2+ data groups/dir
#  gsnapdir="gsodm13xpa gsodm13xtr" ..

icpu=0;
for gdir in $gsnapdir ; do { 
  echo cd $datad/$gdir; 
  cd $datad/$gdir; 
  glist=`ls *.gsnap1.concordant_mult`; 
  for gsfile in $glist; do { 
    gsnapset=`echo $gsfile | sed 's/\.gsnap1.concordant_mult//;'`

    ssuf=$bigone
    outs=$gsnapset-gsnap-$ssuf
    ( echo cat $gsnapset.gsnap0.$ssuf TO  $bindir/samtools view -1 -o $outs.bam -S - ; \
      echo $bindir/samtools flagstat $outs.bam TO  $outs.fstat; ) &
    echo forked;
   
    # Ooops .. this whole loop is 1 fork cmd ..
    ( for ssuf in $smallkeep; do { outs=$gsnapset-gsnap-$ssuf; \
       echo cat $gsnapset.gsnap0.$ssuf TO  $bindir/samtools view -1 -o $outs.bam -S - ; \
       echo $bindir/samtools flagstat $outs.bam TO  $outs.fstat;  \
    } done ) &
    echo forked

    ## innerloop per readfile
    icpu=$(( $icpu + 2 )); if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi
    echo atcpu  $icpu

   } done
   # gsfile

  cd $datad
  echo
} done
# gsnapdir

echo waiting
wait
echo done wait

# # ... oneway: forkeach part ....
# for ssuf in $allkeep; do {
#   outs=$gsnapset-gsnap-$ssuf
#   ( cat $gsnapset.gsnap*.$ssuf | $bindir/samtools view -1 -o $outs.bam -S - ; \
#     $bindir/samtools flagstat $outs.bam > $outs.fstat; ) &
# } done
# wait
#
# #--- otherway for 2 cpu shared que -------
# ssuf=$bigone
# outs=$gsnapset-gsnap-$ssuf
# ( cat $gsnapset.gsnap*.$ssuf | $bindir/samtools view -1 -o $outs.bam -S - ; \
#   $bindir/samtools flagstat $outs.bam > $outs.fstat; ) &
# 
# ## no fork here; ncpu=2
# for ssuf in $smallkeep; do {
#   outs=$gsnapset-gsnap-$ssuf
#   ( cat $gsnapset.gsnap*.$ssuf | $bindir/samtools view -1 -o $outs.bam -S - ; \
#     $bindir/samtools flagstat $outs.bam > $outs.fstat; ) 
# } done
# 
# wait
# 
