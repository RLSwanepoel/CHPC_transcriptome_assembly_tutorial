#! /bin/bash
### env gsnapset=Dman_80_GCCAAT_L006_1  datad=`pwd`/gsodm13cxpa qsub -q normal ...
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=1:ppn=16,walltime=7:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

## change to shared? NOT ON gordon; unavail que; only conmult is big, takes long; revise to run 2+ gsnapset
## 01:32 hr for Dman_80; 30Gb conmult.bam, 2.3Gb conuniq.bam, others <1Gb
## samtools merge for gsnap2genes.sh new split-output
## .. merge subparts to each split-part.bam; skip empties (transloc); 11 valid sparts
## ls -tlh gsodm13cxpa/D*.gsnap7.* 
#.  7.8G Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.concordant_mult
#.  154M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.concordant_uniq
#.  108M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.halfmapping_mult
#.  7.4M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.halfmapping_uniq
#.  6.7M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.nomapping
#.  242M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.paired_mult
#.   17M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.paired_uniq_inv
#.  3.5K Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.paired_uniq_long
#.  1.5M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.paired_uniq_scr
#.  107M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.unpaired_mult
#.  8.9M Mar 19 16:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.unpaired_uniq
#.     0 Mar 19 04:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.concordant_transloc
#.     0 Mar 19 04:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.halfmapping_translo
#.     0 Mar 19 04:31 gsodm13cxpa/Dman_80_GCCAAT_L006_1.gsnap7.unpaired_transloc

bigone="concordant_mult"
keepset="concordant_uniq halfmapping_uniq paired_mult paired_uniq_long"
poorset="halfmapping_mult paired_uniq_inv paired_uniq_scr unpaired_mult unpaired_uniq nomapping"
skipset="concordant_transloc halfmapping_transloc unpaired_transloc"  
smallkeep="$keepset $poorset"
allkeep="$bigone $keepset $poorset"

if [ "X" = "X$gsnapset" ]; then echo "ERR: miss gsnapset=gsnapdataprefix"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapset "; exit -1; fi
bindir=$HOME/bio/bin
cd $datad/

for ssuf in $allkeep; do {
  outs=$gsnapset-gsnap-$ssuf
  ( cat $gsnapset.gsnap*.$ssuf | $bindir/samtools view -1 -o $outs.bam -S - ; \
    $bindir/samtools flagstat $outs.bam > $outs.fstat; ) &
} done

wait

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
## ------- old -------------
## drna=$gsnapset
## dont need this now ; --samheader=part0
# grep '^\@' $drna.gsnap0.concordant_uniq > $drna.hdr.sam
# ( grep -hv '^\@' $drna.gsnap*.$ssuf | cat $drna.hdr.sam - |\
# $bindir/samtools view -1 -o $outs.bam -S - ; $bindir/samtools flagstat $outs.bam > $outs.fstat; )&
##.......
# grep '^\@' $drna.gsnap0.concordant_uniq > $drna.hdr.sam
# cat $drna.gsnap*.{concordant_mult,concordant_uniq,halfmapping_uniq,paired_mult,paired_uniq_long} |\
#  grep -v '^\@' | cat $drna.hdr.sam - | $bindir/samtools view -1 -o $out1.bam -S -
# $bindir/samtools flagstat $out1.bam > $out1.fstat  


