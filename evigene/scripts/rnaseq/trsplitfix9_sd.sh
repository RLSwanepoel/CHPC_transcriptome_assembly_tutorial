#! /bin/bash
### trsplitfix9_sd.sh
### env trsin=xxx.tr  qsub -q batch trsplitfix9_sd.sh
#PBS -N trsplitfix
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=11:55:00
#PBS -o trsplitfix.$$.out
#PBS -e trsplitfix.$$.err
#PBS -V

## run per trasm.tr .. split to ncpu parts; per part: splitjoins.tr, gmap.tr, findcds.tr, best.tr
## FIXME: need run -bestonly on all splits together, split.*.an.gff inputs
##  cat xxx*.an.gff | $sdir/trasmfixup.pl -bestonly -vers $dv -config $config -trfa $name -gff stdin

ncpu=32
export TMPDIR=$HOME/scratchn/tmp
export LC_ALL=C 

dv=j9
config=trsplitbest9sd.config
subd=strandgroups/asmfull

datad=$HOME/scratchn
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
evigene=$HOME/bio/evigene/
sdir=$evigene/scripts/rnaseq

cd $rund
echo "start trsplitfix: `date`"  

## be tidy? put parts into new subdir?
trnam=`basename $trsin | sed 's/\..*//;'`
odir="spl$trnam"
mkdir $odir

# gunzip $trsin.gz
strand=
trstr=`echo $trsin | sed 's/.*fall.tr/fwd/; s/.*rall.tr/rev/; s/.*f.tr/fwd/; s/.*r.tr/rev/;'`
if [ $trstr = "fwd" ] || [ $trstr = "rev" ]; then  strand="-strand $trstr"; fi

optsplit="$strand -nobest -splitjoins"

if [ ! -f $trsin.split.1.fa ]; then
 pindir=`dirname $trsin`
 splitsize=`grep -v '^>' $trsin | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $trsin
fi

trlist=`/bin/ls $trsin.split.*.fa`

i=0; # cpu counter
for trfa in $trlist; do
{
  nam=`echo $trfa | sed 's/.fa//;'`
  rlog=trfix.$nam.log

  echo $sdir/trasmfixup.pl $optsplit -vers $dv -config $config -trfa $trfa
  $sdir/trasmfixup.pl $optsplit -vers $dv -config $config -trfa $trfa >& $rlog &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

## FIXME: need run -bestonly on all splits together, split.*.an.gff inputs
cat $trnam*.split.*.an.gff | $sdir/trasmfixup.pl -bestonly -trfa $trnam -vers $dv -config $config -gff stdin

# # package up results: cat the *.split.*.xxx parts to one file per xxx ?
mv trfix.$trnam*.log $odir/
mv $trnam*.split.* $odir/

echo "end trsplitfix: `date` "
