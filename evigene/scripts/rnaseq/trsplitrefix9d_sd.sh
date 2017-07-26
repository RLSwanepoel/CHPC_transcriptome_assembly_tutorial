#! /bin/bash
### trsplitrefix9d_sd.sh
### env strand=rev splitdir=spltrasmn4f  qsub -q batch trsplitrefix9d_sd.sh
#PBS -N trsplitfix
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=11:55:00
#PBS -o trsplitfix.$$.out
#PBS -e trsplitfix.$$.err
#PBS -V

#.. trsplitrefix9d_sd.sh: redo trasmfixup.pl, after trsplit/trgmap, adding trchimerpick before new genefindcds, overbestgene1
#.. input is tr.fa and trfa.gmap.gff, from last trsplitfix, in spl$trnam/ output dirs,
#.. eg: input -gff spltrasmn4f/trasmn4f.tr.split.19.fa.split.gff -trfa spltrasmn4f/trasmn4f.tr.split.19.fa.split
# ** fix name clashes with last run .. intermediate outputs .an.gff, final best$dv.gff ok?

ncpu=32
export TMPDIR=$HOME/scratchn/tmp
export LC_ALL=C 

# dv=j9
dv=m9
config=trsplitbest9sd.config
subd=strandgroups/asmfull

datad=$HOME/scratchn
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
evigene=$HOME/bio/evigene/
sdir=$evigene/scripts/rnaseq

cd $rund
echo "start trsplitrefix9d: `date`"  

odir=$splitdir
trnam=`echo $splitdir | sed 's/^spl//;'`
##trnam=`basename $trsin | sed 's/\..*//;'`
if [ ! -d $odir ]; then
  echo "ERROR: missing directory: $odir"; exit -1;
fi
cd $odir
# fixme:
config="../$config";

opstrand=
if [ "X$strand" = "X" ] ; then
  if [ -f $trnam*.1.fa.split.fwd.an.gff ]; then strand="fwd"; fi
  if [ -f $trnam*.1.fa.split.rev.an.gff ]; then strand="rev"; fi
fi
if [ ! "X$strand" = "X" ] ; then
  opstrand="-strand $strand"
fi
# trstr=`echo $trsin | sed 's/.*fall.tr/fwd/; s/.*rall.tr/rev/; s/.*f.tr/fwd/; s/.*r.tr/rev/;'`
# if [ $trstr = "fwd" ] || [ $trstr = "rev" ]; then  strand="-strand $trstr"; fi

#. optsplit="$strand -nobest -splitjoins"
optsplit="$opstrand -nobest -chimerapick"

if [ ! -f $trnam*.1.fa.split.gff ]; then
  echo "ERROR: missing input data: $trnam .fa.split.gff"; exit -1;
fi

trlist=`/bin/ls $trnam*.fa.split`

i=0; # cpu counter
for trfa in $trlist; do
{
  i=$(( $i + 1 ))
  nam=`echo $trfa | sed 's/.fa.split//;'`
  rlog=trfix$dv.$nam.log

  if [ ! -f $trfa.gff ]; then
  echo "ERROR: missing input $trfa.gff"; continue;
  fi
  if [ -f $trfa.$strand.an.gff ]; then
  mv $trfa.$strand.an.gff $trfa.$strand.an.gff1  
  #dontcare# mv $trfa.$strand.gff $trfa.$strand.gff1
  fi
  
  echo $sdir/trasmfixup.pl $optsplit -vers $dv -config $config -trfa $trfa -gff $trfa.gff
  $sdir/trasmfixup.pl $optsplit -vers $dv -config $config -trfa $trfa  -gff $trfa.gff >& $rlog &

  #above#i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

## FIXME: need run -bestonly on all splits together, split.*.an.gff inputs
cat $trnam*.split.*.an.gff | $sdir/trasmfixup.pl -bestonly -name $trnam -vers $dv -config $config -gff stdin

# # package up results: cat the *.split.*.xxx parts to one file per xxx ?
##mv trfix.$trnam*.log $odir/
##mv $trnam*.split.* $odir/

echo "end trsplitrefix9d: `date` "
