#! /bin/bash
### qsub -q batch trasmfix9_sd.sh
#PBS -N trasmfix
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=23:55:00
#PBS -o trasmfix.$$.out
#PBS -e trasmfix.$$.err
#PBS -V

ncpu=8
export TMPDIR=$HOME/scratchn/tmp
export LC_ALL=C 

dv=a9
config=trasmbest9sd.config
strand=fwd
# strand=rev
# optrasm="-ver $dv -grepread mismatch:1 -allstrandspan"

subd=strandgroups/asmfull

datad=$HOME/scratchn
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
sdir=$HOME/bio/evigene/scripts/rnaseq

cd $rund

trlist=`ls *.tr.gz`
gflist=`ls *.gff.gz`
# trasmb9fall.fwd.gff.gz trasmb9fall.tr.gz     
# dmag2vel4f1xinco.gff.gz dmag2vel4f1xinco.tr.gz 

echo "start trasmfix: `date`"  
i=0; # cpu counter

for gff in $gflist; do
{
  nam=`echo $gff | sed 's/.gz//; s/.gff//; s/.fwd//; s/.rev//;'`
  trfa=$nam.tr.gz 
  rlog=trfix.$nam.log

  echo $sdir/trasmfixup.pl -vers $dv -config $config -gff $gff -trfa $trfa
  $sdir/trasmfixup.pl -vers $dv -config $config -gff $gff -trfa $trfa >& $rlog &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

echo "end trasmfix: `date` "

