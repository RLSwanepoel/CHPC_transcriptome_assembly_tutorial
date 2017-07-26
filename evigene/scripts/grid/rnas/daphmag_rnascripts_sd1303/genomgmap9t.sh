#! /bin/bash
### env fasta=pathto/myset.fa workd=path/to/work qsub -q normal genomgmap9t.sh
#PBS -N genomgmap9
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=9:55:00
#PBS -o genomgmap.$$.out
#PBS -e genomgmap.$$.err
#PBS -V

ncpu=32

# trestles
scrad=$HOME/scratcht
genod=$scrad/chrs/daphmag
bindir=$HOME/bio/gmap1206/bin

gdb=$genod/genome/gmap12
dgenome=dmagna20100422assembly; gtag=dmag10

if [ "X" = "X$workd" ]; then echo "env workd=path/to missing"; exit -1; fi
outdir=$workd/outf
cd $workd/

##for trmap: gopt="--nosplicing -n 9 -S"; osuf=outns
# gopt="-n 9 -S"; osuf=out
gopt="--suboptimal-score=2 --min-intronlength=29 -n 9 -S"; osuf=out9
# gopt="-n 4 -f 2"; osuf=gff

if [ "X" = "X$fasta" ]; then echo "env fasta= bad"; exit -1; fi

echo "START `date`" 
mkdir $outdir

#... START LOOP estin
for estin in $fasta; do {
  if ! test -f "$estin" ; then continue; fi

  dest=`basename $estin .tr | sed 's/\.fasta//; s/\.fa//; s/\.cds//;'`
  outf=$dest-$gtag.gmap.$osuf
  
  i=0; while [ $i -lt $ncpu ]; do { 
  $bindir/gmap $gopt -D $gdb -d $dgenome --part=$i/$ncpu $estin > $outdir/$dest.$gtag.gmap$i.$osuf & 
  i=$(( $i + 1 ))
  }
  done
  wait
  
  cat $outdir/$dest.$gtag.gmap*.$osuf > $outf
  gzip --fast $outf

} done
#... END LOOP estin

echo "DONE `date` "

