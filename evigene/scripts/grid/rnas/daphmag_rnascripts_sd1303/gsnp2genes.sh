#! /bin/bash
### env trdb=pathtogmapdb readset=fastq/xxx*_1.fastq.gz  qsub -q normal gsnp2genes.sh
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=1:ppn=32,walltime=29:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

ncpu=32
npart=32
jstart=0
maxmem=1000000000

if [ "X" = "X$trdb" ]; then echo "ERR: miss trdb=pathtogmapdb"; exit -1; fi
if [ "X" = "X$readset" ]; then echo "ERR: miss readset=fqs/xxx*_1.fastq.gz"; exit -1; fi

sopt="-n"; ssuf="nsort"; # name sort
# sopt=""; ssuf="sort"; # ref/loc sort

#trestle#  
datad=$HOME/scratcht
workd=$datad/chrs/daphmag
## make param..
rund=$workd/rnas

bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1206
#new.bad#gmapd=$HOME/bio/gmap1303
## new opt 1303:  --maxsearch=1000 >= npaths ; smaller == speedup 
## GONE opt:  --pairexpect=450 --pairdev=50 
## new: SEGfaults .. no clue why..  NO GOOD cant get rid of segfaults

# *** need env params here ..
## gdbpath=$workd/genome/gmap12
## dgenome=dmagna20100422assembly; gtag=dmag2
gdbpath=`dirname $trdb`
dgenome=`basename $trdb`
## daphmag5xall13f_okall
gtag=`echo $dgenome | sed 's/_.*//;'`

#5 for genome map, .sam out, with introns, no snps print
optsnp=""
## dmag fq4: helski rna ; was 475..
optreads="--pairexpect=450 --pairdev=50"
#new.bad#optreads=""
## dmag fq5:  nodame rna
# optreads="--pairexpect=375 --pairdev=50"
## genome -N; transcript -N 0
## optgenome="-N 1  " 
#o# optgenome="--gmap-mode=none --npaths=200 " 
#n# optgenome="--gmap-mode=none --maxsearch=400 --npaths=200 " 
optgenome="--gmap-mode=none " 
snapopt1="$optgenome $optsnp $optreads" 

cd $rund

# tnam=`basename $trdb .tr`
reada=($readset)
rnam=`dirname ${reada[0]}`
outdir=gso$rnam
# outna=$odir/$tnam-$rnam

mkdir $outdir

#... START LOOP rnain : limit to suffix .fq, .fastq, sequence.txt ..

for rnain in $readset ;  do 
{ 
  cd $rund/
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/.fastq//; s/.fq//; '`
  outna=$drna-$gtag

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fastq/2.fastq/; '`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    ## elif $paironly .. fail
    inset=$rnain
  fi

  # i=1,ncpu, j=jstart,jend  subset of npart
  i=0;  
  j=$jstart;

  echo "CMD: $gmapd/bin/gsnap $snapopts -A sam -D $gdbpath -d $dgenome $inset TO $outdir/$drna.gsnap0.samu "

  ##? change output to --split-output=$outdir/$drna.gsnap$i  ??

  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A sam  --part=$j/$npart \
    -D $gdbpath -d $dgenome --split-output=$outdir/$drna.gsnap$i $inset &

    # OLD# -D $gdbpath -d $dgenome $inset > $outdir/$drna.gsnap$i.samu 
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  wait

#.....................
  ##...... sam sort/merge........
  ## FIXME: want sort -n option for namesort; ALSO merge -n

samok=0;
echo "SKIP sam sort/merge"
if [ $samok = 1 ]; then
  samset=$drna
  cd $outdir
  i=0; while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`

    if  test -s $isamu ; then  
      ( $bindir/samtools view -u -S $isamu | $bindir/samtools sort $sopt  -m $maxmem - $inam.$ssuf; )&
    else
      echo "missing $isamu"; 
    fi 

    i=$(( $i + 1 ))
  }
  done
  wait

  if [ -s $samset.gsnap1.$ssuf.bam ]; then
        ## ?? outamem here, 78 GB vmem, 65 GB mem; -1 = gzip --fast
  $bindir/samtools merge -1 $sopt $outna.bam $samset*.$ssuf.bam
  $bindir/samtools flagstat $outna.bam > $outna.fstat
  #later# $bindir/samtools index $outna.bam 
  fi
fi

} 
done

#... END LOOP rnain

