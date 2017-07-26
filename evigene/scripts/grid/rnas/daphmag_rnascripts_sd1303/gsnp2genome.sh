#! /bin/bash
### env fastdir=fastq/ qsub -q normal gsnapset.sh
#PBS -A ind114
#PBS -N gsnp4
#... -l nodes=1:ppn=32,walltime=29:55:00
#PBS -l nodes=2:ppn=16,walltime=29:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

ncpu=32
npart=32
jstart=0
# jend=$(( $jstart + $ncpu ))

#trestle#  datad=$HOME/scratchn
datad=$HOME/scratchg
workd=$datad/chrs/daphmag
bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1206
rund=$workd/rnas
gdb=gmap12
gdbpath=$workd/genome/

# run5 on genome
snpd=snp10
dgenome=dmagna20100422assembly; gtag=dmag2
dgenosize=$dgenome.chr_size.txt

#5 for genome map, .sam out, with introns, no snps print
optsnp=""
# "--use-snps=$snpd "
## dmag fq4: helski rna 
# optreads="--pairexpect=475 --pairdev=50"
## dmag fq5:  nodame rna
optreads="--pairexpect=375 --pairdev=50"
optgenome="-N 1  " 
snapopt1="$optgenome $optsnp $optreads" 

cd $workd/rnas/
if [ "X$fastdir" = "X" ]; then
  fastdir="fastq"
fi

## tie outdir to fastdir or not?
outs=`basename $fastdir`
outdir=$workd/rnas/gsnout$outs
mkdir $outdir

#... START LOOP rnain : limit to suffix .fq, .fastq, sequence.txt ..

#xx# for rnain in $fastdir/*{_1.fastq,.1.fastq}.gz ;  do 
for rnain in $fastdir/*_1.fastq.gz ;  do 
{ 
  cd $rund/
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/.fastq//; s/.fq//; '`
  outna=$drna-$gtag

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fastq/2.fastq/; '`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  ## elif $paironly .. fail
  else 
    inset=$rnain
  fi

  # i=1,ncpu, j=jstart,jend  subset of npart
  i=0;  
  j=$jstart;
  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A "sam"  --part=$j/$npart \
    -D $gdbpath/$gdb -d $dgenome $inset > $outdir/$drna.gsnap$i.samu &
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  wait

  ##...... sam sort/merge........

  samset=$drna
  cd $outdir
  i=0; while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $gdbpath/$dgenosize $isamu | $bindir/samtools sort - $inam.sort; )&
    # .sort may fail#  rm $isamu 
  
    i=$(( $i + 1 ))
  }
  done
  wait

  $bindir/samtools merge $outna.bam $samset*.sort.bam
  $bindir/samtools flagstat $outna.bam > $outna.fstat
  $bindir/samtools index $outna.bam 
} 
done

#... END LOOP rnain

