#! /bin/bash
### env outdir=gsnoutfq4inb drna=Dman_13_ACAGTG_L006_1  qsub -q normal gsnp2merge.sh
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=1:ppn=32,walltime=3:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

## finish up last partial gsnap samu set: sam.sort,merge

ncpu=32
jstart=0

datad=$HOME/scratchn
workd=$datad/chrs/daphmag
bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1206
rund=$workd/rnas
gdb=gmap12
gdbpath=$workd/genome/

# run5 on genome
dgenome=dmagna20100422assembly; gtag=dmag2
dgenosize=$dgenome.chr_size.txt

#. optsnp="" optreads="" optgenome="-N 1  " 
#. snapopt1="$optgenome $optsnp $optreads" 

cd $workd/rnas/
# below# cd $outdir
#. if [ "X$fastdir" = "X" ]; then fastdir="fastq" fi

## tie outdir to fastdir or not?
#. outs=`basename $fastdir`
#. outdir=$workd/rnas/gsnout$outs
#. mkdir $outdir

# for rnain in $fastdir/*_1.fastq.gz ;  do 
# { 
  # cd $rund/
  # if ! test -f "$rnain" ; then continue; fi
  
  # drna=`basename $rnain .gz | sed 's/.fastq//; s/.fq//; '`
  outna=$drna-$gtag

  ##...... sam sort/merge........

  samset=$drna
  cd $outdir
  i=0; 
  while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $gdbpath/$dgenosize $isamu | $bindir/samtools sort - $inam.sort; ) &
    # later rm isamu
 
    i=$(( $i + 1 ))
  }
  done
  wait

  $bindir/samtools merge $outna.bam $samset*.sort.bam
  $bindir/samtools flagstat $outna.bam > $outna.fstat
  $bindir/samtools index $outna.bam 

# } done
#... END LOOP rnain

