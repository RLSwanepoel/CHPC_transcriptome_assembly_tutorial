#! /bin/bash
### qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

ncpu=8

# dv=4nal; fasta5/all_nodame.dnorm30
# dv=4nal; inpe=fasta5/all_nodam*.dnorm30.fa2.gz
dv=5nun1; inpe=nodamx_allunc 
# inpe2=soap5nun1/nodamx_allunc.2.fq

kset="31 25 29 21 75 65 55 45"
# F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

## trestles
datad=$HOME/scratcht
workd=$datad/chrs/daphmag/rnas
bindir=$HOME/bio/soaptrans
rund=soap$dv
oname="dmag$dv"

cd $workd/
# mkdir $rund
## check for existing rund/pereads.fa
if [ ! -f $rund/$inpe.1.fq ]; then
 echo "missing $rund/$inpe.1.fq"; exit -1
fi

cat > $rund/rdconfig  <<EOT
max_rd_len=110
[LIB]
rank=1
asm_flag=3
avg_ins=350
q1=$inpe.1.fq
q2=$inpe.2.fq
EOT

cd $rund

for kmer in $kset; do {
 odir=sod$kmer
 outname=so${oname}.k$kmer
 mkdir $odir
if [ $kmer -lt 32 ]; then
  $bindir/SOAPdenovo-Trans-31kmer all -s rdconfig -o $outname -K $kmer $sopt
else 
  $bindir/SOAPdenovo-Trans-127mer all -s rdconfig -o $outname -K $kmer $sopt
fi

 mv $outname* $odir/
 if test -f $odir/$outname.contig ; then
  rm $odir/$outname.{preArc,vertex,edge.gz,Arc,updated.edge}
  gzip --fast $odir/$outname.{readOnContig,readInGap,ctg2Read,scafSeq,contig}
 fi

} done

