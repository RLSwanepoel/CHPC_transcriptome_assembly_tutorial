#! /bin/bash
### env inpe1=pathto/reads.1.fq qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

ncpu=6

##.. ~13Gb 1,2.fq ; soap.k31 ~50G
# dv=5xun1; inpe=hiskx2l47_unc
# dv=5xun2; inpe=hiskx2l36_unc
#.. obs ins sizes hiskx: 200..250; add 200 read pair = 425
# n=1519766; ains=288; mdins=265; n=1482344; ains=283; mdins=261; n=1456546; ains=224; mdins=203

## pergroup unc.[12].fq data sets: fq4xun/
# xinb3_bn.unc.1.fq xinb3_bx.unc.1.fq xinb3_ca.unc.1.fq xinb3_copt2.unc.1
# xinb3_co.unc.1.fq xinb3_cr.unc.1.fq xinb3_fi.unc.1.fq xinb3_pa.unc.1.fq

if [ "X" = "X$inpe1" ]; then echo "ERR: miss env inpe1=pathto/reads.1.fq"; exit -1; fi
#below# if [ ! -f $inpe1 ]; then echo "ERR: miss env inpe1=$inpe1 "; exit -1; fi

inpe2=`echo $inpe1 | sed 's/1.fq/2.fq/'`
gp=`basename $inpe1 1.fq | sed 's/\..*//;'`
dv=5$gp

kset="31 25 29 21 75 65 55"
## drop hi kmers? failing often..
# F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

## trestles
datad=$HOME/scratcht
workd=$datad/chrs/daphmag/rnas
bindir=$HOME/bio/soaptrans
rund=soap$dv
oname="dmag$dv"

cd $workd/
mkdir $rund
cd $rund

if [ ! -f $inpe1 ]; then echo "ERR: miss env inpe1=$inpe1 "; exit -1; fi
if [ ! -f $inpe2 ]; then echo "ERR: miss inpe2=$inpe2 "; exit -1; fi

cat > rdconfig  <<EOT
max_rd_len=110
[LIB]
rank=1
asm_flag=3
avg_ins=425
q1=$inpe1
q2=$inpe2
EOT

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

