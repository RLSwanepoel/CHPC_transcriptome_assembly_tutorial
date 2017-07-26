#! /bin/bash
### env trdb=pathtogmapdb readset=fastq/xxx*_1.fastq.gz rnam=yyy rund=path/to  qsub -q normal gsnp2genes.sh
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=3:ppn=16,walltime=39:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

## new params: readset fails/last only ; need readseq="fq/xx*1.fq fq/yy*1.fq" not {curlies}
#  env readset=fq4x/{Dman_42,Dman_64,Dman_65}*_1.fastq.gz rnam=xfi trdb=trdb/xxx rund=`pwd` qsub ..
## takes 12h .. 15+hr for one readset, 48 cpu

ncpu=48
npart=48
jstart=0

if [ "X" = "X$trdb" ]; then echo "ERR: miss trdb=pathtogmapdb"; exit -1; fi
if [ "X" = "X$readset" ]; then echo "ERR: miss readset=fqs/xxx*_1.fastq.gz"; exit -1; fi
#gordon#
if [ "X" = "X$rund" ]; then rund=$HOME/scratchg/chrs/daphmag/rnas; fi
cd $rund

# keepset="concordant_mult,concordant_uniq,halfmapping_uniq,paired_mult,paired_uniq_long"
# poormap="halfmapping_mult,paired_uniq_inv,paired_uniq_scr,unpaired_mult,unpaired_uniq,nomapping"
# skipset="concordant_transloc,halfmapping_transloc,unpaired_transloc"  # these are empty probably
#
# #---- split-output parts ----- samsort some of these, not others ----
#  8.9G  gsnap9.concordant_mult 	Yes
#  9.0M  gsnap9.concordant_transloc	No/Yes?
#  208M  gsnap9.concordant_uniq		Yes
#  208M  gsnap9.halfmapping_mult        No/Yes?
#  9.0M  gsnap9.halfmapping_transloc
#   21M  gsnap9.halfmapping_uniq	Yes maybe ?
#   31M  gsnap9.nomapping		NO
#  175M  gsnap9.paired_mult		Yes/No ?
#   28M  gsnap9.paired_uniq_inv
#  9.0M  gsnap9.paired_uniq_long	Yes maybe?
#   11M  gsnap9.paired_uniq_scr
#   98M  gsnap9.unpaired_mult		No/Yes ?
#  9.0M  gsnap9.unpaired_transloc
#   21M  gsnap9.unpaired_uniq
## unpaired are pairs mapping b/n 2 tr:
# 14460:98619	65	dmag5xun1velvk65Loc2112t1	82	40	101M	dmag5xun3cavelvk81Loc464t8
# 14460:98619	129	dmag5xun3cavelvk81Loc464t8	42	40	63M38S	dmag5xun1velvk65Loc2112t1
# paired_uniq_inv,_scr probably for keepset : paired to same tr; _long means 1000bp inner span
#-----------

bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1206

#new.bad#gmapd=$HOME/bio/gmap1303
## new opt 1303:  --maxsearch=1000 >= npaths ; smaller == speedup 
## GONE opt:  --pairexpect=450 --pairdev=50 
## new: SEGfaults .. no clue why..  NO GOOD cant get rid of segfaults

gdbpath=`dirname $trdb`
dgenome=`basename $trdb`
gtag=`echo $dgenome | sed 's/_.*//;'`

## add sam header controls:
#  --no-sam-headers               Do not print headers beginning with '@'
#  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q
## ^^ this for split-output still bad need to grep out '^@' for merge.
#
#5 for genome map, .sam out, with introns, no snps print
optsnp=""

## dmag fq4: helski rna ; was 475..
optreads="--pairexpect=450 --pairdev=50"

#new.bad#optreads=""
## dmag fq5:  nodame rna
# optreads="--pairexpect=375 --pairdev=50"
## genome -N 1; transcript -N 0
## optgenome="-N 1  " 
#o# optgenome="--gmap-mode=none --npaths=200 " 

optgenome="--gmap-mode=none " 
snapopt1="--sam-headers-batch=0 $optgenome $optsnp $optreads" 

cd $rund

readfiles=`ls $readset`
if [ "X" = "X$rnam" ]; then
  reada=($readset)
  rnam=`dirname ${reada[0]}`
fi

outdir=gso$rnam
mkdir $outdir

# for rnain in $readset ;  do 
for rnain in $readfiles ;  do 
{ 
  cd $rund/
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/\.fastq//; s/\.fq//; '`
  outna=$drna-$gtag

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fastq/2.fastq/; '`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    echo "ERR: missing pair file2 $rnain2"; continue;
  fi

  i=0;  j=$jstart;

  echo "CMD: $gmapd/bin/gsnap $snapopts -A sam -D $gdbpath -d $dgenome $inset TO $outdir/$drna.gsnap0.samu "
  ##* change output to --split-output=$outdir/$drna.gsnap$i

  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A sam  --part=$j/$npart \
    -D $gdbpath -d $dgenome --split-output=$outdir/$drna.gsnap$i $inset &

   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  wait

} 
done

#... END LOOP rnain

