#!/bin/bash
### env trset=transcripts.fa.gz qsub -q shared maketraa.sh
#PBS -N maketraa
#PBS -l vmem=16gb,nodes=1:ppn=2,walltime=5:55:00
#PBS -o maketraa.$$.out
#PBS -e maketraa.$$.err
#PBS -V

## note this is not parallelized, will run quickly on 1 cpu, minimal mem.; 
## cd-hit does have parallel cpu ppts
# mason.indiana.edu config
module add cd-hit
export evigene=/N/u/gilbertd/Mason/projectd/bio/evigene
# use env trset=xxx instead of below trset=$*

# maketraa.sh : find best proteins (ORFs) in transcript set, and cluster to non-redundant longest aa set
# usage: maketraa.sh transcripts.fa.gz
#
# FIXME: use aa flag utrbad/utrpoor here, 
#  1. segregate: ok.aa, poor.aa , 
#  2. cd-hit -i ok.aa -o ok_cd.aa
#  3. cd-hit-2d -i ok_cd.aa -i2 poor.aa -o poor_cd.aa
## fixme2: filter poor.aa also by smallsize
## fixme3: cd-hit -c 0.95 not 0.90 default
## fixme4: need cd-hit on poor.aa before cd-hit-2d ok - poor
## Fixme5: no utrbad/poor annot, parse aalen=236,23%,partial for CDS%

## should these be in relation to trsize? ie very bad = pcds<10..20; bad=pcds<20..40, ..
#def: $MINAA= 30; $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
AMINBAD=200
AMINPOO=100
AAMIN=40

function usage() {
  echo maketraa.sh : find best proteins, ORFs, in transcript set, and cluster to non-redundant longest aa set
  echo usage: maketraa.sh transcripts.fa.gz
  echo part of http://eugenes.org/EvidentialGene/
  exit -1
}

## FIXME: fixed paths no good, check ENV
if [ "X$evigene" = "X" ]; then
evigene=/bio/bio-grid/mb/evigene
fi

cdhitapp=`which cd-hit`
cdhit2dapp=`which cd-hit-2d`
if [ "X$cdhitapp" = "X" ]; then
cdhitapp=/bio/bio-grid/mb/bin/cd-hit
cdhit2dapp=/bio/bio-grid/mb/bin/cd-hit-2d
fi

if [ ! -x $evigene/scripts/cdna_bestorf.pl -o ! -x $cdhitapp ]; then
  echo "ERR: need path to evigene/scripts/cdna_bestorf.pl or cd-hit "; usage;
fi

## add option to input .aa instead of .tr and skip cdna_bestorf..
#above# trset=$*
# input may be .tr or .tr.gz

for faz in $trset; do { 
  if [ ! -f $faz ]; then
   echo "ERR: missing input transcripts: $faz"; usage;
  fi
  # run in tr directory?
  tdir=`dirname $faz`
  #NOT basename, keep path# nam=`basename $faz .gz | sed 's/\.fasta//; s/\.fa//; s/\.tr//;'`
  nam=`echo $faz | sed 's/\.gz//; s/\.fasta//; s/\.fa//; s/\.tr//;'`
  if test -f $nam.aa ; then  /bin/mv $nam.aa $nam.aa.old; fi
  cd $tdir

  $evigene/scripts/cdna_bestorf.pl -nostop -minaa=$AAMIN -cdna $faz > $nam.aa

  #specialcase# perl -pi -e'if(/^>/) { s/\si=\d+;//; s/\s+flag=.+$//; s/\s+path=.+$//; }' $nam.aa

  cat $nam.aa | env mpoo=$AMINPOO mbad=$AMINBAD nam=$nam perl -ne\
'BEGIN{ $mpoo=$ENV{mpoo}; $mbad=$ENV{mbad}; $nam=$ENV{nam}; 
open(OK,">$nam.ok.aa"); open(BAD,">$nam.poor.aa"); }
if(/^>/) { $bad=(/utr(poor|bad)/)?1:0; ($al)=m/aalen=(\d+)/; $skip=0; 
if($al>0) { ($ap)= m/aalen=$al,(\d+)/; 
if($bad) { $skip=1 if((/utrbad/ and $al<$mbad) or (/utrpoor/ and $al<$mpoo)); }
elsif($ap) { $bad=1 if($ap<60); $skip=1 if(($ap<=33 and $al<$mbad) or ($ap<60 and $al<$mpoo)); } } }
if($skip){ } elsif($bad) { print BAD $_; } else { print OK $_; }'

  $cdhitapp -c 0.95 -d 0 -i $nam.ok.aa -o ${nam}.ok_cd.aa >& $nam.cd1.log

  ## should 2d use ok_cd.aa or ok.aa ?
  $cdhitapp -c 0.9 -d 0 -i $nam.poor.aa -o $nam.poorcd1.aa  >& $nam.cd2a.log

  $cdhit2dapp -d 0 -i $nam.ok.aa -i2 $nam.poorcd1.aa -o ${nam}.poor_cd.aa >& $nam.cd2.log

  /bin/rm $nam.*.bak.clstr $nam.poorcd1.aa*
  # log.cd* ?
  
  cat ${nam}.{ok,poor}_cd.aa > ${nam}_allcd.aa
  grep '^>'  ${nam}_allcd.aa | sed 's/>//; s/ .*//;' > ${nam}_allcd.ids
  $evigene/scripts/prot/aaqual.sh  ${nam}_allcd.aa

  echo "# aa-quality for ${nam}_allcd.aa : top 1k summary"
  cat ${nam}_allcd.aa.qual | egrep -v '^#|^total' | sort -k2,2nr | head -1000 | env nam=$nam perl -ne \
'($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; 
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1]; 
print "#$ENV{nam}\t  n=$n; aw=$aw; med=$md; min,max=$mi,$mx; sw=$sw; sn=$sn,$an\n"; }'

  TCAT=cat; 
  fan=`echo $faz | sed 's/\.gz//'`; if [ $faz != $fan ]; then TCAT="gunzip -c"; fi
  $TCAT $faz | env idf=${nam}_allcd.ids perl -ne \
'BEGIN{open(F,$ENV{idf});while(<F>){chomp; $ok{$_}=1 if(/\w/);}}if(/^>(\S+)/){$ok=$ok{$1}} print if $ok;' \
    > ${nam}_allcd.tr

  # gzip --fast ${nam}*.aa ${nam}_allcd.tr
  
} done


