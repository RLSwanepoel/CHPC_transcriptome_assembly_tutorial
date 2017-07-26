#!/bin/tcsh
# makeblastscore.sh

# assume blastp is score sorted
set suf=bltab
set evigene=/bio/bio-grid/mb/evigene/

## .. id pattern to skip if any, eg self genes
# set skipho="_ACYPI"

#env blz=bp1-plant8*.blastp.gz makeblastscore.sh  

## fixme: ncbi changes tab to '#' in >defline now, as in
## .. 'ricco:29682.m000589#calmodulin'

foreach bl ( $blz )
 ##set bs=`echo $bl | sed 's/bp1-/bs3-/; s/arp5hum/self/; s/apishymarp/self/'`
 set bs=`echo $bl | sed 's/bp1-/bs3-/; s/plant8/self/; s/arp5hum/self/; s/apishymarp/self/;'`
 set btab=`echo $bl | sed 's/blastp.gz/bltab/'`
 if( -f $btab) continue;
 # echo evigene/scripts/makeblastscore.pl $bl $bs TO  $btab
 # env skipho=$skipho \
 $evigene/scripts/makeblastscore.pl $bl $bs > $btab
end

exit

#cleanids# perl -pi.old  -e's/\#\S+//;' *.bltab

# cat *.bltab | cut -f1,3 | grep -v '     na' | perl -pe 's|homolog=||; s|/[0-9\.]*||; s|,|       |;' | sort -k3,3 -k2,2nr -k1,1 | more

# cat *.bltab | cut -f1,3 | grep -v '     na' | perl -pe 's|homolog=||; s|/[0-9\.]*||; s|,|       |;' | sort -k3,3 -k2,2nr -k1,1 | perl -ne'chomp; ($g,$v,$p)=split"\t";  if($lp eq $p) { $ag.="$g," if($v >= $lv);  } else { putp(); $ag=""; $lg=$g; $lv=$v; } $lp=$p;  sub putp{ print join("\t",$lp,$lv,$lg,$ag)."\n" if($lp); } ' | more

