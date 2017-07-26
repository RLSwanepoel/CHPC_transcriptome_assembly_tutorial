#!/bin/tcsh

foreach samz (sams/*.gsnap.sam.gz)
  set gp=`basename $samz .gsnap.sam.gz`
  echo "# Make velvparts/$gp.pairo.s01361.fa"

  gzgrep scaffold01361 sams/$gp.gsnap.sam.gz | perl -ne \
'($d,$f,$c,$cb,$xf,$cig,$mc,$mb,$mx,$sq)=split"\t"; \
$rv=($f & 0x0010)?1:2; $dr="$d/$rv";  unless($did{$dr}++) { \
$o="u"; $vcig=""; if(m/XS:A:(.)/){ $o=$ors{$1}||"u"; $vcig=" $cig"; } \
BEGIN{ %ors=("+"=>"f","-"=>"r"); } \
print ">$dr/$o $c:$cb:$o $mc:$mb$vcig\t$sq\n"; }' \
| sort -k1,1 | perl -pe's/\t/\n/;' > velvparts/$gp.pairo.s01361.fa0

# patch in strand to u-mate
 cat velvparts/$gp.pairo.s01361.fa0 | perl -ne\
'if(/^>(\S+)/){ $lh=$h; $lo=$o; $h=$1; $ld=$d; $d=$_; \
$h=~s,/./(.)$,,; $o=$1; if($h eq $lh) { if($o ne $lo){ \
if($lo eq "u") { $ld=~s,/u,/$o,; } elsif($o eq "u"){ $d=~s,/u,/$lo,;} } }  \
print $ld,$la if $ld;  } else { $la=$_;  } END{ print $d,$la;}' \
>  velvparts/$gp.pairo.s01361.fa

end

