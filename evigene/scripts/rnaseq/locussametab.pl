#!/usr/bin/perl
# sameloctab.pl

use strict;
use Getopt::Long;

my $MINCDS=20;
my $MINALN=85;

my($eqgene,$alntab,$diffloci,$idtable); # = @ARGV;
$eqgene="kf2pub11_both.all.eqgene";
$alntab="kf2pub11_both.alnsense.tab";
$diffloci="kf2pub11_both.diffloci.eqalt2main2";
$idtable="kf2pub11_both.newids";

my $optok= GetOptions(
  "eqgene=s", \$eqgene, 
  "alntab=s", \$alntab,  
  "diffloci=s", \$diffloci, 
  "idtable=s", \$idtable, 
);

my(%okmap,%poormap,%difflocus,%aasize,%aaref,%garef);

open(F,$alntab) or die $alntab;
while(<F>) {
  my($id,$cov,$pid)=(split"\t")[3,7,8]; my $coi=$cov*$pid/100;  
  $okmap{$id}++ if($coi>=$MINALN);  
  # $poormap{$id}++ if($coi<$MINALN and not $okmap{$id});
  # print "$id\t\n" if($coi<$MINALN and not $okmap{$id});
} close(F);

open(F,$diffloci) or die $diffloci;
while(<F>) {
  my($id,$oid,$othrd,$loc,$cla)=split"\t";
  $difflocus{$id}=$cla if($cla =~ /Df\d/);
} close(F);

open(F,$idtable) or die $idtable;
while(<F>) { next unless(/^\w/); 
  my @v=split"\t"; my($id,$gd,$aq)=@v[0,2,5]; 
  my ($aw)= $aq=~m/(\d+)/;  $aasize{$id}=$aw; 
  my ($arv,$ar)= m/aaref:(\d+),([:\w]+)/;  
  if($ar) { $aaref{$id}=$ar; $garef{$gd}{$ar}+=$arv; } 
} close(F);


open(F,$eqgene) or die $eqgene;
while(<F>) {
  my ($td,$od,$ad,$loc)=split"\t"; (my $tg=$td)=~s/t\d+$//; 
  next unless( $okmap{$td} and not $difflocus{$td} );
  my @ad= grep{ not /^$tg/ } split",",$ad; 
  my @sd= grep /\w/, map{ my($d,$p)=split"/"; (($p=~/^[CI]/ or $p>=$MINCDS))?$d:""; } @ad; 
  @sd= grep { $okmap{$_} and not $difflocus{$_} } @sd;
  foreach my $sd (@sd) { 
    my($md,$ad)= sort ($td,$sd); 
    my($clsame,@sameval)= aaorder($md,$ad);
    print join("\t",$clsame,@sameval)."\n"; # if($clsame =~ /^same/);
    } 
} close(F);

# warn "#stats: ...\n";

sub aaorder {
  my($td,$sd)= @_;
  my($tr,$sr,$trc,$src,$taw,$saw,$v,$diffr,$swap);
  my($tg,$sg)=map{ ($v=$_)=~s/t\d+$//; $v; } ($td,$sd); 
  $tr=$aaref{$td}; $sr=$aaref{$sd}; 
  $trc=$garef{$tg}{$tr}; $src=$garef{$sg}{$sr}; 
  $taw=$aasize{$td}; $saw=$aasize{$sd};
  $diffr=($tr and $sr and $tr ne $sr)?1:0;  
  $swap=($src>$trc)?1:($src<$trc)?0:($saw >$taw)?1:0; 
  
  my @val=($sd,$td,$saw,$taw,"$sr,$src","$tr,$trc");
  @val=@val[1,0,3,2,5,4] unless($swap);  
  if($diffr) { return ("notsamer",@val);  } 
  elsif($swap){ return ("samelocsw",@val); } 
  else { return("sameloc",@val); } 
}


__END__

cat $nam.all.eqgene | perl -ne 'BEGIN{$MINCDS=20;}\
($td,$od,$ad,$loc)=split"\t"; ($tg=$td)=~s/t\d+$//; @ad=grep{ not/^$tg/ } split",",$ad; \
@sd=grep /\w/, map{ ($d,$p)=split"/"; ($p=~/^[CI]/ or $p>=$MINCDS)?$d:""; } @ad; \
foreach $s (@sd) { @d=sort ($td,$s); print join("\t","samelocus",@d)."\n"; }  ' | sort -u \
 > $nam.sameloci.eqall

cat kf2pub11_both.alnsense.tab | perl -ne'($id,$cov,$pid)=(split)[3,7,8]; $coi=$cov*$pid/100; \
$ok{$id}++ if($coi>85);  print "$id\t\n" if($coi<85 and not $ok{$id});' | sort -u | \
 ggrep -v -F -f - kf2pub11_both.sameloci.eqall > kf2pub11_both.sameloci.eqall.goodmap

grep Df0 kf2pub11_both.diffloci.eqalt2main2 | cut -f1 | sed 's/^/samelocus  /; s/$/ /;' | \
 ggrep -v -F -f - kf2pub11_both.sameloci.eqall.goodmap > kf2pub11_both.sameloci.eqall.goodnodf

## FIXME correct sameloci table w/ aasize,aaref : dont join tr of diff aaref ?
cat kf2pub11_both.newids kf2pub11_both.sameloci.eqall.goodnodf | perl -ne\
'if(/^samelocus/) { ($ssl,$td,$sd)=split; ($tg,$sg)=map{ \
($v=$_)=~s/t\d+$//; $v; } ($td,$sd); $tr=$ar{$td}; $sr=$ar{$sd}; \
$trc=$gar{$tg}{$tr}; $src=$gar{$sg}{$sr}; $taw=$aw{$td}; $saw=$aw{$sd}; \
$diffr=($tr and $sr and $tr ne $sr)?1:0;  $swap=($src>$trc)?1:($src<$trc)?0: \
($saw >$taw)?1:0; @val=($sd,$td,$saw,$taw,"$sr,$src","$tr,$trc"); \
@val=@val[1,0,3,2,5,4] unless($swap);  if($diffr) { print \
join("\t","notsamer",@val)."\n";  } elsif($swap){ print \
join("\t","samelocsw",@val)."\n"; } else { print \
join("\t","sameloc",@val)."\n"; } } elsif(/^\w/) { @v=split"\t"; \
($id,$gd,$aq)=@v[0,2,5]; ($aw)=$aq=~m/(\d+)/; ($arv,$ar)=m/aaref:(\d+),(\w+)/; \
if($ar) { $ar{$id}=$ar; $gar{$gd}{$ar}+=$arv; }  $aw{$id}=$aw; }' \
> kf2pub11_both.sameloci.swaps.good