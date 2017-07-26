#!/bin/bash
# gff2aligntab.sh : tabulate alignment annots from gmap.gff and gsplign.gff
# ** gff2aligntab.pl replaces this

echo "# NOTE: use evigene/scripts/ests/gff2aligntab.pl instead."

#..............
# $nam.align.tab columns
# GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon 
#	splice  aalen   offs    aamap   sense   oid     tag
# 
# cat $nam.align.tab | cut -f1-4,7,8,12,18 | grep -v NOPATH | head
# GenomeID        gespan  geor    AQueryID        qlen    cov     nexon   oid
# Scaffold0       5067-10366      -       Funhe2Exy3m032549t1     679     100.0   3       Funhe2Emap3m022605t1
# Scaffold0       25351-26047     -       Funhe2Exy3m069279t1     697     100     1       Funhe2Emap3m032045t1
# Scaffold0       25660-58049     +       Funhe2Exy3m002726t6     5016    88      25      Funhe2E6bm002559t3

## Add Split gene handling, assume Split=1/2 mRNA follow each other.  Split and path same column?
## Fixme; mRNA.Split=2 has "^#s." rowprefix, sometimes
## gmap orig: Funhe2Exx3m073081t1_C2;trg=Funhe2Exx3m073081t1 597 651;path=2/48;chimera=break..;chim1=Scaffold10054:400380-400949:.;

opts="dosplit=1"; #?
GTAG="NoneSuch"; if [ "X$gtag" != "X" ]; then GTAG=$gtag; fi # eg: xxx-kfish2b

ingff=$*
for az in $ingff; do {
  TCAT=cat;
  nogz=`echo $az | sed 's/.gz//;'`; if [ $az != $nogz ]; then TCAT="gunzip -c"; fi
  nam=`echo $az | sed "s/.gz//; s/\.gff//; s/\.gmap//; s/.$GTAG//;"`;

  $TCAT $az | env $opts perl -ne \
'BEGIN{ @k=qw(match qlen cov pid path indels nexon splice aalen offs aamap sense Split oid tag); 
@gs=qw(gescore clen gaps chimera); @hd=grep{ not/Split/ } @k;
print join("\t","GenomeID","gespan","geor","AQueryID","quspan",@hd)."\n"; } 
chomp; s/^#s\.//; @v=split"\t"; ($gid,$src,$typ,$gb,$ge,$gv,$go,$at)=@v[0,1,2,3,4,5,6,-1];
if($typ eq "CDS") { $cw=1+$ge-$gb; $cdsw+=$cw; }
elsif($typ eq "intron"){ $inspl += ($gv>69)?2:($gv>44)?1:0; }
elsif($typ eq "mRNA") { $issplit=(/Split=|chimera=/)?1:0; 
($tid,$tb,$te)=m/(?:trg|Target)=(\S+) (\d+) (\d+)/; 
putg() if($ltid); $inspl=$cdsw=0; 
%at= map{ $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $_ => $v; } ("ID",@k);
$id=$at{ID}; $tid ||= $id; $at{Split}=~s/^C//; 
%ags= map{ if($at=~m/\b$_=([^;\s]+)/) { $_=>$1; } } @gs;
if(/gescore=/) {  $tag="gspl";
  if($at{splice}) { $at{splice}= 2 + int($at{splice}/2); }
  $at{qlen}= $qlen = $ags{clen}||0; $gaps=$ags{gaps}||0;
  ($pcov,$match,$ql2)= $at{cov} =~ m/(\d+).,(\d+).(\d+)/;
  $at{pid}||=99; $at{cov} = $pcov; $at{match}= $match;
} else { $tag="gmap";
  $aamap= $at{aamap}||$at{aalen}||0;
  ($chi)= $id =~ m/_C(\d)$/; if($ags{chimera} or $chi) { $chi||=1; $at{Split}="$chi/2" unless($at{Split}); }
  if($at =~ /aalen=(\d+,\d+[^;\s]+)/) { $aaq=$1; $at{aalen}=$aaq; $at{aamap}=$aamap if($aamap ne $aaq); }
} 
$at{tag}=$tag; ($ltid,$ltb,$lte)=($tid,$tb,$te);
($lgid,$lgb,$lge,$lgv,$lgo,$lat)=($gid,$gb,$ge,$gv,$go,$at);
} END{ $issplit=0; putg(); }
sub putg { 
if($issplit and $ltid eq $tid) { %lat=%at; $cov1=$at{cov}; 
if($gid ne $lgid) { $lat{Split}="C3:$lgid:$lgb-$lge,$lgo,$cov1"; }
elsif($lge < $gb-1000 or $lgb > $ge+1000) { $lat{Split}="C2:$lgb-$lge,$lgo,$cov1"; } 
else { $lat{Split}="C1:$lgb-$lge,$lgo,$cov1"; } return; }
elsif($at{Split} and $lat{Split}) { %att=%lat; 
  map{ $att{$_}+=$at{$_}; } qw(cov match nexon aamap); %at=%att; %lat=(); } 
if($at{Split}) { $at{path}= $at{Split}; }
$at{splice}= 2+$inspl if($inspl>0 and not $at{splice});
$at{aamap}= int($cdsw/3) if($cdsw>0 and not $at{aamap});
@at= @at{@hd}; print join("\t",$lgid,"$lgb-$lge",$lgo,$ltid,"$ltb-$lte",@at)."\n"; } ' \
 | sort -k1,1 -k2,2n -k4,4 > $nam.align.tab

} done

