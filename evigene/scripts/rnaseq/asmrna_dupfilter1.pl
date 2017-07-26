#!/usr/bin/perl
# asmrna_dupfilter.pl

=item about
  
  Filter out hi identity duplicate transcripts from transcript assembly,
  using megablast AFTER identifying best proteins.
  
  Principle is that over-assembling transcript reads with many assembly options
  produces a subset of accurate assemblies in a superset of crappy ones.
  
  Best-protein detection is done on superset (longest orf + cds/utr quality),
  then clustered (cd-hit) to reduce to "best" subset by protein quality.
  
  However cd-hit protein filtering retains nearly identical transcripts assembled
  from same reads, but differ in protein enough to avoid cd-hit cluster.
  
  This step identifies high identity transcripts from bestcd subset, from megablast align,
  then marks/removes the hi-id fraction, keeping longest-protein/longest-tr of each
  hi-id cluster.

=item usage

 $evigene/scripts/rnaseq/asmrna_dupfilter.pl -aa=$pt.aa.count -tr=$pt.tr.count -blast=$pt.mblast\
    > $pt.dupclust


=item inputs

  inputs: 
    protein sizes, transcript sizes, as faCount table (id, size); see aacount for gappy proteins

   tr.megablast output from:
   makeblastdb -in $pt.tr -dbtype nucl -out $pt
   blastn -db $pt -query $pt.tr -task megablast -evalue 1e-19 -perc_identity 99 -dust no -outfmt 7 -out $pt.mblast

  outputs: 
    table of same/subset transcripts, rough equiv of cd-hit-est
    bestset.tr, filtered: 1.bestsetaa.ids only, 2.remove trsame subset (like cd-hit-est, but diff methods)


## FIXME asmrna_dupfilter.pl: check for ID mismatches == too many zeros in aa.count /tr.count / blastn
## eg. tables have diff prefixes like 'litova:vel' vs 'litovavel'

=item revise output

  FIXME: make OUTSPANTAB=1 default, see below alnclass.sh as better/faster classifier.
      -- newer class assignment: cdsidentclass.sh, using %ident (>98) and %align (>50%)
        i.e., shorter alignments of identity are tested as valid criteria for alternate-transcripts
              using genome-mapped transcripts, hi-ident + part-align gives *mostly* alt-tr, some paralogs,
              depending on species & freq of hi-ident paralogs.
              
  update: fold in tested methods for tr class assignment, to replace outclusters
  
=item updated for blat psl input

  blat -out=psl -t=dna -minIdentity=98 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat
  blat -out=psl -t=dna  -minScore=99  -minIdentity=99 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat99
  #  -minScore=99 or such to reduce fragment aligns, 30 default, nmatch - nmiss

=item see also
  
  evigene/scripts/prot/aacount.sh
  evigene/scripts/makeblastscore.pl
  evigene/scripts/rnaseq/asmrna_equalgenes.pl
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use warnings;
use Getopt::Long;

my $debug= 1;
my $OUTSPANTAB= 0;  # make default until replace outclusters()

my $MINALIGN= 90; # change to several levels
my $TINYALN = 50; # ignore less than this
my $MINSIZE = 999; # not used yet
 
my $OK_CDSUTR= 60;
my $BAD_CDSUTR= 30;
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

# my $IDPREFIX= $ENV{idprefix} || ""; # refgene id pattern; need for aaclusters at least,    
my ($aasizes,$trsizes,$blatpsl,$blastab,$outeqtab,$trokids,$trpoorids,$logfile,$head);

my $optok= GetOptions(
  "aasizes=s", \$aasizes, 
  "trsizes=s", \$trsizes, 
  "blastab=s", \$blastab,  
  "blat=s", \$blatpsl, # other input format opts here..?  -informat=xxx
  "MINALIGN=i", \$MINALIGN,  
  "MINSIZE=i", \$MINSIZE,  
  "CDSUTR=i", \$OK_CDSUTR,  
  "OUTSPANTAB!", \$OUTSPANTAB, 
  "debug!", \$debug, 
#  "outeqtab=s", \$outeqtab,   
#  "trokids|trgoodids=s", \$trokids,   
#  "trpoorids|trbadids=s", \$trpoorids,  
#  "logfile=s", \$logfile,  
#  "geneidprefix=s", \$IDPREFIX,  
  );

die "usage:  asmrna_dupfilter.pl -aasize=name.aa.count -trsize=name.tr.count -blast=name.mblast | -blat=name.blatpsl
 opts: -MINALIGN=$MINALIGN -CDSUTR=$OK_CDSUTR percents  -OUTSPANTAB
" unless($optok and $aasizes and $trsizes and ($blatpsl or $blastab));

if($MINALIGN>0) {
  $MINALIGN= 100*$MINALIGN if($MINALIGN<1);
  $TINYALN=$MINALIGN if($TINYALN>$MINALIGN);
}

my(%aasize, %trsize);

if($aasizes) {  
  ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
  if($aasizes =~ /count|qual/) { open(F,$aasizes) or die "FAIL: read $aasizes ..."; }
  else { open(F,"faCount $aasizes |") or die "FAIL: faCount $aasizes ..."; }
  while(<F>) { next if(/^#/ or /^total/); my($id,$al,@ax)=split; 
    ## if(@ax==1) { $al -= $ax[0]; } # gap count
    $aasize{$id}=$al; } close(F); 
}
if($trsizes =~ /^aasize/) {
  foreach my $id (keys %aasize) { $trsize{$id}= 3*$aasize{$id}; }
} elsif($trsizes) {  
  if($trsizes =~ /count|qual/) { open(F,$trsizes) or die "FAIL: read $trsizes ..."; }
  else { open(F,"faCount $trsizes |") or die "FAIL: faCount  $trsizes ..."; }
  while(<F>) { next if(/^#/ or /^total/); my($id,$al)=split; $trsize{$id}=$al; } close(F); 
}

my(%better, %outrows, %outids, %bspans, %didpair);

if($blatpsl) { readblatpsl($blatpsl); } # can detect from input table ...
elsif($blastab) { readblasttab($blastab); }

outclusters() unless($OUTSPANTAB);

#.................................

sub outclusters {

  # $better{$lq}{$lt}++;
  my @ids= sort keys %outids;
  my %sumscore;
  foreach my $id (@ids) { 
    foreach my $jd (@ids) { next if($id eq $jd); 
    my $v= $better{$id}{$jd}||0; $sumscore{$id} += $v;
    }
  }
  
  ## need header:
  # outrow: join("\t",$typ,$lq,"$aq/$wq",$lt,"$at/$wt",$sa,$diffaln,$sm)."\n"; 
  print join("\t",qw(Cluster Score12 Type Qid Qsize Tid Tsize Align Daln Ident))."\n" unless($head++);
      
  my %didid=(); my $cluid=0;
  my @topid= sort{ $sumscore{$b} <=> $sumscore{$a} } @ids;
  foreach my $topid (@topid) { 
    next if($didid{$topid});
    $cluid++; $didid{$topid}++;  
    my $topscore= $sumscore{$topid} || 0; 
    my @nextid=  sort{ $better{$topid}{$b} <=> $better{$topid}{$a} } keys %{$better{$topid}};
    foreach my $nextid (@nextid) {
      next if($didid{$nextid});
      my $nextscore= $sumscore{$nextid} || 0; 
      # my $bscore= $better{$topid}{$nextid} || 0;
      my $outrow= $outrows{$topid}{$nextid};
      unless($outrow) { 
        if($outrow= $outrows{$nextid}{$topid}) { 
          my @r=split"\t",$outrow;
          @r[1,2,3,4]= @r[3,4,1,2]; 
          my $ty=$r[0]; unless($r[0] =~ s/qsub/tsub/) { $r[0] =~ s/tsub/qsub/; }
          $outrow= join"\t",@r;
        }
      }
      $outrow ||= "na\n"; # outrows has \n
      print "cl$cluid\t$topscore,$nextscore\t",$outrow;
      $didid{$nextid}++;
    }
  }
}

# need to do more before output: cluster all same/subset by ids
sub puts { 
  my($lq,$lt,$sa,$sident)= @_;
  return if($lq eq $lt);
  # $sa -= $sm #?? adjust align by mismatches? No
  my($qsame,$tsame,$diffaln,$adiffqt)= (0) x 10;
  my $typ=""; 
  my $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
  my $at= $aasize{$lt}||0; $at *=3;
  my $wq= $trsize{$lq}||0; 
  my $wt= $trsize{$lt}||0;  
  
  map{ $typ="missingsize" unless($_ > 0); } ($aq,$at,$wq,$wt);
  unless($typ) {
    my $paq= 100*$sa/$wq;
    my $pat= 100*$sa/$wt;
    my $cdq= 100*$aq/$wq;
    my $cdt= 100*$at/$wt;
    
    if($paq < $TINYALN and $pat < $TINYALN) { 
    
    } else {
    
    # add more classes here for pALIGN: sim90, sim80, sim70, ..
    # esp. for tiny fragments that mostly fit in bigger tr of good aa qual
    # need also class crappy big-tr of bad utr
    
    # $qsame= $paq >= $MINALIGN; 
    # $tsame= $pat >= $MINALIGN;   
  
    foreach my $pa (90,80,70,60) { # TEST; 
      # this doubles ndups, most at 80% level, but dont know yet if these are useless subset trasm or valid alts/paralogs
      if($qsame==0 and $paq >= $pa) { $qsame=$pa; }
      if($tsame==0 and $pat >= $pa) { $tsame=$pa; }
    }
  
    $diffaln= ($wt > $wq) ? $wq - $sa : $wt - $sa;
    $adiffqt = $aq - $at;
  
    #? check also wcds/wtr for utrbad score?
    # add class: when pCDS < 10%-33% ?, q,t/crappylong, NOT better .. want to trash these cases
    my $qutrbad= ($cdq >= $OK_CDSUTR)?0: ($cdq > $BAD_CDSUTR)?1: 2;
    my $tutrbad= ($cdt >= $OK_CDSUTR)?0: ($cdt > $BAD_CDSUTR)?1: 2;
    my $utrdiffqt= ($qutrbad and not $tutrbad)? -1 : (not $qutrbad and $tutrbad)? +1 : 0;
    
    if($qsame and $tsame) { 
    
      $typ= ($adiffqt==0 and $wq==$wt) ? "same$tsame"
      : ($adiffqt>0) ? "tsubset$tsame" 
      : ($adiffqt<0) ? "qsubset$qsame" 
      : ($utrdiffqt<0)? "qsubset$qsame"  # qutrbad==2 : "qsubsetbad"
      : ($utrdiffqt>0)? "tsubset$tsame"  # tutrbad==2 : "tsubsetbad"
      : ($wq<$wt)? "qsubset$qsame"
      : ($wq>$wt)? "tsubset$tsame"
      : "same$tsame"; 
    }
    ## add bad == 2 class here, for qsame but tutrbad==2, t is bad; vv
    elsif($qsame) { $typ= ($adiffqt <= 0) ? "qsubset$qsame" : "tsubsetlong$qsame"; } #? tutrbad
    elsif($tsame) { $typ= ($adiffqt >= 0) ? "tsubset$tsame" : "qsubsetlong$tsame"; } #? qutrbad
    }
  }
  
  # unless($head) { print join("\t",qw(Type Qid Qsize Tid Tsize Align Daln Ident))."\n"; $head++; }
  if($typ) { 
    my $val= join("\t",$typ,$lq,"$aq/$wq",$lt,"$at/$wt",$sa,$diffaln,$sident)."\n"; 
    # need clustering hash here to pick out best for many q,t subsets
    if($typ =~ /missing/) { $better{$lq}{$lt}=-999999; }
    elsif($typ =~ /tsubset|same/) { $better{$lq}{$lt}++; }
    elsif($typ =~ /qsubset/) {  $better{$lt}{$lq}++; }
    $outids{$lq}++; $outids{$lt}++;
    $outrows{$lq}{$lt} = $val; # if exists $outrows{$lq}{$lt}, *should* be same val, check?
    #NO# print $val;
    } 
}

use constant SPANSUM => 1;

sub putspans {
  my($lq)= @_;
  foreach my $lt (sort keys %bspans) {
    my @sp= @{$bspans{$lt}}; 
    my ($tbit,$taln,$tidn)= (0) x 9;
    foreach my $sp (@sp) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident)= @$sp;
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident;
      }
    # my $mis= $taln - $tidn; # dont need this; replace w/ tidn
    if($OUTSPANTAB) { 
      # add Qsize,Tsize  cds/trlen ?
      my $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      my $at= $aasize{$lt}||0; $at *=3;
      my $wq= $trsize{$lq}||0; 
      my $wt= $trsize{$lt}||0;  
      next if( (100*$taln / _max(1,_min($aq,$at)) ) < $TINYALN); # skip trival matches
      print join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
      } 
    else { puts($lq,$lt,$taln,$tidn); }
  } 
  %bspans=();
}

sub readblasttab {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my $fh;
  if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open(F,$bother) or die $bother; $fh=*F; }
  %bspans=();
  
  while(<$fh>) { 
    unless(/^\w/) { 
    #  if(/^# Query/) { ## Unused info
    #  #Notnow# puts($lq,$lt,$sa,$sm) if($lt); 
    #  ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; }
      next; } 

    my @v= split; 
    my($q,$t,$bits,$aln,$mis,$gap,@bspan)= @v[0,1,-1,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 

    # FIXME: need to parse align parts before summing;
if(SPANSUM) {
    if($lq and $q ne $lq) {
      putspans($lq); %bspans=();
    }
}    
    ## UNUSED now: qd, wq
    if($t eq $q) { $qd=$q unless($qd); $wq=$aln unless($wq); }
    else { 

if(SPANSUM) {
      my $aident= _max(0,$aln-$mis-$gap); # other way to calc: $aident = $pctident * $aln;
      sumblastpart( $q, $t, $bits,$aln,$aident, @bspan); # if($bits > $pMINLOW * $bmax or $bspans{$t}); 
} else {      
      if($t eq $lt) {
        # sumscore( $q, $t, $bits,$aln,$aident, @bspan); 
         $sa += $aln; $sm += $mis+$gap;  
      } else { 
        my $aident= _max(0,$sa-$sm); # other way to calc: $aident = $pctident * $aln;
        puts($lq,$lt,$sa,$aident) if($lt);
        $sa=$aln; $sm=$mis+$gap;  
      } 
}       
    } 
    $lq=$q; $lt=$t;  
  } close($fh);
  
if(SPANSUM) {
  putspans($lq);
} else {
  my $aident= _max(0,$sa-$sm); # other way to calc: $aident = $pctident * $aln;
  puts($lq,$lt,$sa,$aident) if($lt); 
}  
}

sub readblatpsl {
  my ($bother)= @_;
  
  my $ispsl=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my $fh;
  if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  $ispsl=1 ;
  while(<$fh>) { 
    ## $ispsl=1 if m/^psLayout/; # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    if(@v==22) { shift(@v); } # ucsc psl starts with a 'bin' field?
    unless(@v==21 and $ispsl){ die "# error: doesnt look like psl format I know" ; }
    my( $mat, $mis, $rep_matches, 
      $qgapw, $tgapw, $orient,
      $qid, $qsize, $qstart, $qend,
      $tid, $tsize, $tstart, $tend,
      $blocksizes, $qstarts, $tstarts,
      )= @v[0..2, 5,7, 8..16, 18..20];
    $qstart++; $tstart++; # move to 1-origin

## FIXME: allow for blat split rows per query, happens sometimes.., use bspans as for blast
    if($lq and $qid ne $lq) {  
      putspans($lq); %bspans=();
    }

    if($tid eq $qid) { 
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      # my $pmat= 100 * $mat / $qsize;
      # next if($pmat < $TINYALN); # NOT NOW with bspan parts
      my $aident= _max(0,$mat-$mis);  # not sure -gap right there, align should be +gap, ident = match
      my $aln= $mat + $qgapw; # is this right? align span includes query gaps, and mismatches
      my $bits=  $mat; # what? maybe 2*mat ? 2*aident? 
      #** should split blocks in qstarts, tstarts to sum up q/t start,end
      ## but dont need see putspans:       $tbit += $xbit; $taln+= $aln; $tidn+= $aident;
      
      $bspans{$tid}=[] unless(defined $bspans{$tid}); ## [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]; 
      # putspans($qid);  # only 1 pairmatch/span per psl row.
      push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); 
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq);   
}





sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
sub sumblastpart {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0;
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }
  unless($bspans{$t}) { 
    $bspans{$t}=[]; push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident]); 
    return; }
  my $ov=0;
  ## 2011oct overslop pct fix
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  foreach my $sp (@{$bspans{$t}}) {
    my($xb,$xe,$tb,$te,$xbit)= @$sp;
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
  ## add 2011.aug
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident]); }
}





__END__


=item newer cdsidentclass.sh  2013.02.13

if [ "X" = "X$phi" ]; then phi=99; fi   #below 99 gives many more misclass diff locus
if [ "X" = "X$pmid" ]; then pmid=95; fi
if [ "X" = "X$plow" ]; then plow=80; fi
if [ "X" = "X$mina" ]; then mina=50; fi
if [ "X" = "X$cdsalign" ]; then cdsalign=1; fi
if [ "X" = "X$suf" ]; then suf=clid2tab; fi

  cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n |\
  env mina=$mina cdsw=$cdsalign phi=$phi pmid=$pmid plow=$plow  perl -ne \
'BEGIN{ 
$MINAL=$ENV{mina}||50; $CDSW=$ENV{cdsw}; $PHIALN=65; $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; }
chomp; @v=split; ($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)=@v; $isfrag= $aclass="";
if($tc > $qc) { ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } #?? always swap to greater cds?
# if($class{$qd} or $class{$td}) { next; } # maybe not. td-noclass > qd-main other
if($class{$td}) { next; }  
if($CDSW) { $ww=($qc>$tc and $tc>0)?$tc:($qc>0)?$qc:$tc;  $isfrag= ($tc < 0.5*$qc)?"frag":""; } 
else { $ww=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw;  $isfrag= ($tw < 0.5*$qw)?"frag":"";  }
$pid= ($aln<1)?0: int(100*(0.5+$iden)/$aln); $pid=100 if($pid>100);
$pal= ($ww<1)?0 : int(100*(0.5+$aln)/$ww); $pal=100 if($pal>100);
$bestmatch{$qd}="$td,$pid/$pal" unless($bestmatch{$qd});
$bestmatch{$td}="$qd,$pid/$pal" unless($bestmatch{$td});
if($pal < $MINAL) { } #defer: $aclass="noalign$isfrag";  
elsif($tc == 0 and $qc > 0 and $pid >= $PLOW) { $aclass="frag0aa"; } # tc==0 special case "frag0"
elsif( $pid >= $PHI ) { $aclass= ($pal<$PHIALN or $isfrag)?"parthi":"althi"; }
## elsif( $pid >= $PHI ) { $aclass="althi$isfrag"; }
elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }
if($aclass) { $class{$td}=$aclass; $bestmatch{$td}="$qd,$pid/$pal"; $ismain{$qd}++; }
END { foreach $d (sort keys %class) {
 ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; print join("\t",$d,$c,$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"main",$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"noclass",$q,$pal)."\n"; }
}' > $nam.$suf


=item redo using -outspantab (quick vs slow for all above)

$evigene/scripts/rnaseq/asmrna_dupfilter.pl -outspan \
 -aa=$pt.aa.count -tr=$pt.tr.count -blast=sd-slf95-$pt.blastn.gz > $pt-s95.alntab

head locust1best5-s95.alntab
Qid                             Qclen   Qtlen   Tid                        Tclen   Ttlen   Align   Ident   Bits
locust1sop4p3k31loc9483t1       369     370     locust1Svel1K31L7036t1     1857    1858    370     368     673
locust1sop4p3k31loc9483t1       369     370     locust1sop4p1k23loc1146t2  1152    4185    370     369     678
locust1sop4p3k31loc9483t1       369     370     locust1sop4p2k23loc30t2    867     3971    370     368     673
locust1sop4p3k31loc9495t1       1632    4114    locust1sop4p0k31loc11571t1 288     687     307     306     563
locust1sop4p3k31loc9495t1       1632    4114    locust1sop4p1k31loc7134t3  3231    5461    2412    2406    4355

# pick longest aa, class those aligning well as (a) alts, (f) fragments (add utrbad class?)
# alts : aa >=60% of 1st, align < 100%
# frags: aa < 60% of 1st, align > 95% min(slen,tlength) : min to drop utrbad class

# FIXME: .alntab from blastn is missing any No-blast-match tr (or just self-match, which should be there).
# .. add nomatch to alntab via asmrna_dupfilter ?

#.............................................
# TEST: need to validate blastn align classes w/ genome-mapped trasm,
#    how many of each class are true/false same-locus,
#    vs paralog vs artifacts (poor genomap, poor intron agreement, ..)
#  -- cacao, daphmag, maybe fungr test cases; 
# Alternate, test by read-trmap (bowtie), for multi-maps, proper pairs, homogeneous? read cover
# I.e., need way to validate false pos/false neg and true pos for trasm that have some blastn identity.
# See aaeffect data sets, cacao & daphmag
/Users/gilbertd/Desktop/dspp-work/cacao/cacao3d/rnas/genes/
cacao align6f:
  cacao3tri1asmcd.aligny.tab	cacao3estasm_cd.aligny.tab	cacao3sopc11nr.aligny.tab ..
  $caca/rnas/genes/  trgmap, genogmap, introns
dmagalign6f:  classify each trset_allcd w/ asmrna_dupfilter, then count agree/disagree for alts == same-locus, paralog, other
	daphmag2nwb_allcd.aligny.tab	daphmag2vel9h_allcd.aligny.tab daphmag3cuf13th3.aligny.tab  ..
  $dmag/rnas/asmrna3/genes/  trgmap, genomap, introns,
  
# ** PROBLEMS: cacao, daph show 1/2 of these subset aligns are on diff scaffolds.. 
#  hi-id align trs are hard to classify as same/diff locus.  Use blastn mismatch as key? no mismatch?
#  use pident >= 99% as strict same-locus trasm classifier?

#.............................................

# sort by big q-aa, hi-align
# fixme: add utrbad tests, pick utrgood vs bad for (near)same-size-aa ; sort -k6,6n to get smaller Ttlen for same clen/align

#test# cat locust1best5-s95.alntab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n | perl -ne 
#.............................................
#!/bin/bash
# alnclass.sh

alnset=$*
for atab in $alnset; do {
  nam=`basename $atab .alntab`
  if [ ! -f $nam.classtab ]; then 
  echo "# $atab TO $nam.classtab"

# fixme: if td/tw==0 should call subclass frag or missing, NOT set ww=al;
# below patch will class td/tw==0 as frag

  cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n | perl -ne \
'chomp; @v=split; ($qd,$qc,$qw,$td,$tc,$tw,$al)=@v;
if($class{$qd} or $class{$td}) { next; }
$ww=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw; 
$pal=($ww<1)?0:int(100*(0.5+$al)/$ww);
$bestmatch{$qd}="$td,$pal" unless($bestmatch{$qd});
$bestmatch{$td}="$qd,$pal" unless($bestmatch{$td});
$aclass="";
if($tc < 0.5*$qc and $al >= 0.95*$ww) { $aclass="frag"; }
elsif($tc >= 0.8*$qc and $al >= 0.95*$ww ) { $aclass="althi"; }
elsif($tc >= 0.5*$qc and $al >= 0.95*$ww ) { $aclass="althifrag"; }
elsif($tc >= 0.5*$qc and $al >= 0.80*$tw and $al < 0.95*$ww) { $aclass="altmid"; }
elsif($tc >= 0.5*$qc and $al >= 0.80*$ww ) { $aclass="altlo"; }
elsif($tc < 0.5*$qc and $al >= 0.80*$ww ) { $aclass="altlofrag"; }
if($aclass) { $class{$td}=$aclass; $bestmatch{$td}="$qd,$pal"; $ismain{$qd}++; }
END { foreach $d (sort keys %class) {
 ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; print join("\t",$d,$c,$q,$pal)."\n"; }
foreach $d (sort keys %ismain) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"main",$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"noclass",$q,$pal)."\n"; }
}' > $nam.classtab

fi

tabsum=`cut -f2 $nam.classtab | sort | uniq -c | perl -pe 's/\n/, /;'`
echo "# classes $nam: $tabsum"

} done

#.................
# recalc v3:

cut -f2 locust1best5-s95.classtab2 | sort | uniq -c | head
4062 althi
4915 althifrag
1590 altlo
10244 altlofrag
10248 altmid
18662 frag
20007 main
23808 noclass
-----
49721 allsubset
93507 total uniqids (27 are classed as main + altlo/..)
drop: frag, segregate: althi/hifrag, maybe altmid
keep as primary trset: main,noclass,altlo, maybe altmid

# calc v1:
cut -f2 locust1best5-s95.classtab | sort | uniq -c 
  9570 alt
  4127 althi        # keep but mark or segregate alt/althi set
  3231 althifrag    # maybe drop these
   629 altlo        # combine altlo/altlofrag, not distinct enough
  8397 altlofrag    # could be various: paralogs, alts, ..
 20723 frag         # drop these
 -----
 46677 allsubset
 20382 mainclass (has alts + frags)
 62362 unclassed (unique-ish) : add to .classtab ? > n=47089 this way, others?
129421 totalaa

#..recalc
cut -f2 locust1best5-s95.classtab2 | grep -v noclass | sort | uniq -c
 10214 alt   : rename altmid
  4061 althi
  4909 althifrag
  1624 altlo
 10242 altlofrag
 18659 frag
 -----
 49709 allsubset
 29007 mainclass
 23809 unclassed (unique-ish)   
102525 total 

1719 alt,althi have same gene id
 376 frag have same gene id (true alts?)
  42 altlo same gene id

cut -f2 locust1tri1_allcd-s95.classtab2 | grep -v noclass | sort | uniq -c
2588 althi
1090 althifrag
 123 altlo
 662 altlofrag
1990 altmid
1048 frag
----------
7501 allsubset
6882 mainclass
7413 unclassed
21796 total ?? low, 50168 in aa

cat locust1tri1_allcd-s95.classtab2 | egrep '    (frag|althi|alt )' | cut -f1 | \
sed 's/$/        /' | ggrep -v -F -f - locust1tri1_allcd.aa.count > locust1tri1_noalt.aa.qual


#locust1best5.aa          n=1000; aw=2082; med=1722; min,max=1307,14259; sw=2082891; sn=24999,24.9
#locust1best5noalt.aa     n=1000; aw=1784; med=1472; min,max=1144,14259; sw=1784691; sn=23075,23
#... vs ...
#locust1tri1_allcd.aa     n=1000; aw=1367; med=1179; min,max=919,6074; sw=1367083; sn=0,0
#locust1tri1_noalt.aa     n=1000; aw=1331; med=1152; min,max=894,6074; sw=1331403; sn=0,0



Eg.,
grep locust1Svel1K23L10030t  locust1best5.aa.count
locust1Svel1K23L10030t1 64      0       64,99%,partial  193   << frag
locust1Svel1K23L10030t4 795     0       795,95%,partial3        2510  << mainclass, best prot
locust1Svel1K23L10030t8 707     0       707,99%,partial 2122   << 2nd mainclass
locust1Svel1K23L10030t11        421     0       421,40%,complete-utrbad 3102
  ^^ has lowish align to locust1Svel1K23L10030t4, 1083/(2510,3102); probably bad assembly.
  
grep locust1Svel1K23L10030t locust1best5-s95.classtab
locust1Svel1K23L10030t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L26321t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L37361t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L51443t1 frag    locust1Svel1K23L10030t8
locust1sop4p0k31loc57894t1      frag    locust1Svel1K23L10030t8
locust1sop4p4k31loc57872t1      frag    locust1Svel1K23L10030t8
locust1tri1loc1401455c0t1       frag    locust1Svel1K23L10030t8
locust1tri1loc240526c0t1        altlofrag       locust1Svel1K23L10030t8
locust1tri1loc892983c0t1        frag    locust1Svel1K23L10030t4
locust1tri1loc897709c0t1        frag    locust1Svel1K23L10030t4
locust1tri1loc904858c0t1        frag    locust1Svel1K23L10030t4
locust1vel5k35Loc6409t1 frag    locust1Svel1K23L10030t4


=item FIXME for crappy long utrbad big.tr

whitefly1vel5k35Loc1t58595 << "locus" with 100k's of alternates must be crappy bin of all repetitive junk
  tr is 112k long but only 6000 longest orf used ..
  -- should skip such crappy big tr, classify in output list
  -- look for pCDS < 10%?
  
grep whitefly1vel5k35Loc1t58595 ../../../tsaevg/trsets/whitefly*.{tr,aa}.count 
whitefly1best3.tr.count:whitefly1vel5k35Loc1t58595       112482  31960   23096   23417   32074   1935    4896
whitefly1best3.aa.count:whitefly1vel5k35Loc1t58595       1869    3       1872,4%,complete-utrbad 112482
flamingo2.% grep whitefly1vel5k35Loc1t58595  whitefly*.dupclust | wc -l
     110


=item prelim tests

# duptrfilter.sh
# messy now; somewhat akin to cd-hit-est but 
# need to use blast-align instead of clustering, to get same tr w/ minor breaks
# or would cd-hit-est do what is needed? prob not. test?
#   this is not right yet; DROPPED LARGEST prots ***
#   problem: locust1vel5k55Loc697t2 = locust1sop4p1k23loc4964t1 longer tr, much shorter 4000 aa ***

# inputs: bestset.aa or bestsetaa.ids from cd-hit alltr.aa; all-source.tr
# outputs: bestset.tr, filtered: 1.bestsetaa.ids only, 2.remove trsame subset (like cd-hit-est, but diff methods)

ncbin=$bg/mb/ncbicpp/bin

pt=whitefly1best3vel1k_cd90
pt=locust1best51k

# s1: collect tr from bestcds.ids, limit to long tr? or do all as prep for TSA submit (drop crap)
# $d=~s/whiteflyvel/whitefly1vel/; s/whiteflyvel/whitefly1vel/; << do id change elsewhere.

gzcat vel?trs/*.tr.gz | cat bestof3/whitefly1best3vel_cd90.ids - | env idp=whitefly min=999 perl -ne \
'BEGIN{$IDP=$ENV{idp}; $MINTR=$ENV{min}||999;} if(/^($IDP\S+)\s*$/) { $ok{$1}=1 } else { \
if(/^>(\S+)/) { $d=$1; $ok=$ok{$d}; if(m/len=(\d+)/) { $ok=0 if($1<$MINTR); } } print if($ok); }' \
 > whitefly1best3vel1k_cd90.tr

gzcat vel*trs/locust1vel{4,5,u.all2_cd}.tr.gz  tr{soap,in}/locust*all2_cd.tr.gz | \
cat bestof5/locust1vel15trisop_cd90.ids - |  env idp=locust min=999 perl -ne \
...
> locust1vel15trisop_cd90.tr

# s2: megablast tr x tr; want high ident only

$ncbin/makeblastdb -in $pt.tr -dbtype nucl -out $pt
$ncbin/blastn -db $pt -query $pt.tr -evalue 1e-19 -perc_identity 99 -dust no -outfmt 7 -out $pt.mblast

# s3: critical, identify "same" transcripts.
# this needs adjusting to make sametr.tab from mblast
# need sizes of both q,t tr to judge sameness : same if align(tq) ~= sizeof(t) or sizeof(q)

# orig
# cat $pt.mblast | perl -ne'if(/^# Query/) { ($wq)=m/len=(\d+)/; ($qd)=m/Query:\s+(\S+)/; } elsif(/^(\w+)/) { ($q,$t,$p,$al,$mi,$ga,@v)=split;  unless($t eq $q) { if($t eq $lt) { $sa+=$al; $ss.=$_; $sm+=$mi+$ga;} else { if($sa/$wq > 0.9 and $wq/$sa > 0.9) { print join("\t","same",$lq,$lqw,$lt,$sa,$lqw-$sa,$sm)."\n";  print $ss; } $ss=$_; $sa=$al; $sm=$mi+$ga; } } $lq=$q; $lqw=$wq; $lt=$t;}' > $pt.sametab

# add tr.count input; AND aa.count here, pick same/subset cases from aasize > trsize ..
cat $pt.mblast | perl -ne\
'if(/^# Query/) { puts($lq,$lt,$sa,$sm) if($lt); ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; } \
elsif(/^\w/) { ($q,$t,$pi,$al,$mi,$ga,@v)=split; \
if($t eq $q) { $qd=$q unless($qd); $wq=$al unless($wq); } \
else { if($t eq $lt) { $sa += $al; $sm += $mi+$ga; $ss.=$_; } \
 else { puts($lq,$lt,$sa,$sm) if($lt); $sa=$al; $sm=$mi+$ga; $ss=$_;  } } \
$lq=$q; $lt=$t; } END{ puts($lq,$lt,$sa,$sm) if($lt); }
sub puts { my($lq,$lt,$sa,$sm)= @_;
my $wq= $sizes{$lq}||-1; my $wt=$sizes{$lt}||-1;
$qsame=($wq<1)? 0: $sa/$wq > 0.9; 
$tsame=($wt<1)? 0: $sa/$wt > 0.9;   
$adif=($wt > $wq) ? $wq-$sa : $wt - $sa;
if($qsame and $tsame) { $typ="same"; }
elsif($qsame) { $typ="qsubset"; }
elsif($tsame) { $typ="tsubset"; } else { $typ=""; }
if($typ) { print join("\t",$typ,$lq,$wq,$lt,$wt,$sa,$adif,$sm)."\n"; } }'\
> $pt.sametab


# s4: critical, create list of drop-tr ids, decide from sametr.tab which to keep
# fix to pick largest aa,aaqual before largest tr, for same trs
grep same $pt.sametab | perl -ne'($ss,$d,$dw,$t,$tw,$ddt,$mi)=split; $dw{$d}=$dw; $ds{$d}{$t}++; \
END { @d=sort{$dw{$b}<=>$dw{$a}} keys %dw; foreach $d (@d) { $w=$dw{$d}; @s=sort keys %{$ds{$d}}; \
foreach $s (@s) { next if($s eq $d); $sw=$dw{$s}; if($sw<=$w) { $drop=1; $drop{$s}++ unless($keep{$s}); } } \
$keep{$d}++ unless($drop{$d}); } foreach $d (sort keys %drop) { print "drop\t$d\n"; } \
foreach $d (sort keys %keep) { print "keep\t$d\n"; } }' > $pt.samedrop

grep drop $pt.samedrop | cut -f2 | sed 's/$/       /' > $tp.dropids

cat $pt.aa.count | egrep -v '^#|^total' | ggrep -v -F -f $pt.dropids - | sort -k2,2nr | head -1000 | env nam=$pt perl -ne '($aw,$nn)=(split)[1,2]; $aw= $aw-$nn;  $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; $md=$aw[int($n/2)]; print "$ENV{nam}\t  n=$n; aw=$aw; med=$md; sw=$sw; sn=$sn,$an\n"; }'

=cut
