#!/usr/bin/perl
# makeblastscore.pl   predictor-ortho.blastp predictor-self.blastp > predictor.blscore
#  cut from annotate_predictions.pl

## FIXME: need to look for protein part scores to combine as one
# melon2.% gzgrep 'EG2ap020688t1' blastp/aphid2pub8d-uparphumbac.blastp.gz | grep E0VM80_PEDHC
# # Query:  EG2ap020688t1 aalen=2955,92%,complete
# EG2ap020688t1   E2BNB5_9HYME    54.16   1852    610     43      1258    2955    1070    2836    0.0     1780 << top score is only 1/2 prot
# EG2ap020688t1   E2BNB5_9HYME    47.74   1175    475     32      1       1126    1       1085    0.0      908 << rest
# ..
# EG2ap020688t1   E0VM80_PEDHC    53.72   1763    615     44      1262    2955    1108    2738    0.0     1692 << best score when summed
# EG2ap020688t1   E0VM80_PEDHC    50.25   1184    446     32      1       1122    1       1103    0.0     1027 <<
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

#** fixme, makeblastscore needs correction for small-overlap bitscore parts; allow overlap < 1% span?
# .. problem esp those long, repetitive prots: dyenin hvy chn, etc.
#gzgrep Nasvi2EG013019t1 aaeval/bp1-apis2-nvit2_evg11d_all.aa.blastp.gz  | head
#Nasvi2EG013019t1        apis2:XP_396548.4   66.41   3168    1033 14      27      3175    12      3167    0.0     4380
#... slight overlap in query (not targ): 10aa / (3175,800)
#Nasvi2EG013019t1        apis2:XP_396548.4   64.60   644     224  1       3164    3803    3269    3912    0.0      905
# bits=5285  == same full bitscore as below equiv, nearly identical model
#
#gzgrep XP_003426444.1 aaeval/bp1-apis2-nvit2_ncbiref.aa.blastp.gz | head
#ncbiref2:XP_003426444.1 apis2:XP_396548.4   66.46   3166    1036 14      24      3175    12      3165    0.0     4380
#ncbiref2:XP_003426444.1 apis2:XP_396548.4   64.30   647     227  1       3176    3818    3266    3912    0.0      905


use strict;

my $KeyHOMOLOG="homolog";  # make config OPTION
my $KeyPARALOG="paralog";
my $KeyINSPLIT="insplit";
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

my $skipho= $ENV{skipho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $keepho= $ENV{keepho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $onlyquery= $ENV{onlyquery} || ""; # species query to keep, or all
$skipho =~ s/[,\s]+/\|/g;
$keepho =~ s/[,\s]+/\|/g;
$onlyquery =~ s/[,\s]+/\|/g;

## pick 1 best target group: TBEST, given ID prefix patt TPRE eg:  '^....'
my $TPRE= $ENV{tpre}||"";
my $TBEST= $ENV{tbest}||""; # $TPRE="" unless($TBEST);
my $pMINLOW= $ENV{pmin} || 0.3; # was 0.3?
my $TALL= $ENV{tall}||0; # all targets >= pmin;
my $pMINBEST2 = 0.75; # dont use TBEST if below 0.5*topscore
my $ADDALIGN= $ENV{align} || $TALL;
my $geneaa= $ENV{aa} ||"";

my( $bother, $bself)= @ARGV;
( $bother, $bself) = ( $bself, $bother) if($bother =~ /self/);

## fixme for alttr : idtag = t[2-n] : drop from paralog=
my $altids= $ENV{altids} ||""; # table of altid  mainid
my $ALTKEY= $ENV{alt} || 't'; ##'t\d+$';
$ALTKEY .= '\d+$' if($ALTKEY =~ /\w/  and $ALTKEY !~ /\W/);

die "usage: makeblastscore genes-ortho.blastp genes-self.blastp > genes.blscore \n"
  unless( -f $bother or -f $bself);
  
my (%bself, %bparalog, %bother, %tother, %balt, %blen, $cat, $bsort, $cmd, %bspans, $lq, $bmax);
      
sub bint { local $_=shift; return (/e\+/) ? int($_) : $_; }

my $haveqlen=0;
if($geneaa) { # warn:  not -f $geneaa
  if($geneaa =~ /count/) { open(AASIZE,$geneaa) or die "FAIL: faCount  aa=$geneaa ..."; }
  else { open(AASIZE,"faCount $geneaa |") or die "FAIL: faCount  aa=$geneaa ..."; }
  while(<AASIZE>) { my($id,$al)=split; $blen{$id}=$al; } close(AASIZE); $haveqlen=1;
}


my %altids=(); my ($faltid,$naltid)=(0,0);
if($altids) {
  open(F, $altids) or die "ERR: altids=$altids bad file of ids\n";
  while(<F>) { if(/^\w/) { my($aid,$mid)=split; $altids{$aid}=$mid; $naltid++; } } close(F);
}

if( $bself and -f $bself ) {
# $cat= ($bself =~ /\.gz/) ? "gunzip -c" : "cat";
# $cmd= join(" ",$cat, $bself, "|");
$cmd= ($bself =~ /\.gz/) ? "gunzip -c $bself |" : $bself;
open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
while(<GSCORE>) { 
 next unless(/^\w/); my @v=split; 
 my($q,$t,$bits,$aln,$mis)= @v[0,1,-1,3,4]; $bits=bint($bits);
 my $qg=$q; my $tg=$t; 
 next if($onlyquery and $q !~ m/$onlyquery/);
 if($altids) { $qg=$altids{$q}||$q; $tg=$altids{$t}||$t; }
 elsif($ALTKEY) { $qg=~s/$ALTKEY//; $tg=~s/$ALTKEY//; }
 
 if($q eq $t) { $bself{$q}= $bits unless($bself{$q}); 
   $blen{$q} += $aln unless($haveqlen); # always == length for self? but for long aa broken blast
 } 
 elsif($qg eq $tg) { $balt{$q}{$t}=1; $faltid++; } # alttr, what? need to keep ids for other match
 else { $bparalog{$q}="$bits,$t" unless($bparalog{$q}); }
 
 } close(GSCORE);
}

if( $bother and -f $bother ) {
# $cat= ($bother =~ /\.gz/) ? "gunzip -c" : "cat";
# $cmd= join(" ",$cat, $bother, "|");
$cmd= ($bother =~ /\.gz/) ? "gunzip -c $bother |" : $bother;
open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
while(<GSCORE>) { 
  unless(/^\w/) { 
   #?? if(/^# Query: (\S+)/) { my $q=$1; my($al)=m/len=(\d+)/; $blen{$q}=$al if($al); $haveqlen++ if($al); }
   next;} 
  my @v=split; 
  # my($q,$t,$bits)= @v[0,1,-1];  #old
  my($q,$t,$bits,$aln,$mis,@bspan)= @v[0,1,-1,3,4, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
  # bspan = ($qb,$qe,$sb,$se);
  next if($onlyquery and $q !~ m/$onlyquery/);

  if($lq and $q ne $lq) {
    ##my($lbits, $lt)= bestscore($lq); 
    my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
    my($tbits, $tt)= split",", $bother{$lq};
    $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
	    or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score
    %bspans=(); $bmax=0; $lq="";
  }

  next if($skipho and $t =~ m/$skipho/);
  next if($keepho and $t !~ m/$keepho/);

  $bits= bint($bits);
  if($q ne $t) { 
    my $aident= _max(0,$aln-$mis); # other way to calc: $aident = $pctident * $aln;
    sumscore( $q, $t, $bits,$aln,$aident, @bspan) if($bits > $pMINLOW * $bmax or $bspans{$t}); 
      # ?? limit # lowscore targets here? careful, 1/2 score * 2 can be best
    $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
    $bmax= $bits if($bits > $bmax);
    }
  elsif($q eq $t and not $bself) {
    $bself{$q}= $bits unless($bself{$q}); 
    $blen{$q} += $aln unless($haveqlen); # always == length for self? but for long aa broken blast
    if($TALL) {
      my $aident= _max(0,$aln-$mis); # other way to calc: $aident = $pctident * $aln;
      sumscore( $q, $t, $bits,$aln,$aident, @bspan); 
      # $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
    }
  } 

  $lq= $q; 
 } close(GSCORE);

  my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
  my($tbits, $tt)= split",", $bother{$lq};
  $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
    or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score

}      

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
sub sumscore {
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

sub bestscore {
  my($q)= @_;
  my($maxb,$maxt,$maxa,$maxi)= (0) x 5;
## 2011.11: best tprefix (pick best of prefixa, but prefixb if no prefixa)
  my(%maxb,%maxt);
  foreach my $t (sort keys %bspans) {
    my @sp= @{$bspans{$t}}; my ($tbit,$ta,$ti)= (0) x 9;
    foreach my $sp (@sp) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident)= @$sp;
      $tbit += $xbit; $ta+= $aln; $ti+= $aident;
    }
   if($TALL) { $tother{$q}{$t}="$tbit,$ti,$ta"; } ##  $bother{$q}="$bits,$t,$aident,$aln"
   if($tbit > $maxb) { $maxb=$tbit; $maxt= $t; $maxa=$ta; $maxi=$ti;}
   if($TPRE) { my ($tp)= $t =~ m/^($TPRE)/;  $tp||="other";
     if($tbit>$maxb{$tp}) { $maxb{$tp}=$tbit; $maxt{$tp}=$t; } # add ta,ti
     }
  }

 if($TPRE) {
    my ($tp)= $maxt =~ m/^($TPRE)/;
    unless($tp =~ m/$TBEST/) { 
      my $maxb1=$maxb{$TBEST}; my $maxt1= $maxt{$TBEST}; 
      if($maxt1) { $maxt=$maxt1; $maxb=$maxb1; }# add ta,ti
    }
 }
 return($maxb, $maxt, $maxa, $maxi);
}

# OUTPUT.....................................
my %ball= map{ $_ => 1 } keys %bother, keys %bself;

my ($ng,$nho)=(0,0);
if($TALL) { # just ho table for all q x t
  print join("\t",qw(Query Source Bits Ident Align Qlen))."\n";
  foreach my $gid (sort keys %tother) {
    foreach my $tg (sort{ $tother{$gid}{$b} <=> $tother{$gid}{$a} } keys %{$tother{$gid}} ){
      my $bia= $tother{$gid}{$tg}; my($bit,$idn,$aln)=split",",$bia; 
      my $qlen= $blen{$gid}||0; # get elsewhere? bself? OR include gid x gid scores in tother for self score?
      print join("\t",$gid,$tg,$bit,$idn,$aln,$qlen)."\n";
    }
  }
  exit; # done
}

#?? report balt faltids ?  $balt{qid}{$tid}
foreach my $gid (sort keys %ball) { 
  $ng++;
  my @out=();
  my $s=$bself{$gid};    my($saln)= ($s=~s/,(\d+,\d+)$//)?$1:"";
  my $b=$bother{$gid};   my($baln)= ($b=~s/,(\d+,\d+)$//)?$1:"";
  my $p=$bparalog{$gid}; my($paln)= ($p=~s/,(\d+,\d+)$//)?$1:"";
  if($s or $b or $p) {
    #fixabove: map{ if(/e\+/){ my($b,$t)=split","; $b= int($b); $_="$b,$t"; } } ($s,$b,$p); # 1.033e+04
    ## FIXME: these can be e-notation:  1.23e+4 > convert to digits
    ## 12jan: val=$maxbits,tag,$maxa,$maxi  now; cut maxa,maxi unless requested
    
    # overbestgenes expects format  "ho=score/maxscore,.."
    if($b) { $b =~ s=,=/$s,= if($s); push @out, "$KeyHOMOLOG=$b";  $nho++;}
    else { push @out,"na"; }
    if($ADDALIGN){ $baln=~s=,=/=; push @out, (($baln) ? "align=$baln" : "na"); }
      
    if($p) { $p =~ s=,=/$s,= if($s); push @out, "$KeyPARALOG=$p"; }
    else { push @out, "na"; }
    if($s and ($b or $p)) {  
      $b =~ s,/.*,,; $p =~ s,/.*,,; $s =~ s,/.*,,;
      my($bb,$bl)=($p>$b)? ($p,"pa") : ($b,"ho"); $bb=int(100*$bb/$s); 
      push @out, "pHOBEST=$bb\%$bl"; 
      } else { push @out, "na"; }
  }
    
  print join("\t",$gid,"na",@out),"\n";
}
