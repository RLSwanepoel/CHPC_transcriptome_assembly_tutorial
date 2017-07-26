#!/usr/bin/perl
# orthomcl_tabulate.pl

=item about

  collection of cmdlines to post-process orthomcl results to summary tables
  and annotated gene family document (ugp.xml) for web search/report
  
=item author

  from EvidentialGene collection
  D. Gilbert, gilbertd at indiana edu  2010..2012
    
=cut


use FindBin;
use lib ("$FindBin::Bin/.."); # assume evigene/scripts/omcl/ << this path
my $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/..";  

use strict;   
use Getopt::Long;
use constant VERSION => '2013.10.08'; # 

## user options
my $DEBUG=(defined $ENV{debug}) ? $ENV{debug} : 1;
my $GTAG= $ENV{gtag} || "ARP";
my $IDPRE= $ENV{idprefix} || "ARP9_G";
my $mcl= $ENV{MCL} || "/bio/bio-grid/mb/mcl9/"; ## findapp('mclcm') and 'mcxdump'
## omcl_mcl2cluster -I $INFLATE is run option .. granularity of grouping (small = narrow groups; 1.5 default input omcl)
my $mcl2INFLATE = 3;
my $MINCOMMON= $ENV{mincommon} || 0; #  || $nspecies - 2;  # $SMIN
  
my $orun="Jan_25";      # get from dirlist
my $phyla="arp11u11";   # from dirlist
my($bpofile, $nspecies, $specieslist, $sppgenes, $omclpath, $logfile);

# my $spl="antc,anth,aphid,apis2ref,bombusimp,bombusterr,daphnia,drosmel,human,trica,wasp"
# my $spl='ACEPH,AECHI,AMELL,CFLOR,HSALT,LHUMI,PBARB,SINVI,bombusimp,bombusterr,wasp'
#^ parse from $orun/parameter.log  SPECIES ..

my $optok= GetOptions( 
  "omclpath|datapath=s", \$omclpath, 
  "IDPREFIX=s", \$IDPRE, "GTAG=s", \$GTAG, # dont need both, make GTAG from IDPREFIX
  "logfile:s", \$logfile,
  "MINCOMMON=i", \$MINCOMMON,  
  "mcl2INFLATE=i", \$mcl2INFLATE,  
##  "ablastab=s", \$aablast,   # option traa-refaa.blastp.tall4 table for asmrna_dupfilter2.pl classifier
##  "CDSBLAST_IDENT=i", \$CDSBLAST_IDENT, "CDSBLAST_EVALUE=i", \$CDSBLAST_EVALUE,  
##  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
##  "tidyup!", \$tidyup, 
##  "dryrun|n!", \$dryrun,
  "debug!", \$DEBUG,
);

die "EvidentialGene orthomcl_tabulate.pl VERSION ",VERSION,"
  convert orthomcl outputs to tables
Usage: ...
  opts: -omclpath=path/to/omclouts/ -MINCOMMON=$MINCOMMON .. -logfile -debug 
" unless($optok); 

MAIN: {
  chdir($omclpath) if($omclpath and -d $omclpath);
  omclout_startup();  # set orun, phyla from dirlist
  ($nspecies, $specieslist, $sppgenes)= omcl_getspecies();
  $MINCOMMON= $nspecies - 2 if($MINCOMMON==0);
  die "#ERROR: nspecies=$nspecies, list=@$specieslist\n" if($nspecies<2);
  
  omcl_gntab();
  omcl_pastebpo();
  omcl_count();
  omcl_gclass();
  omcl_gcommon();

  omcl_avematch(); # eval_orthogroup_genesets.pl, need input blastp.tall4 tables

  omcl_consensusdef();  # needs prepared names/*
  omcl_genegroupdoc();  # needs condef, names, ..
  omcl_groupidtab();

  # omcl_mcl2cluster();  # option or later
        
  omclout_finish(); # gzip big files
}

#-----------------------------------------
sub sysrun
{
  my @CMD= @_;
  if($DEBUG) { warn "#sysrun: ",join(" ",@CMD),"\n"; }
  return system(@CMD);
}

## input gene ids have species_ prefix, omcl adds species: in front; remove dupl prefix  spp:spp_ID

sub omclout_startup
{
  #above# chdir($omclpath) if($omclpath and -d $omclpath);
  opendir(D,"./"); my @fs= readdir(D); close(D);
  ($bpofile)= grep /\.bpo/, @fs;
  $phyla= $bpofile;  $phyla =~ s/.bpo//; $phyla =~ s/_omcl$//; #?
  ($orun) = grep /(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)/i, @fs; # better way? set outdir in omcl run?

  # loggit(0,"input fileset=$bpofile,$orun/all_orthomcl.out");
  die "ERR: missing orthomcl speciesblast.bpo, Date/all_orthomcl.out fileset\n" # loggit(LOG_DIE,...)
    unless( -f $bpofile and -d $orun and -f "$orun/all_orthomcl.out");
  sysrun('perl', '-pi.old', '-e', 's/ \w+:/ /g;', $orun.'/all_orthomcl.out');
}

sub omclout_finish
{
  sysrun("echo rm ${phyla}_omcl_bpo.se ${phyla}_omcl_bpo.idx  $orun/all_orthomcl.out.old");
  sysrun("gzip --fast $orun/all_blast.bbh $orun/all_orthomcl.out $orun/tmp/all_ortho.*");
  sysrun("gzip --fast ${phyla}_omcl.bpo  ${phyla}_omcl.gg  ${phyla}_omclgns2.tab");
}
  
sub omcl_getspecies
{
  my @spp=(); my %sppgenes; my $at;
  open(I,"$orun/parameter.log");
  while(<I>) {
    if(/^#+SPECIES/) { $at=1; }
    elsif($at==1 and /^#/) { $at= 0;} #  last; 
    elsif($at==1 and /^\s+(\w+)\s+(\d+) genes/) { push @spp, $1; $sppgenes{$1}=$2; }
    # if(/^#+PARAMETERS/) { $at="param"; }
    # if(/^#+FILES/) { $at="files"; }
    # if(/^#+START TIME/) { $at="time0"; }
    # if(/^#+END TIME/) { $at="time1"; }
  } close(I);
  my $nspp= @spp;
  warn "#info: nspecies=$nspp, list=@spp\n" if($DEBUG);
  return ($nspp, \@spp, \%sppgenes);
}


sub  omcl_gntab
{
  open(I,"$orun/all_orthomcl.out");
  open(O,">${phyla}_omclgn.tab");
  my( %alts, %didgn);
  while(<I>) { 
    if(/^ORTHOMCL/) { 
    my($om,$gn)=split /:\s+/,$_,2;  my @gn= split" ",$gn; 
    my($og,$ng,$nt)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/;  
    my %didg=(); 
    foreach my $gn (@gn) { 
      my($sp)=split"_",$gn,2; (my $g1=$gn)=~s/[pt\.]\d+$//;  
      $didgn{$g1}++; 
      # $isaltNOT=($didgn{$g1}>1)?1:0; $alts{$sp}++ if($isalt);
      print O "$GTAG$og\t$gn\n"; # unless($isalt); 
      }
    }
  }
  # if(%alts){ print O "#alttr-drops: "; foreach $s (sort keys %alts){ print O "$s=$alts{$s}, "; } print O "\n"; } 
  close(O); close(I); 
}


sub omcl_pastebpo
{
  # sysrun("cat ${phyla}_omclgn.tab ${phyla}*_omcl.bpo | env idprefix=$GTAG $EVIGENES/omcl/pastebpo.pl > ${phyla}_omclgns2.tab");
  # omcl/pastebpo.pl
  my $tag=$GTAG; # $ENV{idprefix} || "ARP";
  my( %gog, $nout,);

  open(I,"${phyla}_omclgn.tab");
  while(<I>){ if(/^($tag\d+)\t(\S+)/) { $gog{$2}=$1; } } close(I);

  # open(I,"cat ${phyla}*_omcl.bpo |");
  open(I,$bpofile) or die "bpofile:$bpofile";
  open(O,">${phyla}_omclgns2.tab");
  while(<I>){
    if(/^\d+;/) { # xxx.bpo
      my @v= split";"; 
      my($g1,$g2,$ev,$pi)= @v[1,3,5,6];
      next if($g1 eq $g2);
      my $og1= $gog{$g1};
      my $og2= $gog{$g2};
      next unless($og1 and $og1 eq $og2);
      print O join("\t",$g1,$g2,$og1,$ev,$pi),"\n";
      $nout++;
      }
  }
  # warn "#pastebpo n=$nout\n" if($DEBUG);
  close(I);

  open(I,"${phyla}_omclgns2.tab");  #? drop this intermediate table? only used here in sub?
  open(O,">${phyla}_omclgn2sum.tab");

  my($g1,$g2,$og,$ev,$pi);
  my($lg,$log,$ng,$sev,$spi,$mpi,$mg,$mev)= (0) x 10;
  sub dumpg{ if($ng>0){
    my $aev=sprintf"%.4g",$sev/$ng; my $api=int($spi/$ng); 
    print O join("\t",$lg,$log,$ng,$aev,$api,$mg,$mev,$mpi),"\n"; } 
    $mpi=$mg=$mev=$lg=$ng=$sev=$spi=0; 
    } 
    
  while(<I>) {
    ($g1,$g2,$og,$ev,$pi)=split; dumpg() unless($g1 eq $lg); 
    $lg=$g1; $log=$og; $ng++; $sev+= $ev; $spi += $pi; 
    if($pi>$mpi){ $mpi=$pi; $mg=$g2; $mev=$ev; } 
    }
  dumpg();
  close(O); close(I);
  warn "#info: pastebpo=$nout\n" if($DEBUG);
}


sub omcl_count
{  
  open(I,"$orun/all_orthomcl.out");
  open(O,">${phyla}-orthomcl-count.tab");

  my($sog,$snt,$sng)=(0) x 3;
  my @spl= @$specieslist;  # $spl=$ENV{spl}; @spl=split",",$spl;
  my %drop; @drop{@spl}= (0) x $nspecies; 
  print O join("\t","OID","Nt","Ng",@spl),"\n";
  while(<I>) {  
    if(/^ORTHOMCL/) { 
    my ($om,$gn)=split /:\s+/,$_,2; my @gn= split" ",$gn; 
    my %spc; @spc{@spl}= (0) x $nspecies; 
    my %didgn=(); 
    foreach my $sg (@gn) { 
      (my $g1=$sg)=~s/[pt\.]\d+$//; # $isaltNOT=($didgn{$g1}++ > 0)?1:0; 
      my($sp,$gn)=split"_",$sg,2; $spc{$sp}++; ## unless($isalt); 
    }
    my($og,$g,$t)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/; 
    my @sc= @spc{@spl}; print O join("\t",$og,$t,$g,@sc),"\n"; 
    $sog++; $sng+=$g; $snt+=$t;
    } 
  }
  close(O); close(I);
  
  $snt= int(100*$snt/$sog)/100; $sng=int(100*$sng/$sog)/100;
  warn "#info: ${phyla}-orthomcl-count.tab  ngr=$sog, taxa/gr=$snt, spp/gr=$sng;\n" if($DEBUG);
}

=item omcl_gclass

# summary gene class table
#? add spp correlations in genes/group ? vxy, vxx, vyy table == distance mat
#? add OrGrp count: number ortho groups species belongs; Orth1 + OrDup-grps
# iskip=8 == human  DANG NO, index from col=0 OID, human skip=11

=cut

sub omcl_gclass
{  
  open(I,"${phyla}-orthomcl-count.tab");
  open(O,">${phyla}-orthomcl-gclass.tab");

  ## FIXed: use %sppgenes to fill in Uniq1, inGene
  #cat ${phyla}-orthomcl-count.tab | env iskip=-1 perl -ne
  my $iskip= -1; # $ENV{iskip};
  my $miss= 1; # $ENV{miss};
  my $onlys= 0; # $ENV{onlys}; # onlyspecies

  my(@hd, @spp, $MIS1);
  my @osp=(); if($onlys) { @osp= split /[,\s]+/,$onlys; }

  print O join("\t",
    qw(species inGene oGene Uniq1 UDup Orth1 OrDup nGroup UniqGrp OrGrp OrMis1 Guniq Gmax Gmin)),"\n";

  my(%sc,%sg,%sgu,%sgh1,%sghp,%sgOrG,%sgMIS,%suni,%smax,%smin);
  while(<I>) {
    my @v=split; my $nv=$#v; 
    if(/^OID/){ @hd=@v; @spp=@v[3..$nv]; $MIS1= @spp - $miss; next; } 
    elsif(/^\D/){ next; } #? ok
    
    my $nt=0; 
    #count instead# my $ntbad=$v[1];  my $ngbad=$v[2]; # fixme: ntbad : count
    # my @s=map{ $hd[$_]."=".$v[$_]; } (3..$nv); # fixme
    my %csp=(); my @s= map{ my $c=$v[$_]; $nt++ if($c>0); $csp{ $hd[$_] }=$c; $hd[$_]."=".$c; } (3..$nv);
    if(@osp) { my $ok=0; map{ $ok++ if( $csp{$_}) } @osp; next unless($ok); }

    my($vm,$smax,$vm2,$smin)= (0) x 10; my $vn=99999; 
    map{ my($s,$v)=split"="; 
      $sc{$s}++ if($v>0); $sg{$s} += $v; 
      $sgu{$s}  += $v if($nt==1); 
      $sgh1{$s} += $v if($nt>1 and $v==1); 
      $sghp{$s} += $v if($nt>1 and $v>1); 
      $sgOrG{$s}++ if($nt>1 and $v>=1); 
      $sgMIS{$s}++ if($v==0 and ($nt>=$MIS1 or ($iskip>=0 and $nt>=$MIS1-1 and $v[$iskip]==0))); 
      if($v<=$vn){ $smin=($v==$vn)?"$smin,$s":$s; $vn=$v;} 
      if($v>=$vm){ $smax=($v==$vm)?"$smax,$s":$s; $vm2=$vm; $vm=$v;}  
    } @s; 
    $suni{$smax}++ if($nt==1); 
    $smax{$smax}++ unless($nt < 2 or $smax =~ /,/);  
    $smin{$smin}++ unless($nt < 4 or $smin =~ /,/);  
  }
  foreach my $s (sort keys %sc) { 
    my($ing,$ong,$u,$h,$d,$c,$ngu,$ngo1,$ngo2,$ngo0,$nog,$ug1)= (0) x 20;
    $ing= $sppgenes->{$s} || 0; 
    $u=$suni{$s}||0; $h=$smax{$s}||0; $d=$smin{$s}||0; $c=$sc{$s}||0; $ong=$sg{$s}; 
    $ngu=$sgu{$s}||0; $ngo1= $sgh1{$s}||0; $ngo2= $sghp{$s}||0; $ngo0=$sgMIS{$s}||0; $nog= $sgOrG{$s}||0;  
    $ug1= ($ing>=$ong) ? $ing - $ong : "na";
    print O join("\t",$s, $ing, $ong, $ug1, $ngu,$ngo1,$ngo2,  $c,$c-$nog,$nog,$ngo0,  $u,$h,$d),"\n"; 
  }      
  close(O); close(I);    
  warn "#info: ${phyla}-orthomcl-gclass.tab;\n" if($DEBUG);
}
 
# perl -pi -e's/bombusimp/bombi/; s/bombusterr/bombt/;s/trica/tribol/; s/apis2ref/apis/; ' ${phyla}-orthomcl-gclass.tab

=item omcl_gclass reformat table this way

Fish orthology gene groups summary (OrthoMCL).
            ---------- GENES -----------    -------- GROUPS --------------  
            nGene Orth1 OrDup Uniq1 UDup    nGroup  OrGrp OrMis1 UniqGrp
            ----------------------------    -----------------------------  
killifish   30000 13958  9397* 2837 3808    17448   16587    107    861
zebrafish   29004 10053 12973  3730 2248    15081   14513    284    568
tilapia     21442 12656  7463   858  465    15207   15086    201    121
stickleback 20875 12669  6327  1560  319    14866   14797    264     69
medaka      19732 12130  5456  1660  486    14145   14015    538    130
tetraodon   19646 11624  5544  2263  215    13725   13649    523     76
-----------------------------------------------------------------------
  Uniq1,UDup  = single-copy and duplicated species-unique genes
  Orth1,OrDup = single-copy and duplicated orthologous genes
  UniqGrp,OrGrp = species-unique and orthologous groups
  OrMis1  = groups missing in species that all other species have

=cut

sub omcl_gcommon
{
  open(I,"${phyla}-orthomcl-count.tab");
  open(O,">${phyla}-orthomcl-gcommon.tab");
  # cat ${phyla}-orthomcl-count.tab | env min=6 perl -ne 
  # my $SMIN= $ENV{mincommon} || $nspecies - 2;
  my(%have, @hd, $ng);
  while(<I>) {
    if(/^OID/){ @hd=split;  next; }  # FIXME:  s/bombusimp/bombi/; s/bombusterr/bombt/;s/trica/tribol/; s/apis2ref/apis/; 
    elsif(/^\D/){ next } 
    my @v=split; my $nv=$#v; my $nt=$v[1]; 
    if($nt>=$MINCOMMON) {  $ng++; 
      for my $i (3..$nv) { my $c=$v[$i]; my $s=$hd[$i]; my $have=($c>0)?"have":"miss"; $have{$s}{$have}++; } 
      } 
  }
  print O "Common gene families presense, ncommon=$ng, min taxa=$MINCOMMON\n"; 
  print O join("\t",qw(Species Have Miss))."\n";
  foreach my $s (sort{ $have{$a}{miss} <=> $have{$b}{miss} or $a cmp $b } keys %have) {  
    my($h,$m)= @{$have{$s}}{qw(have miss)}; 
    print O join("\t",$s,$h,$m)."\n"; 
  }
  close(O); close(I);
  warn "#info: ${phyla}-orthomcl-gcommon.tab;\n" if($DEBUG);
}

=item omcl_avematch

  set karp=$kfish/prot/kfish1ball/omclkf2
  set spp=human
  
  # $EVIGENES/makeblastscore2.pl fish8-xxx.blastp.gz > fish8-all.tall4
  env aa1=fish10main.aa.qual keepho="human:" tall=1 $EVIGENES/makeblastscore2.pl \
    fish10main-fish10main.blastp.gz > fish10main_human.tall4 
  
  .. split fish8-all.tall4 to fish8-{eachspecies}.tall4
  .. OR split blastp by '^species:', run makeblast tall4 on each
  
  $EVIGENES/eval_orthogroup_genesets.pl -nodebug -bitmed -mintaxa 3 \
  -out fish8-$spp.topout2 -tallscore fish8-$spp.tall4 \
  -groupgene $karp/fish8_omclgn.tab -groupcount $karp/fish8-orthomcl-count.tab\
  -geneaa $karp/fish8_omcl.aa.gz 
  
  #........... omclkf2/Aug_25/  : removed alts
  -- slight changes only from alt removal
  -- adding best of kfish alts/main doesnt change these scores; check alt for missed omcl groups
  cat fish8-{kill,medaka,stick,tetra,tilapia,zebr}*.topout2  | $EVIGENES/eval_orthogroup_genesets.pl -intab stdin -digit 0
  
   ...  Common Fish groups ... 
  Geneset  : nGene    Bits    dSize   rGene   rBits   outlierSize
  killifish: 10172    718     9       10172   718     x2Tiny=244 (2.3%)   x2Big=92 (0.9%)
  medaka   : 10172    699     -29     10172   699     x2Tiny=600 (5.8%)   x2Big=20 (0.1%)
  sticklebk: 10172    721     -19     10172   721     x2Tiny=387 (3.8%)   x2Big=22 (0.2%)
  tetraodon: 10172    688     -43     10172   688     x2Tiny=576 (5.6%)   x2Big=19 (0.1%)
  tilapia  : 10172    753     19      10172   753     x2Tiny=60 (0.5%)    x2Big=67 (0.6%)
  zebrafish: 10172    685     11      10172   685     x2Tiny=195 (1.9%)   x2Big=105 (1%)
   ...  All groups ...
  Geneset  : nGene    Bits    dSize   rGene   rBits   outlierSize
  killifish: 16319    574     8       15547   602     x2Tiny=576 (3.7%)   x2Big=311 (2%)
  medaka   : 16319    525     -34     13909   616     x2Tiny=1058 (7.6%)  x2Big=44 (0.3%)
  sticklebk: 16319    560     -26     14624   625     x2Tiny=832 (5.6%)   x2Big=53 (0.3%)
  tetraodon: 16319    515     -42     13580   619     x2Tiny=867 (6.3%)   x2Big=56 (0.4%)
  tilapia  : 16319    597     18      14896   654     x2Tiny=146 (0.9%)   x2Big=167 (1.1%)
  zebrafish: 16319    524     18      14634   584     x2Tiny=311 (2.1%)   x2Big=318 (2.1%)
  human    : 16319    441     45      12752   564     x2Tiny=68 (0.5%)    x2Big=280 (2.1%)
  xenopus  : 16319    429     0       12595   556     x2Tiny=261 (2%)     x2Big=131 (1%)

=item avematch reformat table this way

  Fish species average match to orthology gene groups.
           Common groups(n=10172)        All groups 
  Geneset    cBits  dSize     tBits  rGene   Small outliers
  --------- --------------------------------------
  killifish  718     9         574   15547   3.7%
  medaka     699     -29       525   13909   7.6%
  sticklebk  721     -19       560   14624   5.6% 
  tetraodon  688     -43       515   13580   6.3% 
  tilapia    753     19        597   14896   0.9% 
  zebrafish  685     11        524   14634   2.1% 
  ------------------------------------------------
    cBits = bitscore average to common groups
    tBits = bitscore average to all groups
    dSize = average size difference from group median
    Small outliers : percent species genes < 2sd of median gene size in group 
    rGene = number of gene groups found in species
   
=cut

sub omcl_avematch
{

}


=item omcl_consensusdef

  cat names/*.names ${phyla}_omclgn.tab | \
  env gtag=$GTAG idprefix=$IDPRE $EVIGENES/omcl/arp_condef2.pl \ 
  > ${phyla}_omclgn.consensus_def.txt 
  
  cat ${phyla}_omclgn.consensus_def.txt | env recase=1 count=withID debug=1 \
  $EVIGENES/bestgenes_puban_wasp.pl | sed 's/^$IDPRE//' | sort -k1,1n | sed 's/^/$IDPRE/; s/TE:TE:/TE:/' 
    > ${phyla}_omclgn.consensus_def.rename.txt
  mv ${phyla}_omclgn.consensus_def.txt  ${phyla}_omclgn.consensus_def.txt0
  mv ${phyla}_omclgn.consensus_def.rename.txt ${phyla}_omclgn.consensus_def.txt 

=item pull names

  orthodb6.aa:  gzgrep '^>' $pt.aa.gz | perl -ne\
  '($id,$pid,$uid,$na)=m/>(\S+)\s(\S+)\s+(\w*)\s*(.*)$/; $na=~s/IPR.*$//; print "$id\t$na\n" if($na =~ /[a-z]/);' \
  > ../names/$pt.names
  
  evigene.aa: gzgrep '^>' $pt.aa.gz | perl -ne\
  '($id)=m/^>(\S+)/; ($na)=m/Name=([^;\n]+)/; print "$id\t$na\n" if($na =~ /\w/);' \
  > ../names/$pt.names
  
  refseq.aa:
     
=cut

sub omcl_consensusdef
{

}

=item omcl_genegroupdoc

  cat \
  ${phyla}_omclgn2sum.tab  \
  ${phyla}_omclgn.consensus_def.txt \
  names/*.names  \
  $orun/all_orthomcl.out | \
  env xml=1 date=20120125 clade=HymenoInsectArp title='Arthropod gene group' gtag=$GTAG idprefix=$IDPRE \
  $EVIGENES/omcl/genegroupbpo.pl  
  > ${phyla}_genes.ugp.xml

=item genebriefdoc

	set spt=Nasvi
	set spt=Acyrthosiphon
  cat ${phyla}_genes.ugp.txt | egrep "GeneID|ntaxa|ngene|occur|descript|$spt" | env spt=$spt perl -pe \
  's/^ +//; if(m/similarity:/ and /$spt/){ ($id)=m/iden: (\d+)/; ($d)=m/acc: (\w[^;\s]+)/; $_="$d/$id,"; s/^/$spt: / if $s; $s=0;} \
  else{ $s=1; print "\n" if s/\s*GeneID:\s//; s/\n/ /;} BEGIN{$spt=$ENV{spt};} END{print"\n";} ' \
  >  ${phyla}_genes.ugp_brief.txt

=cut

sub omcl_genegroupdoc
{

}


=item omcl_groupidtab 

  # do this not for 11 taxa (all) but drop human, daphnia?
  cat ${phyla}-orthomcl-count.tab | env nt=11 gtag=$GTAG perl -ne
  '($od,$nt,$ng,@c)=split; print "$ENV{gtag}$od\t\n" if($nt==$ENV{nt});' 
	  > ${phyla}_omcl11.allgr.oids


  ggrep -F -f ${phyla}_omcl11.allgr.oids ${phyla}_omclgn.tab | cut -f2 | sed 's/^/gid /' > ${phyla}_omcl11.allgr.gids
  # use all.gids for restricted distance matrix restricted to common genes.

=cut

sub omcl_groupidtab
{

}


#............. 
# mcl gene group tree info

sub omcl_mcl2cluster
{

  chdir($orun);
  sysrun("$mcl/bin/mclcm", 'tmp/all_ortho.mtx', '-a', '-I $mcl2INFLATE --shadow=vl -te 2'); # >& log.mcm1 # fork
  sysrun("$mcl/bin/mcxdump", '-imx-tree mcl.cone', '--newick', '-o', "${phyla}.newick", '-tab', 'tmp/all_ortho.idx');
  chdir("../");
  
=item mcl newick tree

  # cat ${phyla}_omclgn.tab ${phyla}.newick | perl -ne .. > ${phyla}.newick4
  set newk=${phyla}.newicki2

  ## improve these perls, mcltreesplit.pl
  cat ${phyla}_omclgn.tab  ${phyla}_omclgn.consensus_def.txt  ${newk} |
  env gtag=$GTAG perl -ne 'BEGIN{$gtag=$ENV{gtag};} 
  $p=1; if(/^$gtag\d\D+(\d+)/){ $g=$gtag.$1; s/\(LOC\d+\)//;  s/src=\S+//; 
  $gd{$g}= (m/\s(\S.+)$/)? $1.";" : "";  $p=0; }  
  elsif(/^($gtag\d+)\s+(\S+)$/){$ag{$2}= $1; $p=0;} 
  elsif(m/\(([^\)]+)\)/) { $gs=$1; @g=split",",$gs; %ga=(); $de="";  
  map{ $ag=$ag{$_}; $ga{$ag}++ if($ag); } @g; 
  map{ $de .= $gd{$_} } sort keys %ga; $ga= join ",", sort keys %ga;
  s/\([^\)]+\)/\($ga\)/; s/$/ # $de/ if $de; } print if $p; '
   > ${newk}b
  
  cat ${newk}b | env gtag=$GTAG perl -ne
  '$s=$_; s/\s*\#.*$//; s/[\,\s]+$//; $s=~s/$GTAG/$IDPRE/g; print $s unless($_ eq $ll); $ll=$_;' 
   > ${newk}c
   
  cat ${newk}c | env tag=$GTAG  perl ../omclw/mcltreesplit.pl > ${phyla}.clusters.tab

=cut

=item add clusters to ugp.xml

  cat  ${phyla}.clusters.tab  ${phyla}_genes.ugp.xml |
  perl -ne'if(/^ARC1_/){ chomp; ($c,$gc)=split"\t"; while( $gc =~ m/($GTAG\w+)/g) { $gc{$1}= $gc; } } 
  elsif(/^ARDE_/){next;} else { if(m/<GeneSummary/) { $gc= (m/id=.\w+:(\w+)/) ? $gc{$1} : ""; } 
  elsif( $gc and m,</GeneSummary>,){ print "<related_gene_groups>\n$gc\n</related_gene_groups>\n"; }  
  print; }'
  >  ${phyla}_genes.ugp.xml2

=cut

}

=item add

genegroupvenn/
  waspinsect_genefam_venn.pdf
  venn diagram of genes for 5 species, shared and unique (wasp,bee,ant,beetle,aphid)
  from http://bioinformatics.psb.ugent.be/webtools/Venn/
 
overgroups/   
  table.overgroups.nasonia.txt = list of over/under abundant groups in nasonia with statistics
  table.overgroups.* = same for other species

=cut

__END__

