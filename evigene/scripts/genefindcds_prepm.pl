#!/usr/bin/perl
# genefindcds.pl

=item about

  genefindcds.pl -genes trassembly_genes.gff -dnasequence genome.fasta [ -introns intron_good.gff ]
    > trassmbly_corrected_with_proteins.gff

  finds best orf, adds proteins and CDS to transcript assembly GFF
  when valid introns.gff are provided, used to check/correct for retained introns,
    and reversed-intron errors from assembly join of 2 genes.
    
  ** ASSUMES gene model is mRNA > exon (optionally CDS)
  ** ASSUMES input gff is ordered by gene records (mRNA/exon,CDS all together per ID)

=item author
  
  don gilbert, gilbertd near indiana edu, 2011
  drawn largely from Brian Haas's PASA scripts
  part of EvidentialGene, evigene/scripts/

=cut

##.. 12.07 update to share package

use constant VERSION  => '20120805'; # revised again cdnain vs gff stranding; 20120731 had bad best-strand
 # '20120731';  # utrorf, cdnain fixups
 # '20120706'; # '20120221'; 

use FindBin; 
use lib ("$FindBin::Bin", "$FindBin::Bin/../lib/"); # find Bio::DB::Fasta in evigene/lib/Bio/... from evigene/scripts/this.pl

use strict;
use warnings;
use Getopt::Long;
# use cdna_proteins;
# use Bio::DB::Fasta; # now as get_dna() require  Bio::DB::Fasta, so can use this w/o bioperl  

use constant { kINTRON2SPLICE_OVER=>1, kINTRON2SPLICE_QUERY=>2, 
               kINTRONERROR_OVER=> -1, kINTRONERROR_INSIDE => -3 }; # only want last 2?
use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon

my $debug= 0;
# my $MINID_CDS= 0; # not used
my $MINAA= 30;  # used 
my $MINEXON= 60; # for cut
my $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
# FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
my $ORF_FULLvPART = 0.85;
my $USE_CDSEXONS = 0;
my $USEGOODLEN=1;
my $DO_INFIX= 0; # which default?
my $NODIFCAN=0;
my $REANNOTATE=0;
my $allowalts= 0; # NOT 1 default? : this kills most changes, good ones..
my $pMAXSCORE = 0.05; # 1/20; 1/50?; dont consider intron cuts for this inscore < pMAX * maxvalidscore
my $AA_cdna_GT_genome= 1.75;  # for prot(cdna) > prot(genome) test; option: 1.75 probably too high default; 1.25 better?

my $BINSIZE   = 5000 ; 
my ($overlaps,$passtypes,$dnasequence,$cdnaseq,@input,$intron2splice,
    $itype,$action,$actid,$typeover,$ok,$mark,$nin);
my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $introntypes='intron'; # option?

my $overlaplist= {};
my $fasta_db= undef;
my $stranded=1; # only this opt?
my $samecds= 0; # KEEPSAMECDS; prefer keep same CDS exons but can extend/shorten protein bounds
my $debugin= undef;
my %cdnaseq; my %cdnahead;
my $nostopcodon=0;

# add -output option to file updated genes
# add -cdnaseq option to pick best prot from asmrna transcript, compare to gff-cds
# FIXME: utrorf: works, but now only for same strand as bestorf; need to keep fwd/rev orfs, and retest against final best.
# ... add option here? to split transcript/gff for strong utrorf cases. also deal with 3+orfs/transcript: recursive?



my $optok= GetOptions(
  "introns=s", \$overlaps, 
  "genes=s", \@input,  
  "dnasequence=s", \$dnasequence,  
  "cdnaseq=s", \$cdnaseq,  
  "exontypes=s", \$exontypes, 
  "action=s", \$action, 
  "allowalts!", \$allowalts, 
  "passtypes=s", \$passtypes,  #??
  "minaa|minprot=i", \$MINAA, ## MINID_CDS not used ; replace with minaa|minprot
  "fullorf:s", \$ORF_FULLvPART,  
  "ratiocdnabest=s", \$AA_cdna_GT_genome, 
  "samecds!", \$samecds, 
  "nostopcodon", \$nostopcodon, 
  "goodlen!", \$USEGOODLEN,  # add MINGOOD
  "goodmin=s", \$MINGOOD,  
  "nodifcancel!", \$NODIFCAN, 
  "fixintronerrors!", \$DO_INFIX, 
  "reannotate!", \$REANNOTATE, 
  "CDSEXONS!", \$USE_CDSEXONS, 
  "debug:i", \$debugin, 
  );

die "usage:
  genefindcds.pl -genes trassembly_genes.gff -dnasequence genome.fasta [ -introns intron_good.gff -cdna cdna.fa ]
    > trassmbly_corrected_with_proteins.gff
" unless($optok and ((@input and $dnasequence and -f $dnasequence ) or $cdnaseq));

if(defined $debugin) { $debug=($debugin>0)?$debugin:1; }

## FIXME? add option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
## this is problem:
## protein=MLSHQLLEDSTMMQMKHGLRQGRENICQGSRLLLIGNVLVDNXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
$MINGOOD= $MINGOOD/100.0 if($MINGOOD > 1); 

##was $ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART > 1); 
 # ^^ -full=0 means only full orfs?; -full=1 means never?? change to -full=1 means full only; -full=0 longest part always
$ORF_FULLvPART=1 unless($ORF_FULLvPART =~ /\d/);
$ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART >= 1); 
$ORF_FULLvPART=1 if($ORF_FULLvPART == 0); # means always longest

$intron2splice= ($overlaps) ? kINTRONERROR_OVER : 0; # only choice?

my $nintron= 0;
if($overlaps) {
  my $ovh; 
     if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
  elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
  else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
  die "bad -overlaps=$overlaps" unless($ok);
  
  # $overlaplist= 
  $nintron= collect_overlaps($ovh); close($ovh);
}
my $hasintrons= ($nintron>0 and scalar(%$overlaplist))?1:0;

if($cdnaseq) {
  my $ovh; 
  if($cdnaseq =~ /.gz$/) { $ok= open(OVR,"gunzip -c $cdnaseq |");  $ovh= *OVR; }
  elsif($cdnaseq =~ /stdin|^\-$/) { $ok=1; $ovh= *STDIN; }
  else { $ok= open(OVR,$cdnaseq); $ovh= *OVR; }
  die "bad -cdnaseq=$cdnaseq" unless($ok);
  my $id="none";
  while(<$ovh>) { if(/^>(\S+)/) { $id=$1; $cdnaseq{$id}=""; 
    if(m/ +(\S.+)$/){ my $h=$1; $h =~ s/\s*(len|cf|nt)=\S+//g; $cdnahead{$id}=$h; }
    } elsif(/\w/) { chomp; $cdnaseq{$id}.= uc($_); } } 
  close($ovh);

  cdna_bestorf() unless(@input);
}

foreach my $input (@input) {
  my $inh= *STDIN;
  $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
        : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
        : open($inh,$input);
  die "bad -input=$input" unless($ok);
  
  my ($nchanged,$ngene)= filter_gff($inh);
  warn "#findcds changed=$nchanged, ngene=$ngene\n" if $debug;
}


#..................

# in cdna_proteins.pm
sub bestorf_test
{
  my($bestorf,$nextorf) = @_;
  # need all of these from orfi->xxx ?
  # ($orfproti, $prostart5i, $proend3i, $nextorf)
  my ($bestprot)= orfParts($bestorf);
  my ($nextprot)= orfParts($nextorf);
  my $MinCDS= 3*$MINAA;

       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
  
  my $oki= ( $nextprot =~ /\w/ 
    && ($nextorf->{goodlen} >= $MinCDS)
    && ($nextorf->{goodlen}/$nextorf->{length} >= $MINGOOD)) ? 1 : 0;
    
  if($oki) {
    if( $bestprot ) {
      ## problem here for rev=0,1 missing huge partial rev for short fwd .. but get huge if just rev
      ## also problem case of bestorf.goodlen == nextorf.goodlen, but bestorf adds long NNNN gap
      my $arat = $nextorf->{goodlen} /  $bestorf->{goodlen}; # ok,> MINAA
      my $adiff= $nextorf->{goodlen} -  $bestorf->{goodlen};  
      # my $bestgap= $bestorf->{length} - $bestorf->{goodlen};
      # my $nextgap= $nextorf->{length} - $nextorf->{goodlen};
      # if($aadiff >=0 and $bestgap > 9 and $nextgap == 0)
      
      if( $arat > 1.0 and  $arat < $AA_cdna_GT_genome and $bestprot =~ m/XX/) { ##  $bestorf->{goodlen}/$bestorf->{length} < 0.99
        (my $bp= $bestprot) =~ s/X/./g; 
        # if(length($bp) > length($nextprot) { } # chomp some?
        return(1, $nextprot, $nextorf) if($nextprot =~ m/$bp/);
      }
      
      if( $arat >= $AA_cdna_GT_genome
        or ( $adiff >  0 and ($nextorf->{complete}==3 or $bestorf->{complete} <3) ) 
        or ( $adiff >= 0 and $nextorf->{complete}==3 and $bestorf->{complete} <3 ) )
       { 
        # ($bestprot, $prostart5, $proend3, $bestorf, $isrev, $cdnabest)= 
        #  ($nextprot, $prostart5i, $proend3i, $nextorf, $rev, $cdnain);
        return(1, $nextprot, $nextorf);
       }
    } else {
      #($bestprot, $prostart5, $proend3, $bestorf, $isrev, $cdnabest)=  ## return
      # ($nextprot, $prostart5i, $proend3i, $nextorf, $rev, $cdnain);
      return(1, $nextprot, $nextorf);
    }
  }
  return( 0, $bestprot, $bestorf);
}


# 201207 ** REPLACED by cdna_bestorf.pl and cdna_proteins.pm

sub cdna_bestorf
{
  my @cdnaid= sort keys %cdnaseq; my $ncdna= @cdnaid;
  
  warn "#REPLACED by cdna_bestorf.pl, 201207\n";
  warn "#findcds from cdnaseq:$cdnaseq n=$ncdna\n" if $debug;
  my $ngood=0;
  ## my $MINAA= $MINID_CDS || 40; 
  my $MINTR= $MINAA * 3;
  my @rev=($action =~ /rev/)? (1) : ($action =~ /fwd/)? (0) : (0,1);
  my $fullpart=($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";

# ADD maybe: detect joined/fused genes using 
#  1. cds span < ~ 1/2 trspan, leaving long utr, << annotate these cases, require aa complete/partial5
#  2. check the utr for long cds
 
  foreach my $id (@cdnaid) {
    my $cdnain= $cdnaseq{$id}; 
    my $clen= length($cdnain); 
    my $cdnabest= $cdnain;
    my($orfprot, $prostart5, $proend3, $bestorf, $utrorf, $orflen, $isrev, $aalen, $pcds, $compl)= (0) x 20;
    if($clen < $MINTR) { 
      # warn?
    } else {
    for my $rev (@rev) {
      if($rev==1){  $cdnain= revcomp($cdnain); } 
      my($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= getBestProt($fullpart, $cdnain); # , undef, $oldStart_b,$oldStart_e

      my $change=0;
      ($change,$orfprot,$bestorf) = bestorf_test($bestorf,$bestorfi);
      if($change) { 
        $isrev= $rev; $cdnabest= $cdnain; $utrorf= $utrorfi;
        ( undef, $prostart5, $proend3) = orfParts($bestorf);
      }
      
      }
    } # MINTR
 
    my $crev=($isrev==1)?"-":"+";  
    if($orfprot) {
      $orflen= $bestorf->{length};
      $pcds  = ($clen>0 && $orflen>0) ? int(100*$orflen/$clen) : 0;
      my $u1len= $prostart5 - 1; my $u2len= $clen - $proend3;
      
      if($nostopcodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
      $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
      $aalen-- if(substr($orfprot,-1) eq '*'); # annoyance
      $compl= $bestorf->{complete};
      $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
      ##? not bad if partial? if u1len or u2len == 0
      if($pcds < 35 or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
      elsif($pcds < 60) { $compl.="-utrpoor";  } #?? maybe change to use prostart5 OR protend3 > 35%? 40% ?
      if($utrorf) { my $ol=int($utrorf->{length}/3); $compl.="-utrorf$ol"; }
      $ngood++;
    } else {
      # print join("\t", $id, $clen, 0, 0, 0, 0, "na")."\n" unless($action =~ /fasta/);
    }
    
    if($action =~ /fasta/) {
      my $cdnah= $cdnahead{$id}||"";
      if($action =~ /all/ or $orfprot) {
      $orfprot =~ s/(.{60})/$1\n/g;
      # clen=$orflen/$clen,$prostart5-$proend3; ??
      print ">$id aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; $cdnah\n";
      print $orfprot,"\n"; 
      
      # FIXME for utrorf: find way to split cdna b/n 1st, 2nd orf: midway?
      if($utrorf) {
        my( $uprot, $ustart, $uend) = orfParts($utrorf);
        my $ulen=  $utrorf->{length};  
        my $ualen= length($uprot);  
        my $upcds  = ($clen>0 && $ulen>0) ? int(100*$ulen/$clen) : 0;
        my $ucompl= $utrorf->{complete};
        $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
        $uprot =~ s/(.{60})/$1\n/g;
        print ">$id.utrorf aalen=$ualen,$upcds%,$ucompl; clen=$clen; strand=$crev; offs=$ustart-$uend; flag=UTRprotein\n";
        print $uprot,"\n"; 
        }
        
      if($action =~ /cds|cdna/i) {
        my $tag=($action=~/cds/i)?"cds":"cdna";
        $cdnabest= $bestorf->{sequence} if($action=~/cds/i); # not last
        $cdnabest  =~ s/(.{60})/$1\n/g;
        print ">$id.$tag aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; \n";
        print $cdnabest,"\n";
        }
      }
    } else { 
      print join("\t", $id, $aalen, $pcds, $compl, $clen, $crev, $prostart5, $proend3, $orfprot),"\n"; 
    }

  }
  warn "#findcds from cdnaseq found $ngood / $ncdna\n" if $debug;
}


sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nr,$nchange)= (0,0,0);
  my $printpass=1;
  # $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP);
  my @generec=(); my @otherft=();
  
  ## add own header ?? version?
  my $version=VERSION;
  print <<"EOGFF";
##gff-version 3
# appl: genefindcds
# vers: $version

EOGFF
  
  while(<$inh>){
    unless(/^\w/){ next if(/^(##gff-ver|#n |$)/);  print and next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { print if $printpass; next; } 
    
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nchange += testgene(\@generec, \@otherft) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @otherft=();
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc;         
    } elsif($tb>0 and $te>0) {
      push @otherft, $rloc;         
    }
      
  }
  
  $nchange += testgene(\@generec, \@otherft) if(@generec);
  
  return ($nchange,$ng);
}


sub collect_overlaps
{
  my($gff)= @_;
  my ($nr,$nx)=(0,0);
  while(<$gff>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $tattr ||="";
    
    # if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    #^ drop passtypes for mrnatypes,exontypes
    next unless($typ =~ /^($introntypes)$/);

# ?? only Introns here? need splice-site calcs ; no ID= for introns...
    $nr++;
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$nr; }

    #? only kINTRONERROR_OVER and kINTRONERROR_INSIDE
    my($inb,$ine)= ($tb,$te); # full intron span, keep
    # if($intron2splice == kINTRON2SPLICE_OVER or $intron2splice == kINTRONERROR_OVER) 
    if(1) { # always intron here
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
      ($tb,$te)= ($s1b,$s1e);
      my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $rloc); }
      ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
    }  
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
    my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $rloc); } # $generec
  }
  
  warn"#collect_overlaps n=$nr\n" if $debug;
  return $nr;# return \%overlaps;
}



sub putgene
{
  my ($generec,$otherft,$flags)= @_;
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  my $cc= ($flags and $flags =~ /skip=/)? "#x." : "";
  $otherft ||= [];
  foreach my $ft (@$generec, @$otherft) { 
    if(ref $ft) { my @v= @$ft; 
      $v[8]=~s/$/;$flags/ if($flags and $v[2] eq "mRNA");
      print $cc.join("\t",@v[0..8])."\n" if(@v>4); 
      }
    }
  print "\n"; #?
}



# in cdna_proteins.pm
# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub _sort_over { # @[b,e] min-b, max-e
  return ($a->[0] <=> $b->[0]) || ($b->[1] <=> $a->[1]);
}

sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

sub _sortloc { # _sortgene
  #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

sub _sortinloc1 { 
  #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}
sub _sortinloc { 
  #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || (abs($a->[3] - $b->[3])>90 && ($a->[3] <=> $b->[3])) # begin1 <
      || ($b->[5] <=> $a->[5]) # score > ; should be start+=50 first then score
      || ($a->[3] <=> $b->[3]) # begin <
      || ($a->[4] <=> $b->[4]) # end <
      ;
}


sub intronoverlaps
{
  my ( $ref,$tb,$te,$to, $inmaxscore, $validonly)= @_;    # exon here
  my ( %didid,@overs,@ovok);
  return 0 unless($overlaplist->{$ref});
  $validonly ||= 0;
  
  my($tb1,$te1, $tb2, $te2)=(0) x 4;
  if(1) { # always ($intron2splice == kINTRONERROR_OVER)
    # input == exon, over= intron, test if overlap is splice end or internal
    ($tb1,$te1, $tb2, $te2)= 
      ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
  }
  
  $inmaxscore ||= 0; 
  my $pmaxscore= $inmaxscore * $pMAXSCORE;
  my($nover)= (0); ## NEED MAX inscore over all exons not just 1; global? reset per gene
  
  my($linb,$line, $llo,
     $okleft, $okright, $errleft, $errright, $errinside)= (0) x 10;
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    
    # sort, sameorient 1st, so bidir introns are less problem : NOT working
    @locs= ( (sort _sortinloc grep { $_->[6] eq $to } @locs),  (sort _sortinloc grep { $_->[6] ne $to } @locs) );
    
    foreach my $rloc (@locs) {
      #new# my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
      
      my ($lb,$le, $inscore,$lo,$oid,$inb,$ine)= @{$rloc}[3,4,5,6,9,10,11];  # intron here : 10-11 == full span, 3-4 == splice
       # FIXME now is intron ENDs
       # FIXME: need full intron lb,le and also splice ends
       # FIX-11Dec: use inscore/count to decide if valid cut: if validin.c=1000 and insidein.c=10, skipit.
       
      next if($didid{$oid.$lb.$le}++);
      #NO, count overs# next if($inb < $line && $ine > $linb); # overlap last intron, rev: skip
      # ($linb, $line, $llo)= ($inb, $ine, $lo);
      
      my $over  = ($tb <= $le && $te >= $lb) ? 1 : 0;      
      # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
      # my $inside= ($over and $tb <= $lb && $te >= $le) ? 1 : 0; ## intron inside  
      # $inside=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
 
      # if($inside) {  push @overs, [$lb,$le]; } else
      
      if($over) { #always($intron2splice == kINTRONERROR_OVER)   # test for splice reversed errors
        ($linb, $line, $llo)= ($inb, $ine, $lo);
        
        my $ok=0; # note: lb,le here are one splice end span of intron: 3 bp?
        my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;

        if($tb1 <= $le && $te1 >= $lb) { if($samestrand) { $ok=1; $okleft+=$inscore; } else { $ok= kINTRONERROR_OVER; $errleft+=$inscore; } } # $errt="rev" unless $samestrand
        elsif($tb2 <= $le && $te2 >= $lb) { if($samestrand) { $ok=2; $okright+=$inscore; } else { $ok=kINTRONERROR_OVER; $errright+=$inscore; } } # -1 or -2?
        elsif(($tb + 2*SPLICE <= $inb && $te - 2*SPLICE >= $ine) and $samestrand) { $ok= kINTRONERROR_INSIDE; $errinside+=$inscore; }  #-3
        
#         if($tb1 <= $le && $te1 >= $lb) { $ok= ($samestrand) ? 1 : kINTRONERROR_OVER; } # $errt="rev" unless $samestrand
#         elsif($tb2 <= $le && $te2 >= $lb) { $ok= ($samestrand) ? 2 : kINTRONERROR_OVER; } # -1 or -2?
#         elsif(($tb + 2*SPLICE <= $inb && $te - 2*SPLICE >= $ine) and $samestrand) { $ok= kINTRONERROR_INSIDE; }  #-3
                  # ^^ same strand inside  == retained intron err
        #old# elsif(($tb + 2*SPLICE <= $lb && $te - 2*SPLICE >= $le) and $samestrand) { $ok= kINTRONERROR_INSIDE; } 
        # else what?
        
        if($ok < 0) {
          $nover++;
          # $enderr |= $ok; # 1,2 or 3=both
          if($inscore < $pmaxscore) { } ## $inmaxscore * $pMAXSCORE)  # skip weak cases
          elsif($ok == kINTRONERROR_INSIDE) { push @overs, [$inb,$ine, $ok]; } # need error type also: retained vs strand-err
          else { push @overs, [$lb,$le, $ok]; } # need error type also: retained vs strand-err
          #? need ok1, ok2 for both ends of exon to say if 1 end is ok?
        } elsif( $ok > 0) {
          $nover++;
          if($inscore > $inmaxscore)
            { $inmaxscore= $inscore; $pmaxscore= $inmaxscore * $pMAXSCORE; }
          # $endok |= $ok; # 1,2 or 3=both; want count of each end ok?
          #? push @overs, [$lb,$le, $ok]; # UPDATE 11Dec, see below 
          push @ovok, [$lb,$le, $ok, $inscore]; #? return which end of exon is supported by intron (or both)
          # push @spliceok, ($ok == 2) ? (($to eq "-") ? $tb : $te) : (($to eq "-") ? $te : $tb);
        }
        # next;
      } 
        
      # push @overs, [$lb,$le] if ($over);
      }
  }
  
  ## add for annotations, no changes
  if($validonly >= 2) {  ## is allowalts synonym of this now?
    # my $nok=@ovok; my $nerr=@overs; 
     ## dont want count of introns here, want to know if either/both ends of exon are supported or not
     ## and if there is internal intron kINTRONERROR_INSIDE
     
    my $okOrErr = ($nover==0) ? 0 : ($okleft >= $errleft and $okright >= $errright and $errinside < $inmaxscore) ? 1 : -1;
    # my $inannot = ($nover==0) ? 0 : join ",", $okOrErr, $okleft,$okright, -$errleft,-$errright,-$errinside;
    # exon annot:  ;inok=[1/-1],20,10,-2,-20,-3; want all this data?  or inok=ngood,-nbad | inok=-nbad,ngood
#     my $inannot= ($okOrErr == 0) ? 0 
#         : ($okOrErr>0) ? join ",", ($okleft+$okright), (-$errleft-$errright-$errinside) 
#         : join ",", (-$errleft-$errright-$errinside),($okleft+$okright);
    my $inannot= ($okOrErr == 0) ? 0 
        : ($okOrErr>0) ? join ",", $okleft.'l',$okright.'r', -$errleft.'l',-$errright.'r',-$errinside.'i'
        : join ",", -$errleft.'l',-$errright.'r',-$errinside.'i',$okleft.'l',$okright.'r';
    return($okOrErr, $inannot); ##, $okleft, $okright, $errleft, $errright, $errinside); 
    }
    
  return \@ovok if($validonly);
  return 0 if(@overs > 1 and $allowalts); # dont force cut where we want alts and they exist
   
  # ** FIXME: this makes bad cuts, likely where 2+ alt-introns are overlapped
  my @opens=();
  if(@overs) {
    my $okover= scalar(@ovok); #** Need ovok in overs to let bidir good superceed bad/rev introns
    my($bb,$be)= ($tb,$te); my $errlast= 0; my($llb,$lle)=(0,0);
    #OFF. @overs= sort _sort_over @overs; # NOT NOW, see above sort
    
    foreach my $ab (@overs) {
      my($lb,$le,$errcode)= @$ab;
      ## .. for $errcode == kINTRONERROR_OVER, strand err at exon splice; need to chop off entire exon ?
      # next if($lb >= $llb and $le <= $lle); # skip inside alt-intron
      next if($lb < $lle and $le > $llb); # skip overany alt intron
      
      ## .. for $errcode == kINTRONERROR_INSIDE
      if($le < $bb) {  }
      elsif($lb <= $bb && $le > $bb) { $bb= $le+1; }
      elsif($lb < $be) {  #  && $le > $be
        my($b1,$e1)= ($bb,$lb-1);
        if( 1 or $errcode < 0) { push @opens, [$b1,$e1,$errcode,$okover]; }
        $bb= $le+1; 
        } 
      elsif($lb > $te) { last; } #?
      ($llb,$lle)= ($lb,$le);
      $errlast= $errcode;  #?? if($errcode==kINTRONERROR_INSIDE or $errlast==0); # upd 11dec
      last if($bb >= $te);
    }
    if($bb < $te) { push @opens, [$bb,$te,$errlast,$okover]; } # add end point
    return \@opens;
    # return (\@opens, \@spliceok);
  } else {
    return 0; ## [[$tb,$te]];
    # return ( [], \@spliceok);
  }
}


sub intron_error_cut
{
  my($exongff)= @_;
  
  my $xdebug=($debug>1)?"/xdebug":"";
  
  my @exnew=(); 
  my $changenote=""; 
  my $changed= 0; my $xi=0;  my $droperr=0;
  my $inmaxscore=0;
  foreach my $ex (@$exongff) {
    my $valids = intronoverlaps( @{$ex}[0,3,4,6], 0, 1 ); 
    if(ref $valids and @$valids > 0) {
      foreach my $xloc (@$valids) {
        my($xb,$xe,$errcode,$score)= @$xloc;
        $inmaxscore= $score if($score>$inmaxscore);
        }
    }
  }
  
  foreach my $ex (@$exongff) {
  
    # ** FIXME: this makes bad cuts, likely where 2+ alt-introns are overlapped
    my $exonfix = intronoverlaps( @{$ex}[0,3,4,6], $inmaxscore ); 
    $xi++;
    # result= [array] of [$bb,$te,$errcode] ; errcode == kINTRONERROR_OVER/bad splice, kINTRONERROR_INSIDE/retained in

    # 2011Dec: DAMN, reverse(@$exonfix) for strand "-" : see below
    # 2011Dec: DAMNN, got dupl xcut locations, all with ? xdrop=-1/1009/12.2

# orig:
# scaffold_2      caca11r39cuf8   mRNA    6247942 6251164 225     +       .       ID=caca11r39cuf8_Gsc2g4225t1;
# scaffold_2      caca11r39cuf8   exon    6247942 6248192 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=1;
# scaffold_2      caca11r39cuf8   exon    6248434 6248780 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=2;
# scaffold_2      caca11r39cuf8   exon    6248876 6251164 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;
#.. after xcut ..  intronfix=1,xcut:-3/792/3.1,xdrop:-1/498/3.2,
#  # Thecc1EG007301t1	rna8b:r8caca11r39cuf8_Gsc2g4225t1  << BAD xcut, DUPL exons,CDS; from overlapped introns
# scaffold_2	caca11r39cuf8	exon	6247942	6248192	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=1;
# scaffold_2	caca11r39cuf8	exon	6248434	6248780	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=2;
# scaffold_2	caca11r39cuf8	exon	6248876	6249667	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;;xcut=-3/792/3.1
#   ^^^^^^ bad exon1? or good
# scaffold_2	caca11r39cuf8	exon	6248876	6251164	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;;xdrop=-1/498/3.2
#   ^^^^^^ bad exon2 -- xb should be after 6249667
    
    if(ref $exonfix and @$exonfix > 0) {
      my $fixerr=0;
      my @exadd;
      my $xj=0;
      my($lxb,$lxe,$didx)=(0,0,0);
      foreach my $xloc (@$exonfix) {
        my($xb,$xe,$errcode,$hasok)= @$xloc;  $xj++;
        my $xw= 1+$xe-$xb;
        $fixerr++ if($xw < $MINEXON);
        
        if($errcode == kINTRONERROR_OVER) {
          # drop exon span? but annotate what?
          my $xan="xdrop=$errcode$xdebug/$xw/$xi.$xj";
          $changenote .= "$xan,";
          # careful: dont drop 1st w/ intron error .. keep as last exon before error?
          # exons are 5' > 3' order and prot/CDS probably only in 5' before err
          # .. but doesnt look wrong w/o this.
         if($didx) { # upd 11dec
            my @xx= @$ex;  $xx[3]= $xb; $xx[4]= $xe;
            $xx[8] =~ s/$/;$xan/;
            push @exadd, \@xx;
            ($lxb,$lxe)=($xb,$xe); $didx=0;
          } elsif($droperr==0 and $hasok) {
            my @xx= @$ex;  #? $xx[3]= $xb; $xx[4]= $xe; #< change to this or not?
            ## lxe >> $xx[3]= $xb; $xx[4]= $xe;
            $xx[8] =~ s/$/;$xan/;
            push @exadd, \@xx;
           }
          $droperr++;
          
        } elsif($errcode == kINTRONERROR_INSIDE) {
          my $xan="xcut=$errcode$xdebug/$xw/$xi.$xj";
          $changenote .= "$xan,";
          my @xx= @$ex;  $xx[3]= $xb; $xx[4]= $xe;
          $xx[8] =~ s/$/;$xan/;
          push @exadd, \@xx;
          ($lxb,$lxe)=($xb,$xe); $didx++;
        }
      }
      
      if($fixerr>0) { 
        push @exnew, $ex;
      } else {
        $changed++; push @exnew, @exadd;
      }
      
    } else {
      push @exnew, $ex;
    }
  }
  
  $changenote =~ s/=/:/g;
  $changenote= "intronfix=$changed,$changenote" if($changed);
  # return $badintrons == long, no support
  return ($changed, $changenote, @exnew);
}


# in cdna_proteins.pm; revised
sub getBestProt
{
  my($ptype, $cdna)= @_;
  ## my $oflags=($ptype =~ /dropnnn/i)?"dropnnn":"";
  my ($longorf,$longfull,$orfs) = getAllOrfs($cdna,"fwd");  # ,$oflags
  ## my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  return getBestProtOfOrfs($longorf,$longfull,$orfs, @_);
}
  
  # getBestProt fix this to return longest full prot, and longest partial (if longer)
  # .. test which is best.
  # FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
  #   my $longorf= $longest_orf_finder->get_longest_orf($cdna); # == hash
  # FIXME3: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf

sub getBestProtOfOrfs
{
  my($longorf,$longfull,$orfs,
     $ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my($orfprot,$prostart5,$proend3)=("",0,0,"");
  my ($utrorf,$utrosize)=(undef,0);
   
  if(ref($longorf)) {
    # $orfprot= $longest_orf_finder->get_peptide_sequence();
    # ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();

    my $lookmore= ($ptype =~ /long/)?0:1;
    if($samecds and $oldStart_b > 0) { # may not be right yet.
      my($samestartorf);
      if($ORF_FULLvPART <= 0.8) {
      ($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
      } else {
      ($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
      }
      if(ref $samestartorf) { $longorf= $samestartorf; $lookmore=0; } #NOT# else { return (undef); } # not found == no change, here
    } 
    
    if($lookmore and $longorf->{complete} < 3 and ref($longfull) ) {
      my $keylen=($USEGOODLEN)?"goodlen":"length";
      my $lsize= $longorf->{$keylen}; # was {length}
      my $fsize= $longfull->{$keylen};
      
      # FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
      if($fsize >= $ORF_FULLvPART * $lsize) {
        if(ref($exongff)) {
        my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff( $exongff, orfParts($longorf));
        my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff( $exongff, orfParts($longfull));
        $longorf= $longfull if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
        } else {
	      $longorf= $longfull;
	      }
      }
    }
    
    ($orfprot,$prostart5,$proend3)= orfParts($longorf);
    ## ($orfprot,$prostart5,$proend3)= ($longorf->{protein},$longorf->{start}, $longorf->{stop});

    ($utrorf,$utrosize)= getUtrOrf($longorf, $orfs, $cdna);  
  }

  return($orfprot, $prostart5, $proend3, $longorf, $utrorf);
}

# in cdna_proteins.pm
sub getUtrOrf
{
  my($longorf, $orfs, $cdna)= @_;
  my ($utrorf,$utrosize)=(0,0);
  my $lsize= $longorf->{length};
  my $lgood= $longorf->{goodlen};
  my $cdnalen= length($cdna);
  my $MINUTRORF=300; # global
  return($utrorf,$utrosize) unless(($cdnalen - $lsize >= $MINUTRORF));  # test even if lsize/cdna > 60% ?  
  if($cdnalen and $lsize/$cdnalen < 0.66) { 
    my($lb,$le)= ($longorf->{start},$longorf->{stop});
    foreach my $orf (@$orfs) {
      my($ob,$oe,$ogood,$osize)= ($orf->{start},$orf->{stop},$orf->{goodlen},$orf->{length},);
      if(($ob > $le or $oe < $lb) and ($ogood>=$MINUTRORF or $ogood > 0.5*$lgood) and $ogood>$utrosize) {
        $utrorf= $orf; $utrosize= $ogood;        
        }
      }
    }
  return($utrorf,$utrosize);    
}
  


# in cdna_proteins.pm; revised as getCDSgff2
sub getCDSgff
{
  my($exons,$orfprot,$prostart5,$proend3,$cdnalen) = @_;
  # FIXME: need trlen= length(cdnain) for -cdna, and/or use gmap qlen= tag
  $cdnalen ||= 0;
    
  my ($cds5,$cds3)=(0,0);
  my @cds= ();
  my @utr= ();
  ## for phase; need reverse @exons
  my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3)= (0) x 10;
  
  $cdna1= 0; # was 1; # offby1 at end?
  $nt_length= 0; # $prostart5 % 3; #??
  
  # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
  # .. must make next exon(if exists) part of CDS stop
  
  foreach my $exon (@$exons) {
    my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
    
    ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
    my $xd= abs($xend3 - $xend5); # ?? +1 for width
    
    $cdna1++; # add 1 here, not end loop 
    my $cdna2= $cdna1 + $xd; # ?? +1 for width
    # ** offby1 here ?? YES, need <=, >= to get full CDS stop,start split by intron; see below cdna1=cdna2+1
    #OLD.if($cdna1 < $proend3 and $cdna2 > $prostart5) 
    if($cdna1 <= $proend3 and $cdna2 >= $prostart5) 
    { # overlap
                
      my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
      my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      
      my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
      my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
  
      my $elength = 1 + abs($c3 - $c5);
      $nt_length  += $elength;
      $inc3        = $nt_length % 3;
      $inc5        = ($elength - $inc3) % 3; # only care about this one
      # $frame       = ($c5 + $inc5) % 3;
      if ($inc5 == -1) { $inc5 = 2; }
      
      my $phase= $inc5; # is this right?
      
      my($cb,$ce)= ($c5 > $c3) ? ($c3,$c5): ($c5,$c3); #? rev patch
      my $rloc= [$ref,$src,"CDS",$cb,$ce,".",$xo,$phase,"Parent=$gid",$gid]; 
      push @cds, $rloc;
     
      $cds5=$c5 if($cdna1 <= $prostart5);
      $cds3=$c3 if($cdna2 >= $proend3);
      
    } elsif(1) { # $addutr
      my $d5= ($cdna1 >= $proend3) ? 0 : $proend3 - $cdna1; # pos
      my $u5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      my $d3= ($cdna2 <= $prostart5) ? 0 : $prostart5 - $cdna2; # neg
      my $u3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;

      my($ub,$ue)= ($u5 > $u3) ? ($u3,$u5): ($u5,$u3); 
      my $up= ($cdna1 < $prostart5) ? "five" : ($cdna2 > $proend3) ? "three" : "odd";
      if($cdna1 < $prostart5) { $nu5++; } elsif($cdna2 > $proend3) { $nu3++; }
      my $rloc= [$ref,$src, $up."_prime_utr",$ub,$ue,".",$xo,0,"Parent=$gid",$gid]; 
      push @utr, $rloc;
    }
    
    ##$cdna1= $cdna2+1;  # is this off-by-1 now? yes, dont +1 here, do above
    $cdna1= $cdna2;  
  }       
  
  
  my $trlen= ($cdnalen>$cdna1) ? $cdnalen : $cdna1; # if($cdnalen> $cdna1 or > 0) ??
  my $aalen=length($orfprot); 
  $aalen-- if(substr($orfprot,-1) eq '*');
  my $clen= $aalen * 3; # can be off by -1,-2 here. 
  my $ap=int(100 * $clen/$trlen);

## .. add this prot qual test, also test orig prot.. esp for augustus X inner stops
##        if($istop < 1 and $id =~ /AUG/) { $istop= index($aa,'X'); }
##        if($istop > 0 and $istop < $al-1) { $astat="partialinner"; }
  
  my $mattr="cxlen=$clen/$trlen;aalen=$aalen,$ap%";
  if($orfprot) {
    my $p5= (substr($orfprot,0,1) eq 'M')?0:1;
    my $p3= (substr($orfprot,-1,1) eq '*')?0:1;
    my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
    $mattr.= ",$prostat;protein=$orfprot";
  }
  
  $mattr.= ";cdsoff=$prostart5-$proend3"; #? as per ;utroff=$ustart-$uend
  $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
  ## mattr keys: cxlen,aalen,protein,utrx
  
  # FIXME: resort @cds by loc, not reversed. : let caller do
  # @cds = sort _sortgene @cds;

  ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
  ## ($cdslong, $attrL, $pcodeL, $maxutrL)
  return (\@cds,$mattr,$ap, _max($nu5,$nu3), \@utr); 
}



sub getcdna {
  my($exons, $asexons, $expand, $expend)= @_;
  my $cdna= ""; my $cdnalen=0; my @asexons=();
  my $lstop= 1;   
  my $doexp=(defined $expand && $expand>0)?1:0;
  if($doexp) { $expend=$expand unless($expend); }
  my $nx1= @$exons - 1;
  foreach my $j (0..$nx1) {   # @$exons
    my $ft= $exons->[$j];
    my($ref,$start,$stop,$strand,$phase)= @{$ft}[0,3,4,6,7];
    my $rev=($strand eq "-")?1:0;
    # if($samecds and $phase>0 and $cdnalen==0) { if($rev) { $stop-=$phase; } else { $start+=$phase; } }
    # ^^ not here, use phase in get_orfs ..
    $cdnalen += 1+$stop-$start; # add even if no dna
    my($xstart,$xstop)= ($start, $stop);
    if($doexp) {
       my $nstart= ($j<$nx1)? $exons->[$j+1]->[3] : 999999999;
       my $xp=($j==0)? $expend : $expand;    $xstart= _max(1, _max($lstop,$xstart - $xp));
       $xp=($j==$nx1)? $expend : $expand; $xstop=  _min($nstart, $xstop + $xp);
    }
    my $exondna  = get_dna( $dnasequence, $ref, $xstart, $xstop);
    $lstop= $stop;
    next unless($exondna);
    $exondna = uc($exondna); # always?
    $exondna = revcomp($exondna) if($rev);  # $exondna = reverse $exondna;  $exondna =~ tr/ACGTacgt/TGCAtgca/;
    $cdna .= $exondna;  
    push @asexons, $exondna if($asexons);
  }
  $cdnalen= length($cdna) if($cdna);
  return($cdna,$cdnalen,\@asexons);
}  


sub stripOldAnnot
{
  my ($mrna,$oldtags)= @_;
  return unless(ref $mrna and $mrna->[8]);
  $oldtags= [qw(cxlen aalen protein cdnabest cdnaorf aautrlen utroff utrprot ocds oaaln inqual intronfix xcut xdrop)]
    unless(ref $oldtags);
  my($an)= $mrna->[8]; my $oldan= $an;
  my $olds= join('|', @$oldtags);
  if( $an =~ s/;($olds)=[^;\n]+//g ) { $mrna->[8]= $an; } # or not?
  return ($an,$oldan); 
}


sub testgene
{
  my($generecIN, $geneother)= @_;
  my $addattr="";
  my $changed= 0;
  
  # xFIXME : check/remove existing CDS; compare to new 
  # FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
  # FIXME3: -samecds not much use in test w/ augustus calls.
  # FIXME4: stop changes that replace all CDS with completely new CDS, happens for transposon-spans
  #        .. not quite -samecds, but -nocompletelydifferentcds
  # FIXME5: -samecds needs to use phase0 : partial5 from AUG uses this (eg. at scaf startpos=1)
  # FIXME6: chimera, multimap IDs have 2ndary tag: _C[12] _G[2..n] ; for cdnain fix this
  
  
  my @oldCDS = grep{ $_->[2] eq "CDS" } @$generecIN;
  my @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
  my($mrna)  = grep{ $_->[2] eq "mRNA" } @generec;
  my @exongff= grep{ $_->[2] eq "exon" } @generec;

  unless($mrna and @exongff > 0) {
    if($USE_CDSEXONS and $mrna and @oldCDS > 0) { # or caller can s/CDS/exon/ easily, or duplicate like this
      foreach my $oc (@oldCDS) { my @doc= @$oc; $doc[2]="exon"; push @exongff, \@doc; }    
      push @generec, @exongff;
    } else {
      putgene($generecIN, $geneother,"err=Missing-mrna-exon"); # flag="err=Missing-mrna-exon"
      return 0;
    }
  }
  
  # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  my $gstrand= $mrna->[6]; # FIXME for "." ; add CDS bestpro strand (always +??)
  # *** FIXME 2011Dec : 
  my $fixstrand=0;
  my(@genefwd, @generev);
  if($gstrand eq ".") {
    $fixstrand=1; # need to test bestpro both strands
  }
  
  my ($changed0,$changed1)=(0,0);
  for( my $ifix= 0; $ifix <= $fixstrand; $ifix++) { 
    # BIG loop testing both ways, save protfwd, protrev
    $changed=0; $addattr=""; # clear for loop step
    
    if($fixstrand) {
      $gstrand=($ifix == 1) ? "-" : "+";
      ## must clone generec for save @genefwd, @generev
      @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
      @generec= map { my @xnew= @$_; \@xnew; } @generec; # clone all
      map{ $_->[6]= $gstrand } @generec;
      ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
      @exongff= grep{ $_->[2] eq "exon" } @generec;

      if( $ifix == 1 ) {  @generev= @generec; } 
      else {  @genefwd= @generec; }
    }  
  
  my $geneid= ($mrna->[8] =~ m/ID=([^;\s]+)/)? $1 : ""; # make one up?
  (my $geneidfix= $geneid) =~ s/_[CG]\d+$//; # chimera/splitgene _C[12] and multimap _Gnnn id tags
  
  my $oldprot= ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";
  @exongff= reverse @exongff if($gstrand eq "-");
  
# FIXME5: -samecds needs to use phase0 : partial5 from AUG uses this (eg. at scaf startpos=1)
  my ($oldStart_b,$oldStart_e)= (0,0); # only for -samecds 
  if(@oldCDS and $samecds){ @oldCDS= sort _sortgene @oldCDS; # dang
    my $cstart= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3]; 
    foreach my $ex (@exongff) { my($b,$e)= ($ex->[3], $ex->[4]);
      if($b <= $cstart and $e >= $cstart) { ($oldStart_b,$oldStart_e)=($b,$e); last; }
    }
  }
  
#  FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
#  FIX3? add here option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
#  # my ($inchanged, $infixnote, $badintrons, @exoninfix) = intron_error_cut( \@exongff);
#  # if($badintrons) { putgene( $generecIN,"skip=1;badintron=$badintrons"); return 1; } 
#  FIX? for gaps, chomp off end gaps NNN of getcdna() ? or remove from orfprot
  
  my($cdna,$cdnalen)= getcdna( \@exongff);

  # add samecds opt using oldStart > need start,stop of 1st CDS exon
  my($orfprot, $prostart5, $proend3, $bestorf, $utrorf)= 
    getBestProt("partial", $cdna, \@exongff, $oldStart_b,$oldStart_e);
  
  
#... loop here for intron overlaps >>>>>>>>>>>>>>>>>>>> 
    ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
    # revise : option to annotate errors, but no change to gff.
  my($inoknote, $inok,$inerr,$inzip)=("",0,0,0); # use below with cdnain checks
  
  if($hasintrons and !$DO_INFIX) {
    foreach my $ex (@exongff) {
      my ($okOrErr, $inannot) = intronoverlaps( @{$ex}[0,3,4,6], 0, 2 ); 
      if($okOrErr == 0) { $inzip++; }  
      else {
        if($okOrErr>0) { $inok++; } elsif($okOrErr<0) { $inerr--; } # should this be inerr++ ? or add key?
        # @{$ex}[8] =~ s/$/;inok=$inannot/;  # drop this?
        } #? ## need also add to mRNA summary inqual= 
      }
      # my $inoknote="inqual=$inok/$inerr/$inzip"; # or inqual=9ok,3bad,2none ? or put -err first if > ok; 
      my $incode= int (100 * ($inok + $inerr) / ($inok + abs($inerr) + $inzip)); # can be -
      $inoknote="inqual=$incode," . ((abs($inerr) > $inok) ? "$inerr/$inok/$inzip" : "$inok/$inerr/$inzip");
      $addattr .=";" if($addattr); $addattr .= $inoknote;   

  # FIXME: inqual : overbestgenes1 has better? version; switch to that?
  #   .. add inqual score for DO_INFIX also?
  #      $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
  #      $flags .= "ints=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal") .",$iflag;";
            
  } elsif($hasintrons and $DO_INFIX) {  
   
    my ($inchanged, $infixnote, @exoninfix) = intron_error_cut( \@exongff);
    if($inchanged and @exoninfix) {
  
      if($gstrand eq "-") {   # 2011Dec fix **
        my @xs= reverse sort _sortloc @exoninfix;
        @exoninfix= @xs;
      }
      
      my($cdnafix, $cdnafixlen)= getcdna( \@exoninfix);
      my($fixprot, $fixprostart5, $fixproend5)= getBestProt("partial", $cdnafix, \@exoninfix);
      # test also: poor=(utr5>2 or utr3 >2) from getCDSgff()
      
      ## **?? replace infix when xdrop= error even if fixprot == orfprot
      if( (length($fixprot) > 2 + length($orfprot))
         or ($infixnote =~ /xdebug/)
         or ($infixnote =~ /xdrop=/ and length($fixprot) >= length($orfprot)) ) {
        # annotate this
        $addattr .=";" if($addattr); $addattr .= $infixnote; $changed++;
        ($orfprot, $prostart5, $proend3)= ($fixprot, $fixprostart5, $fixproend5);
        # @generec : need to replace exons in @generec here
        my @oldexon= map { my @xnew= @$_; \@xnew; } @exongff; # clone all
        @generec= grep{ $_->[2] ne "exon" } @generec;
        push @generec, @exoninfix;
        if($debug>1) { foreach my $ex (@oldexon) {
          $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; push @generec, $ex;
          } }
        @exongff = @exoninfix; # preserve oldexon if debug ??
        my($newstart,$newstop)= (0,0);
        foreach my $ex (@exoninfix) {
          my($b,$e)= @{$ex}[3,4];
          $newstart= $b if($newstart == 0 or $b < $newstart);
          $newstop= $e if($newstop == 0 or $e > $newstop);
          }
        $mrna->[3] = $newstart; $mrna->[4] = $newstop;
      }
    }
  } # hasintrons
  
#... end loop for intron overlaps <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#... test input cdna ..............  
  ## FIXME: cdnain is strand-less ; rev for fixstrand ?
  ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
  
  my $cdnain  = ""; 
  my $cdnainlen= 0; my $cdnainNotPrimaryId= 0; my $cdnainIsBest=0;
  my $cdnagfflen= $cdnalen;
  if($cdnaseq) {  # add 2011
    # $cdnain= $cdnaseq{$geneid} || $cdnaseq{$geneidfix}; 
    $cdnain= $cdnaseq{$geneid};
    if( not $cdnain and ($geneid ne $geneidfix) ) { $cdnain= $cdnaseq{$geneidfix}; $cdnainNotPrimaryId=1 if($cdnain); }
   } 
   
  if($cdnain) { 
    $cdnainlen= length($cdnain);

# FIXed: some bug here preventing Much Better cdnain orf from replacing mapped orf
# FIXED: ** userev index(cdna) test was bad; check both orfs; and/or check mrna sense=-1 annot from gmap >> trrev of gff
#   is bestorf_test() ok? yes
#   .. in cases of < 50% mapped cdnain; could be as simple as userev test, may need to do both strands

# FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
#   mapped transcript is chimera at same locus giving clue to mashup.  happens more for velvet low kmer
#   where smaller fragments of 2+ genes are mashed together.   Need to check mRNA chimera tags for clues?

## Dang BUG2: chimera-only? cdnabest prot: protein=*MTGHYYYE..HFEY*  ^* is bogus stopcodon, from where ??
##   also mayb only with aautr ; also -nostopcodon removes * prefix, but doesnt remove all end *
## FIXED below at $oldprot.'*'; 

# FIXME3: chimera, _C1 and _C2 getting same protein for cdnainIsBest, including bad cases of C1=95%, C2 = 5% mapped
#  .. here or later allow only 1 of pair to have cdnainIsBest.. should be based on mapped CDS best match to protein
#     not just longest mapping portion (which often is bogus gene join as UTR).

## FIXME4: strand reversal, $userev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
#      .. or need to reverse strand of output gff
#   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
#   see fixstrand..

## FIXME5: new bug: Now gap-cdna with XXX protein winning over cdnain with full protein, same length. 
#   adjust goodlen test? problem may be index(cdnain,cdna) fails due to XXXX garbage

# new # need test both ** YES, this cures problems
# ** BUT need to know if best cdnain is same/rev of gff cdna .. FIXME4
# .. revert to old, bad; need cdnain strand same as input gff/genome cdna for fully mapped cdna
# .. but test both strands for split-maps (chimera), partial maps (coverage < 60?), and oneexon trs.

    use constant CIN_SAMESTRANDasMAP => 1;
    use constant CIN_BOTHSTRANDS => 0;

    my $cdnain_isrev=0; # doesnt mean samestrand-as-map, but just that cdnain was reversed.
    my $cstrandIsOpposite= 0; # 1 if antisense opposite, 0 == sense same as genome map
    my($orfproti, $prostart5i, $proend3i,$bestorfi, $utrorfi);
    my $cdnaindex= -1;
    # my $cdnainSense=1; # or -1 if opposite of mapped strand
    # from gmap, "sense=-1" means cdnain orient was opposite of mapped strand
    # ?? what of utrorf on rev strand of bestorf ??
    
    ORFCDNAIN: { 
    my $cdnainrev= revcomp($cdnain);
    # replace exon count with valid intron count?
    # my $validgffstrand= ( $inok > $inerr) ? 1 : ($inok+$inerr == 0 and @exongff > 1) ? 1 : 0;
    my $validgffstrand= ( @exongff > 1 ) ? 1 : 0;
    # cancel if (abs($inerr) > $inok) 
    
    # problem here, some chimera cdna parts can index cdnain both ways; ie cdnain is mashup of 1, or 2 near same.
    # .. Or, best of both is reverse of this cdna part, but on other cdna part.. case of rev-genes-join
    # .. cdna part does not map to bestorf part, so cant correct cdsgff from this
    
    my $cstrand_gstrand= 0;
    my $cdnagffgoodlen= $bestorf->{goodlen}; # not cdnagfflen
    my $cdnagood= $cdna; # only for index cstrand_gstrand test
    if($cdnagffgoodlen < $cdnagfflen) {
      $cdnagood =~ s/^N+//; $cdnagood =~ s/N+$//; 
      my $xi= index($cdnagood,"NN"); 
      if($xi>0) { 
        my $xe= rindex($cdnagood,"NN");
        if($xi>$cdnagffgoodlen/2) { $cdnagood= substr($cdnagood,$xe+2); }
        else { $cdnagood= substr($cdnagood,0,$xi); }
      }
    }    
    if( ($cdnaindex= index($cdnain,$cdnagood)) >=0 ) {  $cstrand_gstrand= 1; }
    elsif( ($cdnaindex= index($cdnainrev,$cdnagood)) >=0 ) { $cstrand_gstrand= -1; }
    
    ## here use goodlength tests for cdnain,cdnagff
    my $testboth= 0; #  2=best of both 
    if($fixstrand) { } # DONT test both in this strand loop
    elsif( $validgffstrand and $cdnainlen < 1.5 * $cdnagffgoodlen) { $testboth=0; }
    elsif( $cdnainNotPrimaryId or ($cdnainlen >= 1.5 * $cdnagffgoodlen) or not $validgffstrand) { $testboth= 2; }
    elsif( $cstrand_gstrand == 0 ) { $testboth= 2; } #? or leave 0
    
    if(CIN_SAMESTRANDasMAP) { }   
    elsif(CIN_BOTHSTRANDS) { $testboth=2; }

    my ($longorf,$longfull,$orfs);
    if( $testboth == 2 ) {    
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
      my ($rlongorf,$rlongfull,$rorfs) = getAllOrfs($cdnainrev,"fwd"); # ,"dropnnn"
      ($cdnain_isrev)= bestorf_test($longorf,$rlongorf);
      if($cdnain_isrev) {
        # $cstrand_gstrand= -1 if($cstrand_gstrand == 0);
        $cstrandIsOpposite=1 if($cstrand_gstrand == 1);
        $cdnain= $cdnainrev;
        ($longorf,$longfull,$orfs)= ($rlongorf,$rlongfull,$rorfs);
      } else {
        # $cstrand_gstrand=  1 if($cstrand_gstrand == 0);
        $cstrandIsOpposite=1 if($cstrand_gstrand == -1);
      }
    } elsif( $cstrand_gstrand == -1 ) {
      $cdnain_isrev=1;
      $cdnain= $cdnainrev;
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
    } else {  ##  $cstrand_gstrand == 1
      $cdnain_isrev=0;
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn" 
    }
    
    ($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= 
      getBestProtOfOrfs($longorf,$longfull,$orfs, 
        "partial", $cdnain, \@exongff, $oldStart_b,$oldStart_e);
    }
    

    if($orfproti ne $orfprot) { 
      ($cdnainIsBest,undef,undef) = bestorf_test($bestorf,$bestorfi);
       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
       ## if($orfprot =~ /XXX/) { (my $op=$orfprot) =~ s/X/./g; $cdnainIsBest=1 if($orfproti =~ /$op/); }
       
      # FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
      # add other checks here if cdnain better than genome-mapped cdna .. esp if cdnain == chimeric mapping
      #  versus partly unmapped due to genome gap.
      
      ## reasons cdnain is best: gaps in genome mapping, NNNs or indels, versus cleaner transcript
      ## reasons not best: chimeric mapping over same locus or nearby; reversed introns over longer cdnain mapping, ..
      ## .. maybe not .. leave this off for now.
#             
#       # my $cdnabestExplained= bestorf->goodlen/length < cdnainorf->goodlen/len
      my $cdnainProblems = $cdnainNotPrimaryId or ($hasintrons and $inerr != 0) or ($mrna->[8] =~ /chimera/);
      
      if($cdnainIsBest and $cdnainProblems) { 
      
        if($cdnainNotPrimaryId) { ## and $mrna->[8] =~ /chim\d=(\w+[^;\s]+)/
          # cancel if $pctmapped < 25, < 33? <50?
          my $ingood= $bestorfi->{goodlen} || 1;
          my $bgood = $bestorf->{goodlen};
          
          $cdnainIsBest= 0 if($bgood/$ingood < 0.25); # bad test.. really need other part's stats here
          ## and/or cancel if mapped CDS are missing or small fraction of cdnain?? need from getCDSgff() or @oldCDS?
          
        }
        
#         #?? this will cancel some good changes along w/ bad; add annot for scoring instead?
#         if($hasintrons and $inerr != 0) { $cdnainIsBest= 0; }
#         elsif($mrna->[8] =~ /chim\d=(\w+[^;\s]+)/) { 
#           # pick out chimera span and see if near this one.. chim1=scaffold00024:99418-99575:.
#           use constant CHIOK_MIN => 25000;
#           my $cloc= $1; my($cr,$cb,$ce)= $cloc =~ m/(\w+):(\d+).(\d+)/;  
#           my($thisr,$thisb,$thise)= ($mrna->[0], $mrna->[3], $mrna->[4]);
#           # .. maybe not always cancel; add flag instead? ..
#           $cdnainIsBest= 0 if($thisr eq $cr and (abs($thisb - $ce) < CHIOK_MIN or abs($thise - $cb) < CHIOK_MIN));
#         }
        # elsif($cdnainNotPrimaryId) {} #? what check?
      }
      
    }
      
    if($debug) {  #  and not $cdnainIsBest; debug set note; add orig cxlen;aalen here if cdnainIsBest?
      my $trd = ( $cstrandIsOpposite ) ? "trAnti": "trSense";
      $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
      my( $uprot, $ustart, $uend) = orfParts($bestorfi);
      my $ulen=  $bestorfi->{length};  
      my $ualen= length($uprot); $ualen-- if(substr($uprot,-1) eq '*');
      my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
      my $ucompl= $bestorfi->{complete};
      $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
      $addattr .=";" if($addattr); 
      $addattr .= "cdnaorf=$ualen,$upcds%,$ucompl,trbest:$cdnainIsBest,$trd"; #? ;aacdnaoff=$ustart-$uend;aacdnaprot=$uprot";  
    }
    
    if($cdnainIsBest) {
      # replace?, notice, may need to change CDS like intronfix : depends on diff
          
      my $aadif= length($orfproti) - length($orfprot); 
      my $gadif= $bestorfi->{goodlen} - $bestorf->{goodlen}; 
      # my $trdif= length($cdnain) - length($cdna); 
      my $trdif= $cdnainlen - $cdnagfflen; # or cdnaingood - cdnagffgoodlen ?
      $cdnalen= $cdnainlen; # change it
      
      my $trdif1=$trdif;
      map{ $_="+".$_ if($_>0); $_= (($_==0)?"eq:":"NE:").$_; } ($aadif,$gadif,$trdif);
      my $aac= $bestorfi->{complete} - $bestorf->{complete}; # complete == 3,2,1,0
      
      # also compare cdna vs cdnain : size, where diff (inside, ends?)
      # mrna-attr: Coverage=40.0;Identity=99.7 < use this to see diff?
      my $nx0= @exongff; # no way to count exons in cdnain here.
      my $trin=0; my $nxeq=0;
      
      ## change cdnain_isrev to cstrandIsOpposite here?
      # my $trd= ($cdnain_isrev) ? "trREV" : "tr";
      my $trd= ($cstrandIsOpposite) ? "trdAnti" : "trdSense"; # "trSens" ?
      
      ## FIXME4: strand reversal, $cdnain_isrev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
      #      .. or need to reverse strand of output gff
      #   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
      #   see fixstrand..
      
      ## fixme, use cstrandIsOpposite; also $prostart5, $proend3, are wrong or exongff needs strand change
      ## BUG, dont change this part strand when it is other chimera part that has bestorf w/ other strand..
      ## drop this for now, not right for cases seen.; use trdAnti note ; other flag for problems?      
#       if($cstrandIsOpposite and not $fixstrand) {
#         if($gstrand eq "-") { $gstrand="+"; } elsif($gstrand eq "+") { $gstrand="-"; }
#         map{ $_->[6]= $gstrand } @generec;
#         ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
#         @exongff= grep{ $_->[2] eq "exon" } @generec;
#         ## .. other problems w/ restrand..  oldCDS ? 
#       }
  

      if($cdna eq $cdnain) { $trd .= "eq"; $nxeq=$nx0; }
      elsif( $cdnaindex >=0 ) { # from index($cdnain|cdnainrev, $cdna)
        $trd .= "$trdif,in$trin";
        $nxeq= ($nx0==1)?1:$nx0-1;
        if($cstrandIsOpposite) {}
      }
      # elsif( ($trin=index($cdnain, $cdna)) >=0 ) { $trd .= "in$trdif,$trin"; $nxeq= ($nx0==1)?1:$nx0-1; }  #nxeq maybe; calc?
      # elsif( ($trin=index($cdnain, revcomp($cdna))) >=0 ) { $trd .="rc$trdif,$trin"; }
      else { 
        $trd .= $trdif; # * check each \@exongff > exseq index cdnain ? report if ends or inner x diff
        my($cdna2, $cdnalen2, $exonseq)= getcdna( \@exongff, 1);
        my ($ix,$le)=(0,1); 
        foreach my $xs (@$exonseq) {
          $ix++; my $xi= index($cdnain, $xs);  my $xr=""; 
          if($xi<0) { $xi= index($cdnain, revcomp($xs)); $xr="c" if($xi>=0); }
          if($xi>=0) { my $xe= length($xs)+$xi; my $xb=$xi+1; $trd.=",xeq$ix:$xr$xb-$xe"; $nxeq++; 
             if( $le>1 and (my $g= $xb - $le) > 0 ) { $trd.= "g$g"; } $le=1+$xe; }
          else { $trd.=",xne$ix"; }
          my $xgap= $xs =~ tr/N/N/;  $trd .="N$xgap" if($xgap>0);
        }
      }

      #check further if $nxeq > 0; xinner change?  count NNN gaps in cdna, cdnain; look for gaps at end if inner match
      if($nxeq > 0) {
        my $gaps = $cdna =~ tr/N/N/;
        my $gapi = $cdnain =~ tr/N/N/;
        if($trdif1>25) { 
          my($cdna2, $cdnalen2)= getcdna( \@exongff, 0, 100, _max(100,$trdif1));
          $gaps = $cdna2 =~ tr/N/N/;
        }
        if($gaps > $gapi) { $trd.= ",tN:$gaps"; }
      }       
      
      my $cov= ($mrna->[8] =~ m/(cov|Coverage)=(\d+)/) ? $2 : 0;
      $trd .= ",tcov:$cov" if($cov); # cov < 99 explains
      my $changenote= "cdnabest=aa$aadif,gd$gadif,dfull:$aac,$trd,nx0:$nx0";
      $addattr .=";" if($addattr); $addattr .= $changenote; $changed++;
      ($orfprot, $prostart5, $proend3, $utrorf)= ($orfproti, $prostart5i, $proend3i, $utrorfi);
    }
  }


  if($orfprot) {
    my($annew,$anold)= stripOldAnnot($mrna) if($REANNOTATE); # also changes mrna
    
    my ($cdsgff, $cdsattr)= getCDSgff( \@exongff, $orfprot, $prostart5, $proend3, $cdnalen);
    # test $prostart5, $proend3 vs old
    
    my $diff=0; 
    $diff=1 if($REANNOTATE or $cdnainIsBest);
    
    my $diffcancel=0; my $oldstat="";  
    if(ref $cdsgff and @$cdsgff > 0) {
      my @newCDS= sort _sortgene @$cdsgff; # do here.
      my @oldsave=();
      # push @generec, @newCDS; # wait till decide if new != old
      #? @generec= sort _sortgene @generec; # sorts by genome loc, not stranded

      ## reversed CDS look wrong at end exons (always both ends? for nc>2)    
      if(@oldCDS > 0) {
        my($oldal, $newal, $oldprostart5, $oldproend3)= (0,0,0,0);
        # compare, report in $addattr
        @oldCDS= sort _sortgene @oldCDS;
        ## my @newCDS= sort _sortgene @$cdsgff;
        my $newprot= ($cdsattr =~ m/protein=([^;\s]+)/) ? $1 : "";
        if (@oldCDS == @newCDS) { $oldstat="ocds=eqn"; } else {  $oldstat="ocds=NEn"; $diff++; }
        
        if($oldprot) { 
          (my $op=$oldprot) =~ s/\*$//; (my $np=$newprot) =~ s/\*$//; 
          $oldal= length($op);
          $newal= length($np);
          my $oldap= ($cdnalen>0) ? int(0.5 + 300 * $oldal / $cdnalen) : 0;

          my $eq=($op eq $np)?1:0; $diff++ unless($eq); 
          my $da= $newal - $oldal; $da="+$da" if($da>0);
          $oldstat .= ($eq) ? ",eqaa" : ",NEaa:$da"; 
          
          my $astat=0;
          if($oldprot =~ /^M/) { $astat |= 1; }
          if($oldprot =~ /\*$/) { $astat |= 2;} # ** AUG proteins lack '*' unless pasa-updated
          elsif($geneid =~ /AUG/) { $astat |= 2; } #  also check for inner X == augustus-fake for stop ?

          my $istop= index($oldprot,'*');
          ## find single 'X' inside prot, but allow this for NNN genome: SFXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQP
          if($geneid =~ /AUG/ and (my $ix=index($oldprot,'X')) >0 ) { 
            if(substr($oldprot,$ix,3) eq "XXX") { } # ignore 2+; really should look at genome dna to decide
            else { $istop= $ix if($istop < 1 or $ix<$istop); }
            }
          # if($istop < 1 and $geneid =~ /AUG/) { $istop= index($oldprot,'X'); }

          # my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
          if($istop > 0 and $istop < length($oldprot)-1) { $astat= "partialinner"; } # innerstop ?
          elsif($astat == 3) { $astat="complete"; }
          elsif($astat == 1) { $astat="partial3"; }
          elsif($astat == 2) { $astat="partial5"; }
          else { $astat="partial"; }
          $astat =~ s/^(.)/o/; # dont confuse w/ new in scans
          $oldstat .= ";oaaln=$oldal,$oldap%,$astat";
          
        } else { $diff++; }

   # FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
        my $partialAtEndOk= 0;
        for(my $i=0; $i<@oldCDS; $i++) { 
          my $oc= $oldCDS[$i];
          $partialAtEndOk=1 if($oc->[3] < 450); # dont know high end to check
          unless($i<@newCDS and $oc->[3] == $newCDS[$i]->[3] and $oc->[4] == $newCDS[$i]->[4]) { 
            $diff++; $oldstat .= ",NE$i";  ## BUG: clone oc
            if($debug>1) { my @doc= @$oc; $doc[2]="oldCDS"; $doc[0] =~ s/^/#/; push @oldsave, \@doc; }
            }
          }
          
       
        do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $partialAtEndOk and $oldstat =~ /artial/ and $newal < $oldal);

   # FIXME4: cancel changes that replace all CDS with completely new CDS, eg for transposon-spans
        $oldprostart5= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3]; 
        $oldproend3  = ($gstrand eq "-") ? $oldCDS[0]->[3] : $oldCDS[-1]->[4]; 
        $diff++ if($diff==0 and ($oldproend3 != $proend3 or $oldprostart5 !=  $prostart5));   # off by -1,-2 bugs   

        #  oldc:   [-----]
        #  newc:           [-----]  ; bad = shifted too much
        #  newc:        [-----]     ; bad "
        #  newc:     [---]          ; ok  = not shifted, same stop
        
        if($diff and $oldal>0) { # and $oldstat !~ /artialinner/
          ## .. bugs here ??
          #my($tb,$te)=($prostart5,$proend3); ($tb,$te)= ($te,$tb) if($tb>$te);
          #my($lb,$le)=($oldprostart5,$oldproend3); ($lb,$le)= ($le,$lb) if($lb>$le);
          my($tb,$te)=($newCDS[0]->[3],$newCDS[-1]->[4]); ($tb,$te)= ($te,$tb) if($tb>$te);
          my($lb,$le)=($oldCDS[0]->[3],$oldCDS[-1]->[4]); ($lb,$le)= ($le,$lb) if($lb>$le);
          if($lb==0 or $le==0 or $tb==0 or $te==0) { $tb=$te=$lb=$le=0; } # bug somewhere ...
          
          # argg; location stats not best here, introns have big effect.
          # .. this cancels most span increases...
          if(($tb < $lb and $te < $le) or ($tb > $lb and $te > $le)) { # CDS shifted along exons
            my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
            my $maxo= 1 + $be - $bb; # ** not abs
            my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
            my $pover= $maxo/$leno; # neg desired: cancel any case of maxo < 0 
            #??do{ $diff=0; $diffcancel=1; } if($pover < 0.33); # cancel change, too little cds similar
            do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $pover < 0); # cancel change, too little cds similar
          }
        }

        $cdsattr .=";$oldstat";
        $changed++ if($diff>0);
      } else {
        $changed++; $diff=1;
      }

      if($diff) {
        push @generec, @newCDS; 
        push @generec, @oldsave; 
      } else {
        push @generec, @oldCDS;  
      }
    }
    # else { $diff++; } ## is this error? cdnainIsBest = cdnain best?
    #? else { push @generec, @newCDS;  } # no @$cdsgff no oldCDS
  
    ## add to addattr: count utrs (and span?) and add utrx= if >2
    # change getCDSgff() to return utr-exons?
    # if($u5>2 or $u3>2) { $un=$u5+$u3; $g[0]=~s/$/;utrx=$u5,$u3/; } 

    ## fixme: diff==0, dont remove old attr, add new or not? 
    #     or make option? sometimes always replace? use ocds=eqn,eqaa as flag for no change
    # .. except for some changes above result in NEaa but are cancelled.
    
    if($utrorf) {
      my( $uprot, $ustart, $uend) = orfParts($utrorf);
      my $ulen=  $utrorf->{length};  
      my $ualen= length($uprot);  
      my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
      my $ucompl= $utrorf->{complete};
      $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
      $addattr .=";" if($addattr); 
      $addattr .= "aautrlen=$ualen,$upcds%,$ucompl;utroff=$ustart-$uend;utrprot=$uprot";  
    }
    
    if($cdsattr and $diff==0) {
      # check prot, add '*' if missing and complete
      # remove prot, new aalen if  $diffcancel
      unless( $mrna->[8] =~ m/;protein=/ or $diffcancel) {  ## should set diff=1 for this
        $diff=1; # $mrna->[8] =~ s/$/;$cdsattr/; 
        } ## add if missing
      if($diffcancel) { 
        my ($st)= $oldstat =~ m/ocds=([^;\s]+)/; 
        my ($nal)= $cdsattr =~ m/aalen=([^;\s]+)/;
        $mrna->[8] =~ s/$/;ocds=cancel:$st,nal:$nal/; 
        
      } elsif(not $nostopcodon and $cdsattr =~ /complete/ and $oldprot =~ /\w/ and $oldprot !~ /\*$/) {
        my $ops= $oldprot.'*';  # this WAS bug for protein=*Maaaa* 
        $mrna->[8] =~ s/;protein=$oldprot/;protein=$ops/; # no end ;
      }
      
    } # elsif()
    
    if($cdsattr and $diff>0) {  
        ## cdsattr keys: cxlen,aalen,protein,utrx
      my @keys= $cdsattr =~ m/(\w+)=/g;
      my $keyp= join('|',@keys);
      $mrna->[8] =~ s/[;]?($keyp)=([^;\n]+)//g;
      $mrna->[8] =~ s/$/;$cdsattr/;
    }
    
    $mrna->[8] =~ s/$/;$addattr/ if($addattr);

  } else { # no orfprot : can happen, but error if oldCDS << need to restore those here ??
    push @generec, @oldCDS if(@oldCDS);
    my $addattr="ocds=cancel:NEn,aatrans-failed";
    $mrna->[8] =~ s/$/;$addattr/;
  }

  if($ifix == 1) {  @generev= @generec; $changed1=$changed; }  # already cloned, now preserve
  else {  @genefwd= @generec;  $changed0=$changed; }
  
  } # BIG loop for fixstrand 

## somewhere add nostopcodon opt:
#  if($nostopcodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
  
  if($fixstrand) {
    my $mrna;
    ($mrna)  = grep{ $_->[2] eq "mRNA" } @genefwd;
    my ($aafwd)= ($mrna->[8] =~ m/aalen=([^;\s]+)/) ? $1 : "";
    my $protfwd= ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";
    ($mrna)  = grep{ $_->[2] eq "mRNA" } @generev;
    my ($aarev)= ($mrna->[8] =~ m/aalen=([^;\s]+)/) ? $1 : "";
    my $protrev= ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";

    my $dlenrev= length($protrev) - length($protfwd);
    my $plenrev= abs($dlenrev) / _max(1, _max(length($protrev) , length($protfwd)));
    my $acfwd= ($aafwd =~ /complete/)?1:0; 
    my $acrev= ($aarev =~ /complete/)?1:0; 
    my $revbest=0;
    if($acfwd ne $acrev and $plenrev < 0.33)  {
      $revbest= ($acrev) ? 1 : 0;
    } else {
      $revbest=($dlenrev > 0)?1:0;
    }
    if($revbest) {  @generec= @generev; $changed=$changed1; } 
    else { @generec= @genefwd; $changed=$changed0; }
  }
  
  putgene( \@generec, $geneother, "");
  return ($changed) ? 1 : 0;
}


sub get_dna {
  my($fasta, $ref, $start, $stop)= @_; #, $fasta_db_ref
  unless( $ref && $stop>0) { warn "need ref-ID, start, stop location\n"; return; }
 
  require Bio::DB::Fasta;  # FIXME: not portable w/o parts of BioPerl !
  ## need also (Bio::DB::SeqI Bio::Root::Root)
  ## (DB_File GDBM_File NDBM_File SDBM_File) : are these core Perl modules?
  my $havedb= ref $fasta_db;
  unless($havedb) {
    my $db = eval { Bio::DB::Fasta->new($fasta); }
      or die "$@\nCan't open sequence file(s). "; # and return;
    $fasta_db= $db;  
    }
  
  my $seq = $fasta_db->seq($ref, $start => $stop) 
      or return; ## die "cant locate seq of $ref:$start-$stop\n";#? and return;
  $seq= $seq->seq if(ref $seq); # is this weird bioperl change here or not
  #?? Die if no seq? likely bad genome fasta lots of errs
  return $seq;
}


## revise longest_orf_finder

# in cdna_proteins.pm;  revised for both strands
sub getAllOrfs {   
  my ($input_sequence, $strands) = @_;   # ,$flags
  $strands ||= "fwd";
  # $flags ||="";
  
  return undef unless ($input_sequence or length ($input_sequence) >= 3) ;
  $input_sequence = uc ($input_sequence); # was lc() change all to uc()
  
#  ## fixme: this screws seq indices: start,stop; instead restrict @starts,@stops
#   if($flags =~ /dropnnn|chomp/i) {
#     $input_sequence =~ s/^N+//;
#     $input_sequence =~ s/N+$//; 
#   }
  
  my (@starts, @stops, @orfs);

  if (1) {  # forward_strand_only();
  @stops  = identify_putative_stops($input_sequence);
  @starts = identify_putative_starts($input_sequence,\@stops);
  @orfs   = get_orfs (\@starts, \@stops, $input_sequence, '+');
  }

  if (@orfs) {  # change here, get both longest complete, longest partial
    ## set in order of decreasing length
    my $keylen=($USEGOODLEN)?"goodlen":"length";
    @orfs = sort {$b->{$keylen} <=> $a->{$keylen} or $b->{complete} <=> $a->{complete}} @orfs;
    
   # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
    my $longest = $orfs[0];    
    my($longfull) = grep { $_->{complete} == 3 } @orfs;
    return ($longest, $longfull, \@orfs);
  } else {
    return undef;
  }
}

# in cdna_proteins.pm;  revised 
sub orfParts
{
  my($longorf) = @_;
  if(ref $longorf) { 
    return($longorf->{protein}, $longorf->{start}, $longorf->{stop});
  } else {
    return("",0,0);
  }
}


# in cdna_proteins.pm;  revised 
sub get_orfs {
  # my $self= shift;
  my( $starts_ref, $stops_ref, $seq, $direction) = @_;
    
  # unless ($starts_ref && $stops_ref && $seq && $direction) { warn "get_orfs: params not appropriate"; }  
  # want only max orf, complete + partial
  # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
   
	my %last_stop_part = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
	my %last_stop_full = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
  my @orfs;
  my $seq_length = length ($seq);
  my $norf=0;
  
  ## $seq .= "###" if($debug); # why're we getting bad prot at partial end? isnt prot but cds range needs adjust for -2,-1
  
  foreach my $start_pos (@{$starts_ref}) {
		my $start_pos_frame = $start_pos % 3;
		# includes partial5: starts at 0,1,2 unless real start
		my $isfull5= ($start_pos+3 <= $seq_length and substr($seq, $start_pos, 3) eq "ATG")?1:0;
		
		foreach my $stop_pos (@{$stops_ref}) {
		  if ( ($stop_pos > $start_pos) && #end3 > end5
				 ( ($stop_pos - $start_pos) % 3 == 0) && #must be in-frame
				   # ($start_pos > $last_delete_pos{$start_pos_frame}) # dgg: bad for partial5 + full
				 ($isfull5 ? $start_pos > $last_stop_full{$start_pos_frame} : $start_pos > $last_stop_part{$start_pos_frame} )
				 ) #only count each stop once.
			{
		    # includes partial3: stops at end- 0,1,2 unless real stop
				
				#dgg.no# $last_delete_pos{$start_pos_frame} = $stop_pos;
				if($isfull5) { $last_stop_full{$start_pos_frame} = $stop_pos; }
			  else { $last_stop_part{$start_pos_frame} = $stop_pos; }
			
			  my $stopplus3= _min( $stop_pos+3, $seq_length); #dgg patch; getting overruns.
				my $orflen = ($stopplus3 - $start_pos);
				my $stopoff = $orflen % 3;
				if($stopoff>0) {  $orflen -= $stopoff; $stopplus3 -= $stopoff; } # do before pull seq..
			  
				my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), $stopplus3);
				my ($start, $stop) = ($direction eq '+') ? ($start_pos_adj, $stop_pos_adj) 
					: (&revcomp_coord($start_pos_adj, $seq_length), &revcomp_coord($stop_pos_adj, $seq_length));
				
				my $orfSeq =  substr ($seq, $start_pos, $orflen); #include the stop codon too.
				my $protein= translate_sequence($orfSeq, 1); ## from Nuc_translator

				if ($protein =~ /\*.*\*/) {
				  print STDERR "#ERR: innerstop\n"; #? last unless($debug);
          #?? last; ## error
				}
				
			  $orflen= length($orfSeq);
				my $aalen = length($protein);
				# $orfoff = $orflen % 3;
				# if($orfoff>0) { $stop -= $orfoff; $orflen -= $orfoff; $orfSeq}  # $orflen > 3*$aalen ; ($orflen % 3 != 0) .. need to adjust stop
				  
				my $isfull= 0;
				$isfull |= 1 if($protein =~ /^M/);
				$isfull |= 2 if($protein =~ /\*$/);
				
        ## flag problem if XXX count > non-X count
        ## OPTION? remove leading/trailing XXX and reset sizes
        my $nxxx= $protein =~ tr/X/X/;
        ##my $isbad= ($nxxx > 0.5*length($protein)) ? 1: 0;
        my $goodlen= $orflen - 3*$nxxx;
        
        $norf++;
        # opt to drop *stopcodon, here or where, 
        
        if($debug>1){ ##my $aalen=length($protein);
          print STDERR "#\n" if($norf==1); # ,glen=$goodlen
          #no,wrong# my $dc=  ($direction eq '+') ? 'fwd' : 'rev';
          print STDERR "#DBG: orf$norf,olen=$orflen,alen=$aalen,at=$start-$stop,aa1=",
            substr($protein,0,9),"..",substr($protein,-3,3),"\n";
        }
				my $orf = { sequence => $orfSeq,
							protein => $protein,
							complete => $isfull,  
							start=>$start,
							stop=>$stop,
							length=>$orflen,
							goodlen=>$goodlen, # use to avoid XXXXXXX* crap
							orient=>$direction,
							};
				push (@orfs, $orf);
				last;   
			}
		}
  }
  return (@orfs);
}



# in cdna_proteins.pm;   
sub identify_putative_starts {
  my ( $seq, $stops_aref) = @_;
  my %starts;
  my %stops;
  foreach my $stop (@$stops_aref) {
		$stops{$stop} = 1;
    }
	
  if(1) {  # ($self->{ALLOW_5PRIME_PARTIALS} || $self->{ALLOW_NON_MET_STARTS}) 
    my $i=0;
    if($seq =~ /^N/) {  # FIXME: skip leading NNN
      my $n= length($seq);
      $i++ while( $i<$n && substr($seq,$i,1) eq 'N');
    }
		$starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
    }
    
  if (1) {  # ! $self->{ALLOW_NON_MET_STARTS}  #Look for ATG start codons.
		my $start_pos = index ($seq, "ATG");
		while ($start_pos != -1) {
			$starts{$start_pos} = 1;
			$start_pos = index ($seq, "ATG", ($start_pos + 1));
		  }
    } 

  my @starts = sort {$a<=>$b} keys %starts;
  return (@starts);
}


# in cdna_proteins.pm;   
sub identify_putative_stops {
  my ($seq) = @_;
  my %stops;
  
  if(1){  # $self->{ALLOW_3PRIME_PARTIALS}
	## count terminal 3 nts as possible ORF terminators.
	my $e = length ($seq);
  if($seq =~ /N$/) { # FIXME: skip trailing NNN
    $e-- while( $e>1 && substr($seq,$e-1,1) eq 'N');
  }
	$stops{$e} = 1;
	$e--; $stops{$e} = 1;
	$e--; $stops{$e} = 1;
  }
  
  # my @stop_codons = @{$self->{stop_codons}};
  # my @stop_codons = &Nuc_translator::get_stop_codons(); # live call, depends on current genetic code.
  my @stop_codons = qw(TAA TAG TGA);
  
  foreach my $stop_codon (@stop_codons) {
    my $stop_pos = index ($seq, $stop_codon);
    while ($stop_pos != -1) {
      $stops{$stop_pos} = 1;
      $stop_pos = index ($seq, $stop_codon, ($stop_pos + 1)); #include the stop codon too.
      }
  }
  my @stops = sort {$a<=>$b} keys %stops;
  return (@stops);
}


# in cdna_proteins.pm;   
sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


# in cdna_proteins.pm;   
sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}

# in cdna_proteins.pm;   
# parts from PASA/PasaLib/Nuc_translater.pm 
use vars qw ($currentCode %codon_table);


# in cdna_proteins.pm;   
sub translate_sequence {
  my ($sequence, $frame) = @_;
    
  $sequence = uc ($sequence); # redundant now
	$sequence =~ tr/U/T/;
  my $seq_length = length ($sequence);
  unless ($frame >= 1 and $frame <= 6) { 
		warn "Frame $frame is not allowed. Only between 1 and 6"; # die?
		return -1; # 
	}
	
	if ($frame > 3) {
		# on reverse strand. Revcomp the sequence and reset the frame
		$sequence = revcomp($sequence);
		if ($frame == 4) {
			$frame = 1;
		}
		elsif ($frame == 5) {
			$frame = 2;
		}
		elsif ($frame == 6) {
			$frame = 3;
		}
	}
	
  # $sequence =~ tr/T/U/; # dont need this; change codon_table
  my $start_point = $frame - 1;
  my $protein_sequence;
  for (my $i = $start_point; $i < $seq_length; $i+=3) {
      my $codon = substr($sequence, $i, 3); # problem here for i>seq_length-3 ?? or in caller getting +1,+2 > true len
      my $amino_acid;
      if (exists($codon_table{$codon})) {
          $amino_acid = $codon_table{$codon};
      } else {
          if (length($codon) == 3) {
              $amino_acid = 'X';
          } else {
              $amino_acid = "";
          }
      }
      $protein_sequence .= $amino_acid;
  }
  return($protein_sequence);
}


# in cdna_proteins.pm;   
BEGIN {
  $currentCode = "universal";
  # stops: TAG  TGA  TAA
  ## add/allow N in silent positions; need other codon table w/ extended AA vals for ambiguous
  %codon_table = (    
  TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L',
  CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',  CTN => 'L',
  ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M',
  GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V',  GTN => 'V',
  TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S',  TCN => 'S',
  CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P',  CCN => 'P',
  ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',  ACN => 'T',
  GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',  GCN => 'A',
  TAT => 'Y', TAC => 'Y', TAA => '*', TAG => '*',
  CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q',
  AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K',
  GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
  TGT => 'C', TGC => 'C', TGA => '*', TGG => 'W',
  CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',  CGN => 'R',
  AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R',
  GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',  GGN => 'G', 
  );
}

1;
