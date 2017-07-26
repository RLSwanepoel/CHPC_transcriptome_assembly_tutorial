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

=item FIXMEs

  2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
  	patch in sumblastpart: $or add to bspans
  	
	## FIXME asmrna_dupfilter.pl: check for ID mismatches == too many zeros in aa.count /tr.count / blastn
	## eg. tables have diff prefixes like 'litova:vel' vs 'litovavel'

=item revise classifier/output

  FIXME: make OUTSPANTAB=1 default, see below alnclass.sh as better/faster classifier.
  
  update: fold in tested methods for tr class assignment, to replace outclusters
  -- newer class assignment: using %ident (>98) and %align (>50%)
          from  aabugs4qual/tsaevg/cdsidentclass.sh

  shorter alignments of identity are tested as valid criteria for alternate-transcripts
  using genome-mapped transcripts, hi-ident + part-align gives *mostly* alt-tr, some paralogs,
  depending on species & freq of hi-ident paralogs.
              
  update: add refblast  input table for added classification, per
            aabugs4/aabugs4qual/tsaevg/classparts.pl
        tall4 format: Refid Trid  Bits Iden Algn Rlen Tlen, where Ref/Tr may be swapped columns

  update: add 2ndary alt-tr table, in aa.qual format?,
      from alt-tr called by trasm soft, but excluded by cdhit/other filtering as not distinct
      these all should have same geneid + tNNN suffix as primary input tr (defined in -aasize -blast trself.blast)
      ID matching only used, plus aa.qual, to decide whether to keep/drop
      - input to classifier may be 2nd OUTSPANTAB, generated from 1st pass of this and 2nd-alts.aa.qual table.

=item add UTR BAD/POOR filters

  problem now is that main class accumulates utrbad/poor genes, many w/ other utrbad alternates,
  so as more trasm sets are added to this filtering, main class (and some of others) increase w/ junk.
  
  use input aasizes/aaqual scores
  when reading cdhit/blast aligns (esp cdhit clusters), where feasible swap/replace UTRBAD top/first gene 
    if equivalent utrok gene is in cluster/align (ie. utrok/utrbad have same prot/cds size, or ok is nearly same),
    esp. if utrok is also aacomplete.    
  per evigene/scripts/prot/aabest3.sh



=item add cd-hit-est input alignments

  maybe replace cd-hit-aa with cd-hit-cds(est) for both aa-based and nt-based equivalence/reduction
  of trasm sets.  cd-hit-aa has problem of lumping paralogs w/ silent subsititutions, want to keep
  those in early filtering.  cd-hit-cds/est gives approx same classes of same-locus vs diff-locus as
  blastn or blat (depending on parameters).
    cd-hit-est -c 0.99 -G 0 -aS $pALIGN -l 150 -d 0 -i cacao11pub3ig.cds
        pIDENT=0.99 is good as w/ others; pALIGN=0.25 .. 0.50
        
  cacao11pub3ig.class3
  # alternate classifier, replace both cd-hit aa and mblast-tr ? 
  cacao11g: 44403 cds, 29451 loci, 7630 loci have alts,
  blastn (99% id, >=50%? align): alt class=16179 alts, 484 diff locus; diffclass: noclass=noalts; altmid=236 alts, 1156 diff
      false-alts: 484; false-loci: 236 of 44403
  cdhit25: 29893 clusters, 730 loci have alts in diff clusters; 295 diff loci in same cluster;
      false-alts: 295; false-loci: 730 of 44403
        -- about same as above self.mblast (same loci?)
  cdhit50: 30057 clusters, 841 loci alts diff cluster; 252 diff loci in same cluster
    -- cd-hit pALIGN=0.33 may be good choice;
     
=item updated for blat psl input

  blat -out=psl -t=dna -minIdentity=98 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat
  blat -out=psl -t=dna  -minScore=99  -minIdentity=99 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat99
  #  -minScore=99 or such to reduce fragment aligns, 30 default, nmatch - nmiss
  : tests show blat and megablast give nearly same results; blat can be very slow vs mblast for large tr set
  
  
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

my $debug= 0;
use constant VERSION => '2013.09.01'; # 08.09'; #08.07; 03.25; 24; # '2013.02.21';

my $OUTSPANTAB= 1;  # make default until replace outclusters()
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

# classifier options
our($AAMIN,$AAPART,$AAMINBAD, $AAMINPOO, $AADUP_IDENT, $OK_CDSUTR, $BAD_CDSUTR, $BAD_GAPS, $TINYALN,
    $IS_CDSALIGN,$ALTFRAG,$PHIALN,$NHIALN,$PHI,$PMID,$PLOW);

$AAMIN =$ENV{aamin}||40;   # for aacomplete, utrok
$AAPART=$ENV{aapart}||100; # for aapartial
$AAMINBAD=200;  # for utrbad class
$AAMINPOO=100;  # for utrpoor class
$BAD_GAPS= 25;  # % gaps in AA
## add UTR BAD/POOR filters, per evigene/scripts/prot/aabest3.sh
# $pAADUP=98
$AADUP_IDENT=98; # option for aacluster ident drops

my $MIN_EQGENE=33; # global/opt..

## $TINYALN=$ENV{mina}||50; 
$TINYALN = 25; # ignore less than this == MINAL
$OK_CDSUTR= 60;
$BAD_CDSUTR= 30; # % CDS/trlen
$ALTFRAG= 0.5;
$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN
#     my $tinyalndiff= ($aln - $ww < $NHIALN) ? 1 : 0 # TEST3 add

$PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80;  # tr-self %identity to classify alt-tr
## <90% ident prob too low to class alt-tr; use pmid=95 plow=90 ??

$IS_CDSALIGN=$ENV{cdsw}||0; # WHICH default?  tralign tests better than cds, but cds for other, eg. cds-cdhits

# my $MINALIGN= 90; # NOT USED NOW; use just TINYALN; change to several levels
# my $MINSIZE = 999; # not used yet

# my $IDPREFIX= $ENV{idprefix} || ""; # refgene id pattern; need for aaclusters at least,    
my ($aasizes,$trsizes,$blatpsl,$blastab,$lastz,$bcdhit,$aablast,$aanames,$aacdhit,$outeqtab,$outclass,
    $dupids,$logfile,$head,$eqgene)= (0) x 20;

my $optok= GetOptions(
  "aasizes=s", \$aasizes, 
  "trsizes=s", \$trsizes, 
  "ablastab=s", \$aablast,   # this is traa-refaa.blastp
  "anames=s", \$aanames,   # variant aablast used also for naming
  "acdhit=s", \$aacdhit,    # this is traa-self.cdhit.clstr

  "blastab=s", \$blastab,    # this is tr-self.blastn, -CDSALIGN for cds-self.blastn
  "lastz=s", \$lastz, # lastz general format; was -blastz option -blast[ab] conflict **
  "blat=s", \$blatpsl, # tr-self.blat;  other input format opts here..?  -informat=xxx
  "bcdhit=s", \$bcdhit, #  tr-self.cdhits.clstr; 
  # ntalign=s, ntformat=blastab|blatpsl|cdhit ..
  "dupids=s", \$dupids, #  

  "eqmap|eqgene=s", \$eqgene,    # 130901: mapping equivalence table, of $evigene/equalgene.pl 
  
  "aligntab|outeqtab=s", \$outeqtab,  ## this should have aliases: -outalntab and -inalntab, or just -aligntab ?  
  "outclass=s", \$outclass,   # 2nd out: outclasstab
  "TINYALN|MINALIGN=i", \$TINYALN,  # pTINYALN ?
#  "MINSIZE=i", \$MINSIZE,  
  "pCDSOK=i", \$OK_CDSUTR,   # CDSOKUTR
  "pCDSBAD=i", \$BAD_CDSUTR, # CDSBADUTR 
  "ALTFRAG|fragment=s", \$ALTFRAG,  # pALTFRAG ?
  "OUTSPANTAB!", \$OUTSPANTAB, 
  "CDSALIGN!", \$IS_CDSALIGN, 
  "debug!", \$debug, 
#  "logfile=s", \$logfile,  ## add??
#  "geneidprefix=s", \$IDPREFIX,  
  );


warn "# EvidentialGene asmrna_dupfilter.pl VERSION ",VERSION,"\n" if($debug); # change to loggit() ?
my $hasbalign= ($blastab or $lastz or $blatpsl or $bcdhit) ? 1 : 0;

die "usage:  asmrna_dupfilter.pl -aasize=name.aa.count -trsize=name.tr.count 
 input rna-align: -blast=name.blastn | -blat=name.blatpsl | -bcdhit=name.cdhit.clstr
 opts: -CDSUTR=$OK_CDSUTR percents  -ablastab=traa-refaa-blast.table  .. more options, see source.
" unless($optok and $aasizes and ($hasbalign or $outeqtab)); ## and $trsizes < dont need if aasizes=aa.qual
#  -MINALIGN=$MINALIGN  -OUTSPANTAB

if($TINYALN>0) {
  $TINYALN= 100*$TINYALN if($TINYALN<1);
  ### $TINYALN=$MINALIGN if($TINYALN>$MINALIGN); # drop MINALIGN for TINYALN
}
$ALTFRAG= $ALTFRAG/100 if($ALTFRAG>1); # prop not percent

my(%aasize, %trsize, %aaqual, %aablast, %aablastref);
my(%aacluster,%aaclustermain); # globals for readAAcdhit
my(%better, %outrows, %validids, %bspans); ## %validids was %outids; 
my(%dupids, %dupfirst, $ndupdrop); $ndupdrop=0;
use constant DUPFILTER1 => 1; # tests ok
use constant DUPFILTER2 => 0;

my $OUTH= *STDOUT;

readSizes();
readDupIds($dupids)   if($dupids);
readAAblast($aablast) if($aablast);  # one or other of aablast,aanames
readAAnametab($aanames) if($aanames);
readAAcdhit($aacdhit) if($aacdhit); # also correctAAcluster()

my($eqgenes,$nineqgene,$neqgene)= ($eqgene) ? readEqualGene($eqgene) : (undef,0,0);
## $eqgenes->{$tid}{$sid} >= $MIN_EQGENE alignment

my $nbalign=0;
if($hasbalign) {
  if($outeqtab) {
    rename($outeqtab,"$outeqtab.old") if( -f $outeqtab );
    open(OUTH,">$outeqtab") or die "write $outeqtab";
    $OUTH= *OUTH;
  }
  
  # allow input of OUTSPANTAB instead of regenerate?
  if($blastab) { ($nbalign)= readblasttab($blastab); }
  elsif($lastz) { ($nbalign)= readlastz($lastz); }
  elsif($bcdhit) { ($nbalign)= readcdhit($bcdhit); }
  elsif($blatpsl) { ($nbalign)= readblatpsl($blatpsl); } # detect from input table ??
  if($outeqtab) { close($OUTH); $OUTH=undef; } # *STDOUT ?
  warn "# readalign: nids=$nbalign to $outeqtab\n" if $debug;
}

# unless($OUTSPANTAB) 
if($outeqtab) {
  my $infile=$outeqtab;   # infile == outeqtab ? STDIN?
  my $insorted=0;
  unless($infile and -f $infile) { 
    die "ERR: missing input align table -aligntab $outeqtab";
  } elsif($infile =~ /\.gz/) { 
    die "ERR: cant use gzipped align table -aligntab $outeqtab";
  }

  # >> set this BEFORE correctAAcluster from read{blast|cd|blat} : $validids{$id}
  if($aacdhit) {
    my $havevalid= scalar(%validids)?1:0;
    $havevalid= readIdsFromAlnTab($infile) unless( $havevalid); 
    correctAAcluster($havevalid); # update %aacluster,%aaclustermain
  }

  $OUTH= *STDOUT; 
  if($outclass) {
    rename($outclass,"$outclass.old") if( -f $outclass );
    open(OUTC,">$outclass") or die "write $outclass";
    $OUTH= *OUTC;
  }
  
  identityclass($OUTH,$infile,$insorted); 
  if($outclass) { close($OUTH); $OUTH=undef; } # *STDOUT ?
}

# if(0) { outclusters($OUTH) unless($OUTSPANTAB); } # DROP this old version

#.................................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub openRead { # in cdna_evigenesub.pm
  my($fna, $nostdin)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") 
  	 : ($fna =~ /stdin|^-/) ? *STDIN : open($hin,$fna);  
  # loggit(1,"ERR: openRead $fna") unless($ok);
	die "ERROR: openRead $fna" unless($ok);
  return ($ok,$hin);
}



sub readSizes {

  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  if($aasizes) {  
    ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
    ## drop faCount? require aa.qual here?

    # ($ok,$inh)= openRead($aasizes,1);
    open(F,$aasizes) or die "FAIL: read $aasizes ..."; 
    # if($aasizes =~ /count|qual/) { 
    # } else { open(F,"faCount $aasizes |") or die "FAIL: faCount $aasizes ..."; } # drop this..

    while(<F>) { next if(/^#/ or /^total/); 
      # my($id,$alen,@ax)=split; 
      my($id,$alen,$gap,$aqual,$tlen)=split; # aa.qual cols; gap is removed from alen
      unless($alen =~ /^\d/) { $nerr++; next; }
      if($aqual) { 
        $aqual .= "-gapbad" if($gap>0 and (100*$gap/($alen+$gap) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
        $aaqual{$id}= $aqual; 
        if($tlen =~ /^\d/) { $trsize{$id}= $tlen; $ntr++; }
        }
      $aasize{$id}=$alen; $naa++; 
      } close(F); 
  }
  
  if($trsizes) {
  if($trsizes =~ /^aaqual/ or $ntr>0) { 
    # got above;now default
  } elsif($trsizes =~ /^aasize|^cdssize/) { # for tr == cds
    foreach my $id (keys %aasize) { $trsize{$id}= 3*$aasize{$id}; }
  } else {  
    ## drop faCount? expected use is aa.qual w/ trsize column
    # ($ok,$inh)= openRead($trsizes,1);
    open(F,$trsizes) or die "FAIL: read $trsizes ..."; 
    # if($trsizes =~ /count|qual/) { 
    # } else { open(F,"faCount $trsizes |") or die "FAIL: faCount  $trsizes ..."; }
    while(<F>) { next if(/^#/ or /^total/); my($id,$al)=split; $trsize{$id}=$al; $ntr++; } close(F); 
  }
  }
 
 warn "# readSizes: naa=$naa; ntr=$ntr\n" if $debug;
 return($naa,$ntr); 
}

# fix from fastanrdb all.cds > all_nr.cds >> hdr has cds-identical ids; should have prefiltered this
sub readDupIds {
  my($infile)= @_;
  my($nids,$ndups,$inh)=(0,undef);
  %dupids= %dupfirst=();
	# ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    s/^>//; # if from fastanrdb
    my  @dupids= split; # grep /\w/ or any other?
    next unless(@dupids>1); #? or record all ids?
    my $firstid= $dupids[0];
    if($dupids{$firstid}) {
      my $nextdup= $firstid;
      $firstid= $dupids{$nextdup};
      shift @dupids;
      push @{$dupfirst{$firstid}}, @dupids;  
    } else {
      $dupfirst{$firstid}= \@dupids;  $nids++; 
    }
    map{ $dupids{$_}= $firstid } @dupids;
    $ndups += @dupids;
  } close($inh);  
  warn "# readDupIds: nfirst=$nids; ndups=$ndups\n" if $debug;
  return ($nids,$ndups);
}

=item readEqualGene

Table of map equivalences, add as adjunct to align tab of blast/lastz, which have omissionn mistakes ~10%?
	$evigene/scripts/equalgene.pl -in kf2mixx.main.gff -over kf2mixx.main.gff > kf2mixx.main.eqgene
* need to change main/over ids to input oids somewhere ..
Mainid                  Oid                     Overlapids
Funhe2Exy3m110884t1     Fungr1EG3m041003t1      Funhe2Exy3m129932t1/C97.72,Funhe2Exy3m110083t1/82.77    10098sc:638046-638339:.
Funhe2Exy3m125954t1     Fungr1EG3m043564t1      na      479sc:197743-197931:.
Funhe2Exy3m091394t1     Fungr1EG3m037315t1      Funhe2Exy3m076422t1/I100,Funhe2Exy3m083869t1/89.90,Funhe2Exy3m044035t1_C1/87.88,Funhe2Exy3m063866t1/87.87,Funhe2Exy3m054191t1_C1/87.83,Funhe2Exy3m059312t1/87.81,Funhe2Exy3m025662t1/87.79,Funhe2Exy3m060328t1_C1/86.87,Funhe2Exy3m050942t1/86.86,Funhe2Exy3m052474t1/86.86,Funhe2Exy3m085750t1/86.86,Funhe2Exy3m050513t1/86.84,Funhe2Exy3m040128t1/86.82,Funhe2Exy3m052598t1/86.78,Funhe2Exy3m052599t1/86.78,Funhe2Exy3m049962t1/85.85,Funhe2Exy3m054654t1/85.85,Funhe2Exy3m123553t1/57.83,Funhe2Exy3m026282t1_C1/53.87,Funhe2Exy3m090395t1_C2/0.97    1967sc:10172-10679:-

=cut

sub readEqualGene {
  my($infile)= @_;
  my($nids,$nov,$inh)=(0,0,undef);
  my %eqgene=();
	# ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
  	next unless(/^\w/);
		my($mid,$oid,$overids,$loc)=split"\t";
		next if($overids eq "na" or not $overids);
    $nids++; my $jov=0;
		my @ovd=split",",$overids;
		foreach my $ov (@ovd) {
			my($ovd,$cx)=split"/", $ov; 
			$cx.=".100" if($cx =~ /^I100/); 
			$cx=~s/^[IC]//; 
			my($ca,$xa)=split /[\.]/,$cx; $ca||=0; $xa||=0;
			# Argument "58^I1802sc:10833-11091:+\n" isn't numeric : split $ov not $_
			# uninitialized value in numeric gt ..
			if($xa > $ca) { $ca=$xa if($ca<$MIN_EQGENE and $xa>25+$MIN_EQGENE); } # what? adjust ca to xa? may be wrong. but CDS mapping can be wrong.
			if($ca >= $MIN_EQGENE) { $eqgene{$mid}{$ovd}= $eqgene{$ovd}{$mid}=$ca; $jov++;}
			## ca == cds-align equivalent from gmap overlap
		}
		$nov++ if($jov);
	}
 warn "# readEqualGene: nequal=$nov/$nids\n" if $debug;
	return (\%eqgene,$nids,$nov);
}


=item readAAblast

  maybe revise this to also use trasm.names table, computed for publicset naming,
  has essentially same info w/ best ref per tr.
  names format:
  TrID   Name   Align_score  RefID  RepID (uniprot)
  Funhe2Exx4m000455t12    CDD: Na_trans_assoc, Sodium ion transport-associated    100%,245/230,1997       CDD:219069      pfam06512
  Funhe2Exx4m000455t12    Sodium channel protein type 5 subunit alpha     100%,2057/2016,1997     RefID:UniRef50_Q14524   RepID:SCN5A_HUMAN
  Align_score = val%,nalign/nref,ntr
  
=cut 

sub readAAnametab
{
  my($aanametab)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  open(F, $aanametab) or die "FAIL: read $aanametab";
  while(<F>) { next unless(/^\w/); 
    my($td,$name,$alnscore,$rd,@more)=split"\t"; 
    $naabl++; 
    # alnscore format: 72%,3270/4555,3282 ; fail if mismatch?
    my($apct,$aln,$refsize,$tsize)= $alnscore =~ m=(^\d+)%,(\d+)/(\d+),(\d+)=;
    next unless($aln>0);
    $rd =~ s/^RefID://; # or CDD: or other
    my $bscore= $aln;
    # local $^W = 0;   # no warnings; no help
    unless($aablast{$td} and $bscore{$td} > $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    }
    
  } close(F); 
 warn "# readAAnametab: naabl=$naabl\n" if $debug;
 return($naabl); 
}


sub readAAblast 
{  
  my($aablast)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  if($aablast =~ /\.names$/) { return readAAnametab($aablast); }
  
  ## precheck format before sort of possibly very large file ..
  use constant nCHECK => 49;
	# ($ok,$inh)= openRead($aablast,1);
  open(F,$aablast) or  die "FAIL: read $aablast ..."; 
  while(<F>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    } elsif($swapids==0) {
      if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
      else {
      if($aasize{$td}) { $swapids= -1; }
      elsif($aasize{$rd}) { $swapids= 1; }
      else { $nblerr++; }
      }
    }
  last if($naabl > nCHECK);
  } close(F);
  if($nblerr>2 or $swapids==0) { 
    die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
      ." $nblerr trids from aasize not found in -aablast=$aablast\n";
  }
  
  $naabl=$nblerr= 0;  
  #?? add -aablastsorted option? sort wastes time if not needed.
  #old# open(F,"sort -k3,3nr $aablast |") or die "FAIL: read $aablast ..."; 
	# ($ok,$inh)= openRead($aablast,1);
  open(F,"sort -k3,3nr -k5,5nr -k1,1 -k2,2 $aablast |") or die "FAIL: read $aablast ..."; 
  
  while(<F>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    }
#    # done already ..
#     } elsif($swapids==0) {
#       if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
#       else {
#       if($aasize{$td}) { $swapids= -1; }
#       elsif($aasize{$rd}) { $swapids= 1; }
#       else { $nblerr++; }
#       }
#       if($nblerr>0 and $naabl > 19) { 
#         die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
#           ." trids from aasize not found in -aablast=$aablast\n";
#       }
#     }
    
    my $bscore= $v[0]; # bits, ident, algn
    # local $^W = 0;   # no warnings; no help
    ## maybe fixme: missing some uniq blastref{rd} here? pull aablastref out of aablast{} ?
    unless($aablast{$td} and $bscore{$td} > $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    
    }
    
  } close(F); 
 warn "# readAAblast: naabl=$naabl\n" if $debug;
 return($naabl); 
}


use constant SPANSUM => 1;

sub putspans {
  my($lq)= @_;
  # fixme: save no-match lq ids:  lq, aq, wq, na.... or self-score ?
  my $nmatch=0;
  # fixme: sort output bspans by tidn/taln? $bspans{$b}->[0]->[4] = xbit
  # fixme? throw away dupids here? should be in same bspans align cluster
  
  # my(%dupids, %dupfirst);
  my $dupfirst=""; 
if(DUPFILTER2) {    
  if($dupids) { $dupfirst= $dupids{$lq} || ""; }
}
  
  foreach my $lt (sort keys %bspans) {
    next if($lt eq $lq); # is lt eq lq allowed here?
    
    ## move this dupid filter before into readblastab, readcdhit ?
if(DUPFILTER2) {    
    if( $dupids and $dupids{$lt} ) {
      if($dupfirst eq $dupids{$lt}) { $ndupdrop++; next; }
      $dupfirst= $dupids{$lt};
    }
}
    
    my @sp= @{$bspans{$lt}};
    my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt,$torient)= (0) x 9;
    foreach my $sp (@sp) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or)= @$sp; # 2013.aug: IS_CDSALIGN add $or
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident; 
      ##$torient+=$or; ## weight by aln so tiny -or dont throw it off?
      $torient += $aln * $or; ## weight by aln so tiny -or dont throw it off?
      }
    # my $mis= $taln - $tidn; # dont need this; replace w/ tidn
    if($OUTSPANTAB) { 
      # add Qsize,Tsize  cds/trlen ?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $at= $aasize{$lt}||0; $at *=3;
      $wq= $trsize{$lq}||0; 
      $wt= $trsize{$lt}||0;  
      ## NO: add here?? $aaqual{$lq} ; $aaqual{$lt}; 
      ## FIXME taln may be rel CDS-len (aq,at) or TR-len (wq,wt)
      my $alnmax= _max(1, ($IS_CDSALIGN) ? _min($aq,$at) : _min($wq,$wt) );
      next if( (100 * $taln / $alnmax ) < $TINYALN); # skip trival matches
      $nmatch++;
      
      # 2013.aug: IS_CDSALIGN add $torient here? only care if $tor < 0; add sign to one of these cols?
      ## cant use taln, will screw up sort -n by maxalign; use tbit, unused now in scoring
      if($IS_CDSALIGN and $torient<0) { $tbit= -$tbit; } #? what followon problems does this cause?
      
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
      $validids{$lq}++; $validids{$lt}++;
      } 
    else { 
      puts($lq,$lt,$taln,$tidn);   $nmatch++; # DROP this old version
    }
  } 

if(DUPFILTER2) {  
  if($nmatch==0 and $dupfirst and $dupids{$lq} and $dupids{$lq} ne $dupfirst) { $nmatch=-1; $ndupdrop++; }
}
  if($nmatch==0) {
    if($OUTSPANTAB) { 
      my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt)= (0) x 9; 
      my $lt="self"; # or lq, or use blast-self scores?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $wq= $trsize{$lq}||0; 
      $at=$aq; $wt=$taln=$tidn=$tbit= $wq; # or aq
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
      $validids{$lq}++; $validids{$lt}++;
    } else {
    
    }
  }
  
  %bspans=();
  return($nmatch); # ,$ndupdrop
}



=item identityclass

  from  aabugs4qual/tsaevg/cdsidentclass.sh
  input from above putspans, sorted by ^Qclen, ^Align, vTtlen:
     qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits) 

  update classes per aabugs4qual/tsaevg/classparts.pl
  .. need blastp/blast table input for this..
  .. need 3-4 final categories:  keep-main, keep-alts, trash-fragments (alt,noalt) 
     .. keep includes partial prot where good enough.
     .. trash-alts maybe (not fragment but not distinct enough)
     .. add poor cut to main,alt,noclass per below a,b qualities.
     
cl_main.alttr and cl_alts, keep if: 
  a. homology unique or best (near best for alt),  or
  b. aasize >= 45% of main and aacomplete? and not aasame-as-main?
  c. main == fragment if aasize < MINAAFULL or (aasize < MINAAPART and partial/utrpoor)

FIXME/Check : aablast score has dropped/poor a set of 653 uniq ref (Nref1) matches vs prior clidtab/class1
  see env swapids=1 refblast=aaeval/refdebl/refde-.tall3 ./classparts.pl trsets/$pt.clid2tab
Partition               Nclass  Nmatch  Nrefa   Nref1   Bits    Iden    Algn    Rlen    Tlen
locust1best5.cl_poor    68122   1331    3485    653     113     61      170     361     389

# # FIXME2: maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; 
# # but refbest/good will keep those w/ homology
# cat $pt.class3 | grep 'drop.main' | sort -k6,6nr | head
# fungrvel4ik53Loc89t53   drop    main    fungrvel3bk43Loc142t24  99/66   1408,27%,complete-utrbad        187.2,DRERI:ENSDARG00000086955
#.. ^only aa equiv is also drop:
#  fungrvel3bk43Loc142t24  drop    althi   fungrvel4ik53Loc89t53   99/66   1373,32%,complete-utrbad        181.4,DRERI:ENSDARG00000086955
#
# fungrvel4k25Loc3571t1   drop    main    fungrvel4k29Loc22240t1  99/100  668,28%,complete-utrbad 0,0
# fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
# .. other DR92 genes:
#  fungrvel3bk35Loc16459t1 okay    main    fungrvel3bk29Loc13518t1 99/92   876,82%,complete        543,DRERI:ENSDARG00000068192,refbest
#  fungrvel3bk29Loc13518t1 okay    althi   fungrvel3bk35Loc16459t1 99/92   836,61%,complete        572,DRERI:ENSDARG00000068192,refbest
#  fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
#  fungrvel4k35Loc25356t1  okay    althi   fungrvel4k25Loc15515t1  100/65  363,99%,partial3        256.3,DRERI:ENSDARG00000068192
#
# fungrvel3bk29Loc5288t1  drop    main    fungrvel4k25Loc5314t2   99/87   565,23%,partial5-utrbad 0,0
# fungrvel4k25Loc3268t8   drop    main    fungrvel4ik53Loc14576t1 100/58  526,24%,complete-utrbad 396,DRERI:ENSDARG00000038737
#...    
    

=cut

sub classifytr {
  my($tid,$cla,$qid,$pal)= @_;

use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, };
use constant NOTPOORBAD => kAATINY + kAADUP; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;

  # change to aablastref{tid} ?
  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  ## FIXME: aamain can be bad/missing; fix before this, need cds/tr align info to know if aamainid is good.
  
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  ## add these utr class quals : this should be in tqual string: aasize, pcds, aaqual
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);
  
  ## main class should use aasize, aapartial and utrbad/poor with AAMINBAD, AAMINPOO
  ## FIXME: pcds bad here? comes from aw - gaps effects; separate gapbad from utrbad
  if($tqual =~ m/utrbad/) { $ispoor |= kAAUTRBAD; } elsif($tqual =~ m/utrpoor/) { $ispoor |= kAAUTRPOOR; }
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  else { 
    if($pcds < $BAD_CDSUTR) { $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }
  
  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    ## check that aamainid is main ??
    $cla.="a2" unless($cla=~/hi1/); # hi1/a2 redundant info
    $tbits.=",aadup:$aamainid";
    $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
  }
  # TEST3 : add? ispoor for cla == althi1 and pid/pal == 100/100  ie identicals?? NO, pal is align/shortlen
  # NO, not this, althi1 can be good alt, need align/mainlen stat also. or use only AADUP score as now.
  
  # CHECKME: adding aablast kept 40k more okay, all althi/ahia2 + 2k parthi
  # .. are these true uniq aablast or just althi aadups ?
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop="";
  if($tbrefbest eq $tid) { $keepdrop.="okay:refbest,"; $tbits.=",refbest"; $ispoor=0; }
  elsif($tbref2 and $tbrscore2 >= $MINBLASTSCORE) { 
    $keepdrop.="okay:refok,"; $tbits.=",refok"; $ispoor=0; }
  elsif($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore) { 
    $keepdrop.="okay:refgood,"; $tbits.=",refgood"; $ispoor=0; } # maybe2 ok ?

  ## maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; but refbest/good will keep those w/ homology
  ## keep largeaa,utrpoor for main, noclass, altlow; but drop smallaa,utrpoor
  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  if($cla =~ /parthi|frag0aa/) {
    $keepdrop.= "drop";
    
  } elsif($cla =~ /althi/) {
    $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
    
  } elsif($cla =~ /main/) {
    # $keepdrop.= "okay"; # always ok if has alts? NO, 
    # Fixme: keeping all "main" class gives ever expanding trsets with added trasm
    # should apply aaqual, aaref criteria to main + its althi; 

## 21feb update: dmag5icd35all_cde35.class4
## class7 maindrops: utrpoor/bad:48144, then partial:19493, gapbad:2866,  397 are aalong+complete
## class6 maindrops: utrpoor/bad:57636, then partial:19815, gapbad:3178,  1075 are aalong+complete
# dmag4vel4ibnk31Loc11665t1       drop    main    dmag4vel4ifik31Loc13364t1       99/99   102,71%,complete        0,0,pflag:4
# dmag4vel4ibnk31Loc1189t4        drop    main    dmag4vel4ipak31Loc21101t1       100/35  364,62%,complete        0,0,pflag:18
## pflag:18 = kAAGAPS + kAAUTRBAD
#..
## class4 maindrops: utrpoor/bad:60742, then partial:19737, but 3000 are aalong+complete, WHY drop?
## class5 maindrops: utrpoor/bad:56971, then partial:18807, but 1788 are aalong+complete (have gaps, so utrbad/poor now)
# dmag4vel4ibnk31Loc10846t1       drop    main    dmag4vel4ipak31Loc13986t1       99/94   116,73%,complete        0,0,pflag:4
# dmag4vel4ibnk31Loc11089t1       drop    main    dmag4vel4ifik31Loc18962t1       99/68   129,70%,complete        0,0,pflag:4
# >dmag4vel4ibnk31Loc10846t1 is aacluster:main, not aatiny, not utrpoor, not kAADUP, not althi;
# .. its alt is poor: partial5-utrbad, 37aa; did AAcluster fix mangle main?  pflag:4 == kAAUTRPOOR, miscalc from $pcds ???
# ** GAPS in aa above; 116aa,73%pcds is wrong if -gaps removed. BUT gaps not removed from nt size, so pcds off by that.
#  dmag4vel4ibnk31Loc10846t1 68aa,48gap,475nt,pcds=3*68/475=43%;  pcds= 3*68/(475-3*48) = 62%
# .. so drop.main/aacomplete is ~correct calc given aagaps : probably dont want utrbad flag here, but a gapbad may be useful
# .. 48gaps/116aa is not good prot; add flag gapbad == kAAGAPS above in readAAqual ?


    # $ispoor = $ispoor & kAATINY; # clear bits 2,4,.. # now above
    $keepdrop.= ($ispoor)?"drop":"okay";  # ??
    
  } elsif($cla =~ /noclass/) {
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low 
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  $tbits.=",pflag:$ispoor"; # DEBUG info
  return (wantarray) ? ($tid,$okay,$cla,$qid,$pal,$tqual,$tbits) : $okay;
}


sub identityclass {
  my($outh, $infile, $insorted)= @_;
  
use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing

## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field
### -aln == reverse align ?? for CDS antisense problem
### -bits == reverse align ?? bits not used here for scoring.. best choice other than adding column
### -aln affects here sort: -k7,7nr; use other field Ident?

#... see above now
#   $TINYALN=$ENV{mina}||50; 
#   $IS_CDSALIGN=$ENV{cdsw}||0; # default using tralign, tests better than cds
#   $PHIALN=$ENV{ahi}||98; # was 65; # use mina?
#   $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; 
#   # $ALTFRAG: add isfrag pct; now 50 (0.5)
  
### infile is putspans() alntab: Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
###  sort ^Qclen, ^Align, vTtlen: cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n 

  my($inh, %class,%bestmatch,%ismain, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
  	# ($ok,$inh)= openRead($infile,1);
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort -k2,2nr -k7,7nr -k6,6n $infile |") or die "reading $infile"; $inh=*IN; }
  }
  
  my($lastd)=("");
  while(<$inh>) { # maybe local table, sorted, or from file
    next if(/^Qid|^\W/); chomp; 
    my @v= split; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
    my($isfrag,$aclass,$ww,$pid,$pal,$samesize)= (0) x 10;
    
		## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field: not aln, iden?, bits?
    ## FIXME: *** 98/100/-sense/PitaEaR000975t24 << buggers, mixes up users of this table.
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    
    # my $antisense=($IS_CDSALIGN and $bits < 0) ? -1 : 0; 
    #old# my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : ""; ## append to bestmatch ??
    my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : "/."; ## append to bestmatch ??
    		## bestmatch="id,99/89/-sense" for antisense?
    		## ^^ is this new flag causing problems?  Besthit is appended /after, parsed where?
    
    $isfrag= $aclass="";
    $samesize=($tc == $qc)?1:0; 
    if($tc > $qc){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?

		if($lastd and $lastd ne $qd) { # check for $ltd missing eqgene 
			if($neqgene>0 and $eqgenes->{$lastd}) {
				my @ted= sort keys %{$eqgenes->{$lastd}}; 
				my $miss=0;
				foreach my $te (@ted) { unless($havepair{$lastd}{$te}) { 
					my $palmap= $eqgenes->{$lastd}{$te}||0; 
					if($palmap) {
						$miss++; 
    				my $pidalnval="99/$palmap";  
    				$bestmatch{$lastd}="$te,$pidalnval" unless($bestmatch{$lastd});
    				$bestmatch{$te}="$lastd,$pidalnval" unless($bestmatch{$te});
    				my $laclass="altmap"; $class{$te}= $laclass; #?? is this going to work
						}  
					} }
				if($miss) { } #.. reclass $lastd ??
			}
		}
		
		$havepair{$qd}{$td}++; $lastd=$qd; # here?
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; # $aclass="noclass"; << should be this, but ??
      $bestmatch{$qd}="$td,100/100" unless($bestmatch{$qd});
      next;  # can lose eqgene if no further td/qd ..
    }

#FIXME: have 2nd perfect matches, e.g qd1 >> td but qd2 == td; need that class to drop dups.
# ** BETTER: Remove these before align/cluster; fastanrdb on .cds, .aa?
#  asmrna5/cdsx/alldmag5x.cds n=4093166; alldmag5x.nrcds n=1782663; 501494 have dups (some many-many)
# .. info is in blastn aligns, not in cdhit clusters (no 2ndary aligns).
# .. use alntab ident column, if cdssize1 = cdssize1 = ident-3 (stop), identicals
# eg dmag4vel4xbxk55Loc9866t1  dmag4vel4xfik55Loc11404t1 dmag4vel4xpak25Loc1083t9 : cds-ident alts of main dmag4vel4ibxk55Loc9119t1
#
    
#below# if($class{$td}) { next; }  # problem here for qd1,qd2 > td, need to keep qd2 in bestmatch{qd}

#       # add here?? $aaqual{$qd} ; $aaqual{$td}; .. maybe do this below, END
#     my($qqual,$tqual)= ($aaqual{$qd},$aaqual{$td});
#     my($qbits,$tbits)= ($aablast{$qd},$aablast{$td});

    ## FIXME: check qw-main pCDS for UTRBAD/POOR; dont call tw frag if qc/qw is UTRBAD/POOR
    ## add these utr class quals
    my $qcds= ($qw>0) ? 100*$qc/$qw : $OK_CDSUTR;
    my $tcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    my $qutrbad= ($qcds >= $OK_CDSUTR)?0:1; # ($qcds > $BAD_CDSUTR)?1: 2;
    my $tutrbad= ($tcds >= $OK_CDSUTR)?0:1; # ($tcds > $BAD_CDSUTR)?1: 2;

# FIXME: when qc == tc, can assign althi1 = althi2 instead of main1 = althi2
    if($IS_CDSALIGN) { $ww=($qc>$tc and $tc>0)?$tc:($qc>0)?$qc:$tc;  $isfrag= ($tc < $ALTFRAG*$qc)?"frag":""; } 
    else { $ww=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw;  $isfrag= ($qutrbad==0 and $tw < $ALTFRAG*$qw)?"frag":"";  }
    ## if($qc < $MINCDS) { $class{$qd}.="tiny"; } elsif($qc < $MINCDSPART and $ispartial<needqual){ ..}
    
    $pid= ($aln<1)?0: int(100*(0.5+$iden)/$aln); $pid=100 if($pid>100);
    $pal= ($ww<1)?0 : int(100*(0.5+$aln)/$ww); $pal=100 if($pal>100);
    ##WrongWayPeachy#my $tinyalndiff= (TEST3 && (($aln - $ww) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal
    my $tinyalndiff= (TEST3 && (($ww - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal

## PROBLEM using eqgenes .. wont get here unless td x qd is in blastn.alntab ; some are NOT.
## check all $eqgenes->{$td} ?? should have self match.
		## 13Sep01: add eqgene map info for missed alignments : pal adjustment
		## $eqgenes->{$tid}{$sid} >= $MIN_EQGENE alignment
		if($neqgene>0) { 
			my $palmap= $eqgenes->{$qd}{$td}||0; 
			if($palmap > $pal) { $pal=$palmap; } ## FLAG it somewhere.. aclass="altmap" ?
		}
		
    
    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});

# Maybe problem here making too many fragment alts: a,b near identical alts of c, b <= a subset,
# but both classed/bestmatch to c main.  b should be dropped as fragment/99equiv of a, but dont see that due
# to b<c precedence.

unless(TEST3) {
    if($class{$td}) { next; }  # TEST3 maybe defer this to afer aclass, change for althi1 ?
}    
    if($samesize and $class{$qd} and $class{$qd} =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td/);
    }
      ## reclassify? althi > althi100 for 2ndary alts; althi100 = unneeded duplicates : TEST
      ##if($class{$td} eq "althi" and $pid>99 and $pal>=99) { $class{$td}= "althi100b"; }
      # ^ this probably wont happen, 2ndary althi, 1st will also be 100 pi/pa
      # ** problem using cdhits clusters here, it doesnt give 2ndary aligns as blastn does
      # .. so partial align to big tr may also be perfect align to smaller, valid alt tr
      ## next;   # problem above for qd1,qd2 > td, need to keep qd2 with bestmatch{qd}
      
    if($pal < $TINYALN) { } #defer: $aclass="noalign$isfrag";  
    elsif($tc == 0 and $qc > 0 and $pid >= $PLOW) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier
    	# TEST3: add min-basediff here to PHIALN tests, so shorty tr dont take over w/ minor base diff
      $aclass .= "1" if($pid >= 99 and ($pal >= $PHIALN or $tinyalndiff));  # "100"; use or not? PHIALN reused: was 99 here
      #old# $aclass .= "1" if($pid > 99 and $pal >= 99);  # "100"; use or not? PHIALN reused: was 99 here
      }
    #old# elsif( $pid >= $PHI ) { $aclass= ($pal<$PHIALN or $isfrag)?"parthi":"althi"; } # Primary classifier
    #older# elsif( $pid >= $PHI ) { $aclass="althi$isfrag"; }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }

if(TEST3) {
    if(my $ocl= $class{$td}) {    
      ## $class{$td}= $aclass if($ocl eq "althi" and $aclass eq "althi1");
      $class{$td}= $aclass if($ocl =~ m/^alt/ and $aclass eq "althi1"); # is this right?
      next;    
    }
}    

    if($aclass) { 
      ## FIXME ismain{qd} >> qd may be alt of other main, .. follow chain here?  keep qd as bestmatch? both?
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";

use constant FINDMAINID => 1;  
## FIXME: this still leaves some alts w/o final main id link; including circular alt1 <=> alt2 links
if(FINDMAINID) {    
      if($class{$qd}) {  #  =~ /alt|part/
        my $qnext= $qd; my $more= ($bestmatch{$qnext})?1:0;
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; $attr="$pidalnval/$qmain"; } #was $attr.="/$qmain";
        ## WARNING: now attr has /-sense or /. before /qmain; 
        ##   evgmrna2tsa2.pl needs to know this field's structure, /qmain at end esp.
      }
}
      
      $ismain{$qmain}++; # $ismain{$qd}++; 
      $class{$td}=$aclass; $bestmatch{$td}="$qd,$attr"; # ",$pid/$pal$antiflag ??
      ## maybe add to prevent 2+ circular alta <> altb <> altc with no main ..
if(TEST3) {
      $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..
}
    }
      
  } close($inh);

  
#END:  add more fields to output: aaqual, aablast, tr,aa sizes?
# FIXME: here, elsewhere create ID main,alt table, with new numeric IDs, old/cur ids, main/alt num
 { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; 
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "main";    
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}

# pre-read all of infile ids into validids, before idenityclass
sub readIdsFromAlnTab {
  my($infile, $insorted)= @_;
  my($nids,$inh)=(0,undef);
  # ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    next if(/^Qid|^\W/); chomp; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= split; 
    
    #?? use dupids here to mark as not validids ?
    if( $dupids ) {
      $validids{$qd}++ unless( $dupids{$qd} and $qd ne $dupids{$qd});    
      $validids{$td}++ unless( $dupids{$td} and $td ne $dupids{$td});    
    } else {
      $validids{$qd}++; $validids{$td}++;
    }
    $nids++;
  } close($inh); # rewind if ref() ???
  return $nids;
}


sub readblasttab {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq,$nids)= (0) x 10;
  # my($ok,$fh)= openRead($bother);
  my $fh;
  if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open(F,$bother) or die $bother; $fh=*F; }
  %bspans=();
  my %dupskipids=(); my $dupskipspan=0;
  
## FIXME: self-only matches need  recording via putspans; should be but dupdrop is droping firstid also

  while(<$fh>) { 
    unless(/^\w/) { next; } 
    #  if(/^# Query/) { ## Unused info
    #  #Notnow# puts($lq,$lt,$sa,$sm) if($lt); 
    #  ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; }
      
    my @v= split; 
    my($q,$t,$bits,$aln,$mis,$gap,@bspan)= @v[0,1,-1,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 

    # FIXME: need to parse align parts before summing;
if(SPANSUM) {
    if($lq and $q ne $lq) { # dupskip ???
      putspans($lq) unless($dupskipspan); $nids++; %bspans=(); $dupskipspan=0; 
    }
} 

  ## FIXME: self-only matches need  recording via putspans
    my $dupskip=0;
if(DUPFILTER1) {    ## see DUPFILTER2
    # checkme: is dupids mainid always in blast set? if not dupskip drops useful items
    # ** DUPID list has ids not in blast/cluster set...
    #   add dupskipids{id}=mainid and check later in validids
    if( $dupids ) { ## and ($dupids{$q} or $dupids{$t})
      # if($dupskipspan) {
      #  $dupskip=1;
      # } elsif
      if( $dupids{$q} and $q ne $dupids{$q} ) { 
        $dupskipspan=$q; ## if($q eq $t); # flag to skip putspan
        $dupskipids{$q}++; $dupskip++;   
      } elsif( $q ne $t and $dupids{$t} and $t ne $dupids{$t}) { 
        $dupskipids{$t}++; $dupskip++; }  
      if($dupskip) { $ndupdrop++; }  # next; below, not here  WRONG now, need to skip self putspan unless this dup is main 
    }
}   
    ## UNUSED now: qd, wq : FIXME save ids of no-match
    if($t eq $q) { } ## { $qd=$q unless($qd); $wq=$aln unless($wq); }
    elsif($dupskip==0 and $dupskipspan == 0) { 

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
  putspans($lq) unless($dupskipspan); $nids++; $dupskipspan=0;
} else {
  my $aident= _max(0,$sa-$sm); # other way to calc: $aident = $pctident * $aln;
  puts($lq,$lt,$sa,$aident) if($lt); 
}  

  if($ndupdrop>0) { $ndupdrop=scalar(keys %dupskipids); }
  warn "# readblasttab: nids=$nids; ndupdrop=$ndupdrop\n" if $debug;
  return($nids);
}

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
## 2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
sub sumblastpart {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0; # 2013.aug : ORIENT problem.  if both here are reversed, or=0; if only 1, or=-1
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or=-1; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or=($or<0)?0:-1; } #was $or--
  unless($bspans{$t}) { 
    $bspans{$t}=[]; push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); 
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
  ## IS_CDSALIGN add $or to bspans, problems?
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); }
}


=item readlastz: lastz align general format

test case:
  aabugs4/tsaevgc/daphmag5xbest5/dmag5xau13c2011f
  pt=dmag5xau13c2011

$evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 \
  -aasize inputset/$pt.aa.qual -acdhit tmpfiles/${pt}_cd90.aa.clstr \
  -lastz tmpfiles/$pt-self97.lastz.gz  \
  -outeq tmpfiles/$pt.alnlztab2 -outclass $pt.trclasslz2

$bg/mb/galn/bin/lastz \
 --identity=100 --coverage=25 --step=10 --seed=match15 --notransition --exact=20 --match=1,5 \
 --ambiguous=n --nochain --nogapped  --format=general 'altset1.okboth.cds[multiple]' altset1.okboth.cds

#score  name1              strand1 size1 zstart1 end1  name2               strand2 size2 zstart2 end2 identity idPct  coverage covPct

779     dmvel4xpak25Loc11378t4   +  1083    0    779   dmvel4xpak25Loc11378t10  +   1083  0    779  779/779 100.0%  779/1083  71.9%
1083    dmvel4xpak25Loc11378t10  +  1083    0    1083  dmvel4xpak25Loc11378t10  +   1083  0    1083 1083/1083 100.0%  1083/1083 100.0%
770     dmvel4xpak25Loc11378t5   +  978     208  978   dmvel4xpak25Loc11378t10  +   1083  313  1083 770/770 100.0%  770/978 78.7%
251     dmvel4xpak25Loc11378t1   +  684     0    251   dmvel4xpak25Loc11378t10  +   1083  399  650  251/251 100.0%  251/684 36.7%
330     dmvel4xpak25Loc11378t1   +  684     354  684   dmvel4xpak25Loc11378t10  +   1083  753  1083 330/330 100.0%  330/684 48.2%
324     dmvel4xpak25Loc11378t7   +  753     429  753   dmvel4xpak25Loc11378t10  +   1083  753  1077 324/324 100.0%  324/753 43.0%
    ...
455     dmvel4xpak25Loc11378t3   +  753     0    455   dmvel4xpak25Loc11378t7   +   753   0    455  455/455 100.0%  455/753 60.4%
753     dmvel4xpak25Loc11378t7   +  753     0    753   dmvel4xpak25Loc11378t7   +   753   0    753  753/753 100.0%  753/753 100.0%
456     dmvel4xpak25Loc11378t1   +  684     222  678   dmvel4xpak25Loc11378t7   +   753   297  753  456/456 100.0%  456/684 66.7%
324     dmvel4xpak25Loc11378t10  +  1083    753  1077  dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%
324     dmvel4xpak25Loc11378t5   +  978     648  972   dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%

=cut

sub readlastz {
  my ($bother)= @_;
  
  my $islzg=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  # my($ok,$fh)= openRead($bother);
  my $fh;
  if($bother =~ /\.gz/) { open( $fh,"gunzip -c $bother |") or die $bother; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open($fh,$bother) or die $bother;  }
  %bspans=();

## TEST: is lastz output file sorted properly? LOOKS OK //name2 should all be grouped, does split-run undo that?
  
  my @hd; # $islzg=1 ;
  while(<$fh>) { 
    if(!$islzg and m/^#score\tname1/ and /\tcoverage/) { $islzg=1; chomp; s/^#//; @hd=split"\t"; } # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    unless(@v==15 and $islzg) { die "# ERR: doesnt look like lastz general format I know: hd=",@hd," val=",@v,"\n" ; }
    
    ## allow subset columns?
    #score  name1   strand1 size1 zstart1 end1 name2   strand2 size2   zstart2 end2  identity idPct   coverage   covPct
    my( $lzscore, ## is lz score = count of ident bases? no,
      $tid, $tor, $tsize, $tstart, $tend,
      $qid, $qor, $qsize, $qstart, $qend,
      $nida, $pid, $ncovb, $pcov)= @v;
    $qstart++; $tstart++; # move to 1-origin

    ## NOTE: my lastz out has name2 as first-order == qid, name1 == tid
    if($lq and $qid ne $lq) {  
      putspans($lq);  $nids++; %bspans=();
    }

    if($tid eq $qid) { 
      #?? add:  $trsize{$qid} = $qsize unless($trsize{$qid});
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      my($aident,$na)= split "/",$nida;
      my($aln,$nb)= split "/",$ncovb;
      my $bits=  $lzscore; # what? maybe 2*aident? 
      
if(1) { ## SPANSUM
      ## *?* Need this; lastz hsp as for blastn can overlap lots ..
      sumblastpart( $qid, $tid, $bits,$aln,$aident, $qstart,$qend,$tstart,$tend); 
} else {      
      $bspans{$tid}=[]  unless(defined $bspans{$tid}); 
      push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]);  
}
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq);  $nids++;   
  warn "# readlastz: nids=$nids\n" if $debug;  # ; ndupdrop=$ndupdrop
  return($nids);
}



=item cd-hit(est) align cluster format

  has align-span, percent ident, enough for use.
  
  cacao11pub3ig_cde25.cds.clstr
  >Cluster 0
  0       16221nt, >Thecc1EG029098t1... *
  >Cluster 1
  0       15408nt, >Thecc1EG019010t1... at 1:15408:1:15411/+/99.98%
  1       15411nt, >Thecc1EG019010t2... *
  >Cluster 2
  0       14343nt, >Thecc1EG007738t1... *
  1       10578nt, >Thecc1EG007738t2... at 1:10574:1831:12404/+/100.00%
  >Cluster 4
  0       12732nt, >Thecc1EG015810t1... at 1:12732:304:13035/+/100.00%
  1       13035nt, >Thecc1EG015810t2... *
  2       12504nt, >Thecc1EG015810t3... at 74:12503:148:12584/+/99.94%
  3       12717nt, >Thecc1EG015810t4... at 1:12717:304:13035/+/99.88%
  >Cluster 15
  0       9804nt, >Thecc1EG034527t1... *
  1       9804nt, >Thecc1EG034527t2... at 1:9804:1:9804/+/100.00%
  2       9804nt, >Thecc1EG034527t3... at 1:9804:1:9804/+/100.00%
  3       7512nt, >Thecc1EG034527t4... at 1:7512:2287:9804/+/99.92%
  >Cluster 21
  0       8751nt, >Thecc1EG006991t1... *
  1       8304nt, >Thecc1EG006991t2... at 2894:8304:3336:8751/+/99.83%
  2       8178nt, >Thecc1EG006991t3... at 2894:8178:3336:8630/+/99.74%
  3       8391nt, >Thecc1EG006991t4... at 3239:8377:3356:8494/+/100.00%

  dmag4vel4xfi_cde60.cds.clstr
  >Cluster 8633
  0       1383nt, >dmag4vel4xfik65Loc51t4... *
  1       1383nt, >dmag4vel4xfik75Loc49t1... at 1:1383:1:1383/+/100.00%
  >Cluster 8634
  0       1383nt, >dmag4vel4xfik65Loc96t9... *
  1       1383nt, >dmag4vel4xfik65Loc96t15... at 1:1383:1:1383/+/99.42%
  2       1383nt, >dmag4vel4xfik65Loc96t18... at 1:1383:1:1383/+/99.86%
  3       684nt, >dmag4vel4xfik81Loc3813t5... at 1:684:463:1146/+/99.12%
  4       864nt, >dmag4vel4xfik85Loc1522t3... at 1:864:379:1242/+/100.00%
  5       864nt, >dmag4vel4xfik85Loc1522t5... at 1:864:379:1242/+/99.31%
  6       774nt, >dmag4vel4xfik91Loc1049t2... at 1:774:373:1146/+/99.22%
  >Cluster 8635
  0       873nt, >dmag4vel4xfik45Loc1t9452... at 1:873:70:942/+/99.89%
  1       873nt, >dmag4vel4xfik45Loc1t9453... at 1:873:70:942/+/99.89%
  2       483nt, >dmag4vel4xfik45Loc1t9455... at 1:483:901:1383/+/100.00%
  3       483nt, >dmag4vel4xfik45Loc1t9479... at 1:483:901:1383/+/100.00%
  4       492nt, >dmag4vel4xfik55Loc412t36... at 1:492:892:1383/+/100.00%
  5       1383nt, >dmag4vel4xfik65Loc127t19... *
  6       240nt, >dmag4vel4xfik65Loc127t21... at 1:240:1144:1383/+/100.00%
  7       1383nt, >dmag4vel4xfik75Loc138t9... at 1:1383:1:1383/+/99.86%

=cut

sub readcdhit {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  # my($ok,$fh)= openRead($bother);
  my $fh;
  if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr,$hasdupid)=(0)x10;
  while(<$fh>) { 
    if(/^>/) { 

    ## move this dupid filter before into readblastab, readcdhit ?
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
      putspans($mainid) if($mainid); %bspans=(); $hasdupid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/) # problem: >(ID\.1)\.\.\.
    } elsif(/^(\d+)/) { 
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/; # problem: >(ID\.1)\.\.\.
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*$/)?1:0;
      $hasdupid++ if($dupids and $dupids{$tid}); # need to read all cluster then drop dups
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push( @{$bspans{$tid}}, [1,$tlen,1,$tlen,$tlen,$tlen,$tlen]); 
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        my($qstart,$qend,$tstart,$tend)= (0) x 4;
        my @pinfo= split( m=/=, $pinfo);
        my @aln= split/:/,$pinfo[0]; #? multiple align segs or 1 only?
        my $aor=$pinfo[1]; # dont need?
        # 249nt, >dmag5vel5xco1k75Loc13984t1... at 249:1:421:669/-/100.00% << revalign
        if(@aln>3) { ($qstart,$qend,$tstart,$tend)= @aln; # what if @aln>4 ?
          ($qstart,$qend)= ($qend,$qstart) if($qstart>$qend); } 
        my $aln= 1+$qend-$qstart;
        my $aident= int($aln * $pi/100);
        my $bits= $aln; # aident? 
        push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); # add tlen?
        $nalt++;
      }
      # $iscdhit=2; 
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
  putspans($mainid) if($mainid); %bspans=();
  warn "# readcdhit: nclust=$ncl, nalt=$nalt, ndupdrop=$ndupdrop, nerr=$nerr \n" if $debug;
  return($ncl+$nalt); 
}


# change mainaa/subaa in %aacluster, using other input info: aaqual/size-gaps, cds/tr align ids
sub aaqualscore
{
  my($mqual)= @_;
  my $mqv=0; $mqual ||="missing";
  if($mqual =~ /complete/) { $mqv += 2; } elsif($mqual =~ /partial[35]/) { $mqv += 1; }
  if($mqual =~ /utrbad/) { $mqv -= 2; } elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  return $mqv; # range is -2 .. +2
}


sub correctAAcluster
{
  my($havevalidids)= @_;
  # % aaclustermain == hash{mainid}[subids]
  # % aacluster == hash{eachid} = "mainid,pctident"
  my $nreset=0; my %didmain;
  foreach my $mid (sort keys %aaclustermain) {
    next if($didmain{$mid}); # probably dont need. 
    my $mainid= $mid;
    my @cids= @ { $aaclustermain{$mid} }; # has all ids incl mainid?
    my $maw = $aasize{$mid} || 0;
    my $mqv = aaqualscore($aaqual{$mid});  #? parse for "aasize,pcds,aaqual"
    ## sort cids by aasize? or need to check thru all?
    # @cids = sort{$aasize{$b} <=> $aasize{$a}} @cids;
    my $reset=0;  my @goodids=(); my $ninval=0; 
    if($havevalidids and not $validids{$mainid}) { $maw=0; $mqv=-9; $mainid=0;  $ninval++; }
    foreach my $id (@cids) {
      ## also check each id is valid for tr/cds align, drop invalids including current mainid
      if($havevalidids and not $validids{$id}) { 
        ## if($id eq $mainid) { $maw=0; $mqv=-9; $mainid=0; } # above now: main could be last id.. fixme
        $ninval++; next; # drop from aacluster
      }
      push @goodids, $id; 
      next if($id eq $mainid);
      my $aw= $aasize{$id} || 0; 
      my $qv= aaqualscore($aaqual{$id});
      if($aw > $maw and $qv >= $mqv - 1) {  # reset main; ignore qv if aw >>> maw ?
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++; 
      } elsif($qv > $mqv and ($aw > 0.98*$maw)) {
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++;       
      }
    }
    
    if( @goodids == 0 and $ninval >= @cids) {
      foreach my $id (@cids) { delete $aacluster{$id}; }
      $aaclustermain{$mid}= [];
      $nreset++; $reset=0; 
    }
    if($reset) {
      foreach my $id (@goodids) {
        my($oldmain,$pi)=split",",$aacluster{$id};
        $aacluster{$id}= "$mainid,$pi";
        }
      $aaclustermain{$mainid}= \@goodids;
      $nreset++;
    }
    $didmain{$mainid}++;
  }
  warn "# correctAAcluster nreset=$nreset\n" if($debug); # got nreset=212165 for nclust=232463 TOO HIGH?
  return($nreset);
}


sub readAAcdhit {
  my ($bother)= @_;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  # my($ok,$fh)= openRead($bother);
  my $fh;
  if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  else { open(F,$bother) or die $bother;  $fh=*F; }
  
  %aacluster=(); # global
  our @cluster=();
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr)=(0)x10;

  # FIXME: use aaqual utrbad/poor for main/alt with pi=100; reset main if utrbad and alt utrok
  sub putclus { 
    my($mainid,$cluster1)=@_; ## our @cluster;
    if($mainid) {  
      $aaclustermain{$mainid}=[] unless($aaclustermain{$mainid});
      foreach my $cl (@$cluster1) { 
      my($id,$pi)=split",",$cl; $aacluster{$id}="$mainid,$pi"; 
      push @{$aaclustermain{$mainid}}, $id;
      } 
    }
  }
  
  while(<$fh>) { 
    if(/^>/) { 
      putclus($mainid,\@cluster) if(@cluster); @cluster=(); $mainid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/)  ## BADD patt ??
    } elsif(/^(\d+)/) {
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1.f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      ## FIXME2: new merge fnum at end of line now;
      my @pmore; ($pinfo,@pmore)= split /\t/, $pinfo;
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      # FIXME: merge.aa.clster bad; lacks some main * ; skip those for now
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push @cluster, "$tid,$pi";
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%   <<<<<
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        push @cluster, "$tid,$pi";  $nalt++;
      }
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
  putclus($mainid,\@cluster) if(@cluster); ### putspans($mainid) if($mainid); %bspans=();
  warn "# readAAcdhit: nclust=$ncl, nalt=$nalt, nerr=$nerr \n" if $debug;
  return($ncl); 
}

sub readblatpsl {
  my ($bother)= @_;
  
  my $ispsl=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  # my($ok,$fh)= openRead($bother);
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
      putspans($lq); $nids++; %bspans=();
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
  
  putspans($lq); $nids++;  
  return($nids);
}







#################### DROP, Obsolete classifier ###################

sub outclusters {  # DROP this old version
  my($outh)= @_;
  
  # $better{$lq}{$lt}++;
  my @ids= sort keys %validids;
  my %sumscore;
  foreach my $id (@ids) { 
    foreach my $jd (@ids) { next if($id eq $jd); 
    my $v= $better{$id}{$jd}||0; $sumscore{$id} += $v;
    }
  }
  
  ## need header:
  # outrow: join("\t",$typ,$lq,"$aq/$wq",$lt,"$at/$wt",$sa,$diffaln,$sm)."\n"; 
  print $outh join("\t",qw(Cluster Score12 Type Qid Qsize Tid Tsize Align Daln Ident))."\n" unless($head++);
      
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
      print $outh "cl$cluid\t$topscore,$nextscore\t",$outrow;
      $didid{$nextid}++;
    }
  }
}

# need to do more before output: cluster all same/subset by ids
sub puts { # DROP this old version
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
    $validids{$lq}++; $validids{$lt}++;
    $outrows{$lq}{$lt} = $val; # if exists $outrows{$lq}{$lt}, *should* be same val, check?
    #NO# print $val;
    } 
}


__END__


=item classifytr new result1

# /bio/bio-grid/aabugs4/bugs/locust/bestof5x/trsets
# locust1best5 primary alt.tr of cl_main; n=108607 - 8545 in cl_alts,poor; 5446 cl_alts have cl_main
# .. need to filter out trival/aberrant alttr .. use trsize, aasize, close to main.tr
  108607 locust1best5.cl_main.altids
   25968 locust1best5.cl_main.altok.tab;  >=45% aasize of main;
   14754 are aa-complete (incl utrbad), not in cl_alts,poor
    5798 are same size as main; same prot?

   55376 locust1best5.cl_alts.ids
   33284 locust1best5.cl_main.ids
   40761 locust1best5.cl_poor.ids

** need to add other alt-set: needs to classify self-blast + aablast
/bio/bio-grid/aabugs4/bugs/locust/bestof5x/trsets
   locust1best5.cl_main.altids

/bio/bio-grid/aabugs4/tsaevg/trsets
      29 Feb  6 20:48 locust1best5.aa.qual -> .
 33933545 Jan 17 13:19 locust1best5.alntab
 31561602 Feb 13 23:18 locust1best5.alntab2
 1455946 Feb  8 14:24 locust1best5.cl_alts.ids
  853266 Feb  8 14:24 locust1best5.cl_main.ids
 1055964 Feb  8 14:24 locust1best5.cl_poor.ids
 12266926 Feb 13 23:18 locust1best5.class2
 5581623 Jan 19 22:57 locust1best5.classtab
 5939924 Feb  6 11:28 locust1best5.clid2tab
 6166194 Jan 14 15:22 locust1best5.tr.count

# test new classifier
pt=locust1best5
$evigene/scripts/rnaseq/asmrna_dupfilter2.pl -outspan -outeq $pt.alntab2 \
-aa=$pt.aa.qual -tr=aaqual \
-blast=sdoutz/sd-slf95-$pt.blastn.gz  \
-ablast=../aaeval/refdebl/refde-$pt.tall3 > $pt.class2

cat $pt.class2 | cut -f2 | sort | uniq -c | head
71657 drop
1729 maybeok
54816 okay

cat $pt.class2 | cut -f2,3 | sort | uniq -c | head -30
5595 drop       althi
2642 drop       altmid
4919 drop       altmidfrag
27183 drop      noclass
31318 drop      parthi   : check drops for aablast, aaqual; altmid may be valid due to part align to keeper
1729 maybeok    parthi   : all rescued by aablast score; fix this for >=90% of best?
16731 okay      althi
3511 okay       altmid
3432 okay       altmidfrag
19763 okay      main
11379 okay      noclass


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
