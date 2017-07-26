#!/usr/bin/perl
# evigene/scripts/rnaseq/asmrna_trimvec.pl

=item about asmrna_trimvec 

  merge of veccutmrna2.pl and parts from evigene/scripts/evgmrna2tsa.pl
  - revised and complex script for mRNA with vector screen, removal and NNN gap trimming
  - detects vectors in mRNA assemblies (using vecscreen/blastn -d UniVec_core data)
  - trims end gaps, per NCBI TSA requirements
  - uses protein/CDS info with ref-protein to prevent spurious/poor screen/trim that damages good mRNA
  - pulled simpler version in evgmrna2tsa, tested in veccutmrna2 to work right

  - recomputes protein,cds of trimmed mRNA and merges for evgmrna2tsa uses;
  - logs ambigious/problem cases for expert inspection while trying to minimize this,
    tries to keep valid protein while trimming vectors/gaps.
    
=cut

use constant VERSION => '2013.06.01'; # 05.28

use FindBin;
use lib ($FindBin::Bin,"$FindBin::Bin/.."); # assume evigene/scripts/rnaseq/thisscript.pl

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_proteins;
use cdna_evigenesub;

# cdna_evigenesub globals:
our $EVIGENES="$FindBin::Bin/..";  
our $EGAPP='egtrimvec';  
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;
#drop# my $UVGAP=  $ENV{uvgap}||0;   ## NO: leave nnn gaps for cut vec, preserve codon%3 size?

## Evigene tr2aacds subdirs: see tidyup
## add for trimvec, mrna2tsa:  trimset?  publicset? tsasubmit/submitset ?
## change trimvec,mrna2tsa output subdir: publicset? adding pubids, main2alt, ...
## separate tsasubmit subdir ..
our @evgdirs = qw(okayset dropset inputset tmpfiles erasefiles publicset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset); # tidyup file sets

our $vecoutdir='publicset'; # trimset now ?


use constant{ HasGapNone => 0, HasGapVector => 1, 
    HasGapEnd5 => 2, HasGapEnd3 => 4, HasGapMaxSpan => 8, HasGapTooManyXs => 16, }; # hasBadGaps

my $UniVecDB= $ENV{UniVec} || "UniVec_Core"; # Not UniVec
my $MAXGAP=$ENV{maxgap}|| 15; # NCBI decides, changes..
my $ENDGAP=$ENV{endgap}|| 20; ## was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
my $MINSIZE_NCBI=200; # NCBI
my $MINSIZE= $ENV{min} || 150; # lower to check cuts w/ valid homology?
my $GAPSOK=1; # 2012-dec .. for now, NCBI may change again (2013-may?)
my $NCPU= 1; 
my $tidyup=0;

my(%vecscreen);
my($genenames,$cdnaseq,$trclass,$vecscreenf, $logfile, $VECSUF,$outsuf,$namesuf);  ## ,$sufaa,$sufcds,$sufuncut

## ARGS
## mixup: locust1all5asm.locust1all5asm.mrna.vector.tab2
$VECSUF= "vector.tab"; #fixed in  vecscreen: makename($cdnaseq,".vector.tab");
$outsuf= $ENV{outsuf}|| "uvcut"; # was uvcut.mrna; 
$namesuf= $ENV{namesuf}|| "names";

my $optok= GetOptions(
  "mrna|cdna=s", \$cdnaseq,
  "class|trclass=s", \$trclass,
  "vectors|vecscreen=s", \$vecscreenf,  
  "names|genenames=s", \$genenames, ## ? allow for 2 files: myspecies.namerefids + allrefprot.names
  "logfile:s", \$logfile,
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  "ENDGAP=i", \$ENDGAP,  
  "NCPU=i", \$NCPU,## "MAXMEM=i", \$MAXMEM,  
  "dryrun|n!", \$dryrun, 
  "tidyup!", \$tidyup, 
  "DEBUG!", \$DEBUG, 
  );

die "EvidentialGene asmrna_trimvec VERSION ",VERSION,"
mRNA vector screen/removal and gap trimming
  - detects vectors in mRNA transcripts using vecscreen/blastn -d $UniVecDB 
  - trim end gaps per NCBI TSA requirements
  - uses protein/CDS info with ref-protein to prevent spurious/poor screen/trim that damages good mRNA
  - recomputes protein,cds of trimmed mRNA and merges for evgmrna2tsa uses;
  - logs ambigious/problem cases for expert inspection while trying to minimize this,
    tries to keep valid protein while trimming vectors/gaps.
  
Usage: asmrna_trimvec.pl -mrna mrna.fasta OR -class name.trclass ...
opts: -genenames mrna.names -log  -debug
    -NCPU=$NCPU -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP\n"
  unless($optok and ($cdnaseq or $trclass));     

$tidyup= 1 unless($dryrun||$DEBUG); # default on unless debug|dryrun ?
my $GAPSMAX = ('N') x $MAXGAP;

openloggit($logfile,$cdnaseq||$trclass);
loggit(1, "EvidentialGene asmrna_trimvec.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

our $APPvecscreen= findapp("vecscreen");  
our $APPtraa2cds= findevigeneapp("prot/traa2cds.pl"); # move to cdna_evigenesub for get_mRNA
our $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions?
#-------------------------------------

##.... REWRITE HERE for single input .mrna + vector.tab + genenames ..........
#cdna_evigenesub: my( %genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames); 

=item steps for asmrna_trimvec 

  0. input.mrna of evigene tr2aacds 
    - requires okayset/name.mrna file creation (revcomp -strand); not done in tr2aacds

  q0? should outputs be relocated to new subdir? eg. tsasubmit/ or leave in okayset/
    - probably leave in okayset/ replacing old versions
    
  1. vecscreen() : run ncbi c-- vecscreen, alternately c++ blastn (later)
        using UniVec_Core (or other) db
  
  2. mrna_trimvec() : remove vector spans in mRNA, and end gaps, recalc orfs of trim set
  
  3. update_mrna_fileset : should this be located in folder or okayset/ ?
      ^^ move to evgmrna2tsa, so can make final pubset merging trimset + okayset + pubids + annot ..
      -- fixme : okboth.aa, .cds may not exist where okboth.mrna is input
      
  Next: evigene/scripts/evgmrna2tsa.pl  now looking for these outputs as input.
  
=cut

sub MAIN_start {}
MAIN: {
  loggit(0, "BEGIN with input=",$cdnaseq||$trclass,"date=",`date`);

  my($mrnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,$vecoutdir);
    loggit(0, "get_evgtrset=",$mrnaseq,$trpath,$trname); ## facount($cdsseqnr)
    loggit(LOG_DIE, "Missing -mrna",$mrnaseq) unless($mrnaseq and -s $mrnaseq);

	unless($genenames) { #? put in get_evgtrset
		my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
		loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
	}

  ($vecscreenf)= vecscreen($mrnaseq,$vecscreenf||"");

  ## mrnaseq = input .mrna; outf, outuncut = updated .mrna << DROP outuncut for update_fileset
  my($ntrim, $trimids, @trimfiles)  ## $outf, $outaa, $outcds
    = mrna_trimvec($mrnaseq,$vecscreenf,$genenames);

  ## FIXME: add tidyup: remove or move to tmpfiles 
  ##   vecscreen ncpu "_vecsplit/" if ok: -s mrna.vecscreen.tmp or -s mrna.vector.tab
  ## merge new uvcut files + old uvuncut:
  # $nam  .uvcut.mrna + .uvuncut.mrna >> replace input.mrna and flag action 
  # ditto: .uvcut.aa .uvcut.cds ; need pull uncut .aa,.cds from okayset/*
  # need file? of ids for uvcut, uvuncut ?

#... defer update_mrna_fileset() to mrna2tsa, push @trimfiles to trimset/ dir 
use constant DEFER_UPDATE_FILESET => 1;

  my $upstatus=0; my($upfiles,$uptemp,$tmpfolder); 
if(DEFER_UPDATE_FILESET) {
	push @tmpfiles, @trimfiles;  ## $uptemp = \@trimfiles; 
	$tmpfolder="trimset";
	$upstatus=1;
	$tidyup= 1; #? always
	
} else {
  ## FIXME in mrna2tsa?  publicset/mrna,aa,cds need >pubid not >oid 
  ## .. rewrite fasta hdr again? or make pubids before this? 
  ## need some adjustments for DROPped mrna to mainalt,pubid 
  
  # NOW: mrnaseq in $vecoutdir, get path for others from it .. but need okayset param
  ($upstatus, $upfiles, $uptemp)  
    = update_mrna_fileset($trpath, $mrnaseq, 'trimvec_done', $trimids, @trimfiles) if(@trimfiles>0); 
    # = update_mrna_fileset0($trpath, $mrnaseq, $ntrim, $trimids, @trimfiles) if($ntrim>0);#  $outf, $outaa, $outcds
	push @tmpfiles, @$uptemp if(ref $uptemp);
	#? push @publicset, @$upfiles if(ref $upfiles);
	$tmpfolder="tmpfiles";
}
  
  if( $tidyup and $upstatus > 0 ) { ## not: and -s $upfiles->[2] 
    ## fixme: $fn may have path: okayset/old.mrna ... fixme2: this only works when curdir = main evgr path, may not be
    # sub tidyup{ my($tod,@td)= @_;  mkdir($tod); foreach my $fn (@td)  
    #  { if(-f $fn){ (my $tfn=$fn)=~s,^\w+/,,; runcmd("mv $fn $tod/$tfn");} } 
    #}

    tidyupFileset($tmpfolder,@tmpfiles);  
    tidyupFileset("publicset",@publicset) if(@publicset);  #? not used
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}

#...................................
#... test loop
# my @inmrna= @ARGV;
# foreach my $mf (@inmrna) {
#   (my $nam = $mf) =~ s/\.mrna.*//;
#   my $vtab="$nam.mrna.$VECSUF";
#   my $genenamef="$nam.$namesuf"; # what? evg path is above okayset/mrna, evgmrna2tsa.pl gets it
#   mrna_trimvec($mf,$vtab,$genenamef);
# } # end inmrna

#---------------------------------------------------------------------------


# ($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,$outdir);
# from  evigene/scripts/evgmrna2tsa.pl:getmRNA() : should this be here?

sub get_evgtrset {
  my($trclass,$cdnaseq,$outdir)= @_;
  my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);
  my $notokay=0;
  
  if($cdnaseq) { 
    $notokay=1; # dont look in okayset/? look in $outdir now?
    $trclass= makename($cdnaseq,".trclass") unless($trclass); 
  }
  
  if($trclass) { # dont require this exists. just trpath/okayset
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # BADDDDD
    else { $trname=$trpname; }
    $trpath ||= '.';  
     
    ## ?? fixme for update_mrna_fileset : merge okay.aa, okalt.aa,.cds also
    #my $okpath= "$trpath/$outdir" if($outdir); # need to check both, 3 paths?
    my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?
    ($cdnaseq)= getmRNA($okpath,$trname,$outdir) if(!$cdnaseq and -d $okpath);
  }
  
  return($cdnaseq,$trpath,$trname,); ## $sradatah);
}

# sub getOkFileset ## moved to cdna_evigenesub.pm
# sub openRead # moved to cdna_evigenesub.pm

sub getmRNA_OLD   ## moved to cdnasubs.pm ?
{
  my($okpath,$trname,$pubdir,$ADDutrorf)= @_;
  my($cdnaseq)=(""); # == mrna (oriented), not cdna.parts.tr
  
	use constant ALSOMAKE_AACDS => 1;
  $ADDutrorf= 1; # what? always check for okayset/*.utrorf.mrna ?
  
  #? FIXME? suffix .tr may not be used: .cdna? .fa? .fasta? ...
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna  

  my($oktr,$alttr,$okd)= getOkFileset($okpath,$TRSUFFIX);
  my @okd= @$okd;
  my @pubd=();
  if($pubdir and -d $pubdir) { my($pubd)= getFileset($pubdir,$TRSUFFIX); @pubd= @$pubd; }
  
  ## another bad call: #egr: get_evgtrset= publicset/locust1all5asm.p4.mrna . locust1all5asm
  ##   instead of publicset/locust1all5asm.mrna .. my test version mistake..
  ## Ugh! bad call: #er2g: get_evgtrset= ./okayset/pogonus1all3.mrna0.ann.txt . pogonus1all3
  ## need grep /\.mrna$|\.mrna.gz$|.mrna.fa$/ ??? 
  
  my ($trf);
  ($trf)= grep /$trname\.(mrna$|mrna.gz$)/, (@pubd, @okd); # drop okd here?
  unless($trf) { ($trf)= grep /\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not?
  if($trf and -s $trf) { $cdnaseq= $trf; }
  else { 
    my $okall=0;
    my $cdnatmp= ($pubdir) ? $pubdir : $okpath;
    $cdnatmp .= "/$trname.mrna";  
    mkdir($pubdir) if($pubdir and not -d $pubdir);

## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??

    loggit(0,"Make mRNA $cdnatmp from okayset transcripts");
    my($okaa) = grep /.okay\.aa$|.okay\.aa.gz$/, @okd;  
    my($altaa)= grep /.okalt\.aa$|.okalt\.aa.gz$/, @okd; 
    # #.. problem here .aa.gz not .aa; should makename() do optionally?
    
    if($oktr and -f $oktr and -f $okaa) { 
      my $err= runcmd("$APPtraa2cds -trout -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    if($alttr and -f $alttr and -f $altaa) { 
      my $err= runcmd("$APPtraa2cds -trout -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    
    if($ADDutrorf and $okall > 0) {
 			my($okin) = grep /.utrorf.mrna$|.utrorf.mrna.gz$/, @okd;  
   		if($okin) { runcmd("cat $okin >> $cdnatmp"); loggit(0,"add $okin to $cdnaseq"); } # err check? loggit?
   		else { $ADDutrorf=0; } # dont do .aa,cds
    }
    $cdnaseq= $cdnatmp if(-s $cdnatmp);
    loggit(1,"ERR: No-make mRNA $cdnatmp from $oktr,$okaa + $alttr,$altaa") unless($okall>1);
    
    # FIXmaybe: also make pubdir/.aa,.cds along with .mrna ? see hassle in update_mrna_fileset 
    if(ALSOMAKE_AACDS and $pubdir) {
  		my($ok,$hin,$hout,$fout,$okin,$altin,$utrin);
  		foreach my $suf (".aa",".cds") {
				$fout= makename($cdnatmp,$suf);  ## $cdnaseq
				($okin) = grep /.okay$suf$|.okay$suf.gz$/, @okd;  
				($altin)= grep /.okalt$suf$|.okalt$suf.gz$/, @okd; 
				($utrin)= ($ADDutrorf) ? grep(/.utrorf$suf$|.utrorf$suf.gz$/, @okd) : (); 
				if($okin and $altin and not -f $fout) {
					$ok= open($hout,'>',$fout); 
					($ok,$hin)= openRead($okin);  while(<$hin>){ print $hout $_; } close($hin);
					($ok,$hin)= openRead($altin); while(<$hin>){ print $hout $_; } close($hin);
					if($utrin) { ($ok,$hin)= openRead($utrin); while(<$hin>){ print $hout $_; } close($hin); }
					close($hout);
					}
				}
    }
    
  }
  return($cdnaseq);  # , $aaseq, $cdsseq    
}

  
## see tr2aacds.pl:asmdupclass_fileset
## moved to cdna_evigenesub.pm, expecting getmRNA:ALSOMAKE_AACDS
sub update_mrna_fileset_OLD
{
  my($trpath, $inmrna, $ntrim, $trimids, @trimfiles)= @_; # @trimfiles = $trimmrna, $trimaa, $trimcds + more?
  my $upstatus=0;
  # outputs now should go to pubdir, but use inmrna path.
  my($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles; # hash instead? for below %fset

## ........ **#@&@!*% this file name wrangling is a big waste of time ..........
## ........ use fewer naming choices ???  no okdir for pubdir set
  
  # cdnaseq = input .mrna; trimmrna, outuncut = updated .mrna
  # FIXME: drop outuncut creation, dont need w/ this update fileset,
  #   instead create idfile of cuts == logfile? :  oid, newstats,  oldstats
  # may be suf.gz instead .. should makename() check this?
  #FIXME: okboth.aa, .cds may not exist where okboth.mrna is input;
  #  .. may exist as .okay.aa, .okalt.aa parts
  # see  evigene/scripts/evgmrna2tsa.pl:getmRNA() : should this be here?
  
  my $flagtrimvec= makename($inmrna,".trimvec_done"); 
  my $aaseq = makename($inmrna,".aa"); ## using this below as fname patt
  my $cdsseq= makename($inmrna,".cds"); 

  my($okaa,$altaa,$okcds,$altcds,$okdir,$pubdir);
  (my $mrnapath= $inmrna) =~ s,/[^/]+$,,; ## THIS MAY BE WRONG NOW .. publicset/ vs okayset/
  my $okpath="$trpath/okayset";
  ($pubdir)= getFileset($mrnapath); # getOkFileset
  ($okdir) = getFileset($okpath);

## FIXME: BADDD grep file ****
#egr: uvcut: nin=75302, nochange=58599, ncut=16672, ndrop=20, nbadcut=11 to publicset/locust1all5asm.uvcut.mrna
#egr: update_fileset.aa upd=16683, same=0,
#   ./okayset/locust1all5asm.aa.ids,0 + publicset/locust1all5asm.uvcut.aa 
#  >./okayset/locust1all5asm.aa.ids <<< WRONG FILE name
#egr: update_fileset.cds upd=16683, same=59013, ./okayset/locust1all5asm.okay.cds.gz,./okayset/locust1all5asm.okalt.cds.gz + publicset/locust1all5asm.uvcut.cds >publicset/locust1all5asm.cds
#egr: update_fileset.mrna upd=16683, same=58599, publicset/locust1all5asm.mrna,0 + publicset/locust1all5asm.uvcut.mrna >publicset/locust1all5asm.mrna
  
  ## drop okalt here? see above getmRNA/ALSOMAKE_AACDS, expect/require pubdir/aa,cds?
  my $aapatt=basename($aaseq);
  my($aaseq1) = grep /$aapatt$|$aapatt.gz$/, (@$pubdir); #NO: ,@$okdir);  # drop aaset path for grep
  ($okaa) = grep /.okay\.aa$|.okay\.aa.gz$/, @$okdir;  
  ($altaa)= grep /.okalt\.aa$|.okalt\.aa.gz$/, @$okdir; 
  if($aaseq1) {
    $aaseq=$aaseq1; $okaa= $aaseq; $altaa=0;
  } elsif($okaa and $altaa) {
    #done above# $aaseq = makename($inmrna,".aa"); # failed empty aaseq from above
  } else { } # fail?
  
  my $cdspatt=basename($cdsseq);
  my($cdsseq1) = grep /$cdspatt$|$cdspatt.gz$/, (@$pubdir); # NO: ,@$okdir);  
  ($okcds) = grep /.okay\.cds$|.okay\.cds.gz$/, @$okdir;  
  ($altcds)= grep /.okalt\.cds$|.okalt\.cds.gz$/, @$okdir; 
  if($cdsseq1) {
    $cdsseq=$cdsseq1; $okcds= $cdsseq; $altcds=0;
  } elsif($okcds and $altcds) {
    #done above# $cdsseq= makename($inmrna,".cds");
  } else { } # fail?
  
  return ( 0, $inmrna, $aaseq, $cdsseq) if( -f $flagtrimvec); #check for flag-file made below  
  return (-1, $inmrna, $aaseq, $cdsseq) if($ntrim < 1 or ! -s $trimmrna);

  my $upmrna  = makename($inmrna,".mrna_upd"); 
  my $upaaseq = makename($inmrna,".aa_upd"); # failed empty aaseq from above
  my $upcdsseq= makename($inmrna,".cds_upd");  # failed empty cdsseq from above
  
  $aaseq =~ s/.gz$//; $cdsseq =~ s/.gz$//; # output fupname
  (my $outmrna= $inmrna)=~ s/.gz$//;
  my %fset= (
            # $nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname
    mrna => [ 0, 0, $inmrna, 0, $trimmrna, $upmrna, $outmrna ],
    aa   => [ 0, 0, $okaa, $altaa, $trimaa, $upaaseq, $aaseq],  #?? fix here for .okay.aa + .okalt.aa ?
    cds  => [ 0, 0, $okcds, $altcds, $trimcds, $upcdsseq, $cdsseq],
  );
  
  ## need ids from mrna_trimvec
  ## open input old & newtrim files, open out newmerge files
  ## foreach input fasta id, print if old & not trimid{id}, print if new & trimids->{id}

  ## FIXMEd: fails on aaseq, cdsseq; outmrna ok : open(.gz) bug
  ## FIXME-maybe: old aahdr has evgclass tags:  evgclass=main,okay,match:locust1sop4p4k23loc30t48,pct:99/96;
  ## .. copy to new??  in update_mrna_fileset ?
  
  foreach my $suf (sort keys %fset) { # is inmrna always .mrna ?
    my($ok,$hin,$hup,%keptids);
    my($nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
    if(-s $fin and -s $ftrim) {
    	## pull okayset/{okay,okalt}.$suf, skipping trimids{id}
      $ok= open($hup,'>',$fup); # unless($ok)...
      ($ok,$hin)= openRead($fin); # $ok= open($hin,$fin); # FIX: may be .gz 
      $ok=0; while(<$hin>){ 
      	if(/^>(\S+)/) { my $d=$1; $ok=($trimids->{$d})?0:1;  do{ $keptids{$d}++; $nsame++;} if $ok;} 
      	print $hup $_ if($ok); }
      close($hin);
      if($fin2) { 
      ($ok,$hin)= openRead($fin2); #$ok= open($hin,$fin2); # unless($ok)...
      $ok=0; while(<$hin>){ 
      	if(/^>(\S+)/) { my $d=$1; $ok=($trimids->{$d})?0:1; do{ $keptids{$d}++; $nsame++;} if $ok;} 
      	print $hup $_ if($ok); }
      close($hin);
      }
      
    	## pull trimset/uvcut.$suf, check? for trimids{id} and/or collect above kept ids
      ($ok,$hin)= openRead($ftrim); # $ok= open($hin,$ftrim); 
      $ok=0; while(<$hin>){ 
      	if(/^>(\S+)/) {  my $d=$1; $ok=($trimids->{$d} and not $keptids{$d})?1:0; $nup++; } 
      	print $hup $_; } 
      close($hin);
      close($hup);
      $fset{$suf}->[0]= $nup; $fset{$suf}->[1]= $nsame;  # $fset{$suf}= [$nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname];
      $upstatus++ if($nup>0);
    } else {
      ## error, warn/loggit, skip?
      loggit(1, "ERR update_fileset.$suf empty $fin or $ftrim"); 
    } 
  }
  ## end merge loop

  my @outfiles;
  if($upstatus == 3)  { # $upstatus == 3 or > 0?
    ## rename input files to input.old, better: input.untrim .notrim? .pretrim?
    ## rename newmerge files to input
    foreach my $suf (sort keys %fset) { 
      my($nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      if(-s $fupname) { rename($fupname,"$fupname.untrim"); push @tmpfiles, "$fupname.untrim"; }
      rename($fup,$fupname); push @outfiles, $fupname;
      push @tmpfiles, $ftrim; # keep for checking
      loggit(0, "update_fileset.$suf upd=$nup, same=$nsame, $fin,$fin2 + $ftrim >$fupname"); 
    }
    runcmd("touch $flagtrimvec"); ## touch flag-file that new input has uvtrim results ..
  } else {
    loggit(1, "ERR update_fileset missing status=$upstatus/3");  # list fupnames?
    foreach my $suf (sort keys %fset) { # sort: aa,cds,mrna
      my($nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      push @outfiles, $fup;
      loggit(1, "ERR update_fileset.$suf upd=$nup, same=$nsame, $fin,$fin2 + $ftrim >$fup"); 
    }
  }
  return ($upstatus, @outfiles); # $upaaseq, $upcdsseq,$upmrna, 
  ## return what?
}


sub mrna_trimvec 
{  
  my($mrnaf,$vtab,$genenamef)= @_;
  my %trimids=();
  
  unless ( -f $mrnaf ) { loggit(1, "ERR: uvcut missing mRNA $mrnaf\n"); return -1; }  
  unless ( -f $vtab  ) { loggit(1, "ERR: uvcut missing $VECSUF $vtab\n"); return -1; } 
  #FIXME: keep & return id list of cut,uncut for update_mrna_fileset

  my($hin, $houtf, $houtaa, $houtcds, $hidlist) = (undef) x 9; # file handles replace hard names
  #was IN OUT AAOUT CDSOUT OUTUNCUT, $houtuncut
  my $nam= makename($mrnaf,""); # (my $nam = $mrnaf) =~ s/\.mrna.*//; 
  ## (my $cutsuf= $outsuf) =~ s/\.mrna//; ## or just use outsuf?
  my($outf,$outaa,$outcds,$outids)= map{ "$nam.$outsuf.$_" } ("mrna", "aa", "cds", "ids");
  ## change these to hash of handle/filenames by suffix keys= ("mrna", "aa", "cds", "ids")
  my @outf=($outf, $outaa, $outcds, $outids);
  my @outh=($houtf, $houtaa, $houtcds, $hidlist);
  loggit(0, "uvcut $mrnaf with $vtab to $outf\n"); 

#   my($sufaa,$sufcds,$sufids,$sufuncut); # dont need all these vars
#   $sufaa=$outsuf; unless($sufaa=~s/\.mrna/.aa/){ $sufaa.=".aa"; } # dont need globals here
#   $sufcds=$outsuf; unless($sufcds=~s/\.mrna/.cds/){ $sufcds.=".cds"; }
#   $sufids=$outsuf; unless($sufids=~s/\.mrna/.ids/){ $sufids.=".ids"; }
#   $sufuncut=$outsuf; unless($sufuncut=~s/cut/uncut/){ $sufuncut.=".uncut"; }
#   $outf="$nam.$outsuf"; $outaa="$nam.$sufaa"; $outcds="$nam.$sufcds"; $outids="$nam.sufids"; 
#   #DROP# $outuncut="$nam.$sufuncut";
  
  my($namgot,$namin)= parse_genenames($genenamef);  # do in caller?
  loggit(0, "uvcut names $genenamef n=$namgot\n"); 
  
  my $nvecid= readVectab($vtab);  
  loggit(0, "uvcut tr with vector n=$nvecid\n"); 

  # NO/FIXME: add OUTUNCUT for non-uvector set, later combined
  # add gapclean() / trimNNNends() ? from evgmrna2tsa.pl:putseq()
  # .. needs to adjust CDSoffset, maybe CDS if partial with NNN in ENDGAP; easier to recalc prot as w/ uvcut
  
  my $ok=0; 
  # if( $mrnaf =~ /\.gz$/ ) { $ok= open($hin,"gunzip -c $mrnaf |"); } else { $ok= open($hin,$mrnaf); }
  ($ok,$hin)= openRead($mrnaf);
  if($ok) { 
    for(my $i=0; $i<@outf; $i++) { $ok= open($outh[$i],'>',$outf[$i]); last unless($ok); } 
    ($houtf, $houtaa, $houtcds, $hidlist)= @outh; # for below... fixme
    # $ok= open($houtf,'>',$outf); $ok= open($houtaa,'>',$outaa); $ok= open($houtcds,'>',$outcds); $ok= open($hidlist,'>',$outids); 
    }   # DROP# $ok= open($houtuncut,'>',$outuncut); 
  unless($ok) { loggit(1,"ERR: bad files in:$mrnaf out:$outf,.aa,.cds,.."); return -1; }

## FIXME TOO MANY 'uvcut=nocut' in .uvcut, from hasBadGaps() == 8 == has $GAPSMAX, now allowed...
## FIXME1: 1 case of single 'N' at end5 now is skipped. this came from utrorf.
## SEQ_INST.TerminalNs  N at beginning of sequence: RhipulEGm009616t2 

  my($id,$fa,$hd,$hasVecOrNNN,$didput,$nmin,$nput,$nskip,$nbadcut,$nnochange)=(0)x9;
  while(<$hin>) {
    if(/^>(\S+)/) { my $d=$1; $nmin++;
      if($id) { 
        $hasVecOrNNN |= hasBadGaps($fa);
        ## PROBLEM skipping HasGapMaxSpan : CDShasTooManyXs or not GAPSOK
        unless($hasVecOrNNN == HasGapNone) {  ##  or $hasVecOrNNN == HasGapMaxSpan
          $didput= putVecOrNNN($houtf, $houtaa, $houtcds, $hidlist, # \@outh,
                $id, $hd,$fa,$hasVecOrNNN); 
          if($didput == -99) {  
            $nnochange++; # special means $nnochange++
          } else {
            $trimids{$id}=1; # regardless of $didput return
            $nput++ if($didput>0); $nbadcut++ if($didput<0); $nskip++ if($didput==0); 
          }
        } else { 
          ##note: $hd,$fa have \n endlines here
          $nnochange++; #DROPuncut# print $houtuncut $hd,$fa; 
        }
      } 
      $hasVecOrNNN= ($vecscreen{$d})? HasGapVector: HasGapNone; 
      $id=$d; $hd=$_; $fa=""; 
      } 
    elsif(/\w/) { $fa .= $_; } # dont chomp now, for OUTUNCUT/ASIS
  } 
  
  if($id) { 
    $hasVecOrNNN |= hasBadGaps($fa);
    unless($hasVecOrNNN == HasGapNone) {  ##  or $hasVecOrNNN == HasGapMaxSpan
      $didput= putVecOrNNN($houtf, $houtaa, $houtcds, $hidlist,  # \@outh,
                $id,$hd,$fa,$hasVecOrNNN); ## add outfile handles as param? keep local?
      if($didput == -99) { 
        $nnochange++; # special means $nnochange++
      } else {
        $trimids{$id}=1; # regardless of $didput return
        $nput++ if($didput>0); $nbadcut++ if($didput<0); $nskip++ if($didput==0); 
      }
    } else { 
      $nnochange++; #DROPuncut# print $houtuncut $hd,$fa; 
    }
  } 

  close($hin); for(my $i=0; $i<@outh; $i++) { $ok= close($outh[$i]); }
  # close($houtf); close($houtaa); close($houtcds); close($hidlist); 
  #DROPuncut# close($houtuncut);
  #?? write %trimids to file? should be in .log .. maybe not.

	##add this: 
	if($nbadcut>0) { (my $po=$outids)=~s/ids/PROBLEMCUT/; runcmd("grep PROBLEMCUT $outids > $po"); push @outf, $po; }
  
  my $nalltrim= $nput + $nskip + $nbadcut; # badcut also counts
  my $trimlog="nin=$nmin, nochange=$nnochange, ncut=$nput, ndrop=$nskip, nbadcut=$nbadcut to $outf";
  loggit(0, "uvcut: $trimlog"); 
  return($nalltrim, \%trimids, @outf); ## $outf, $outaa, $outcds, $outids); # , $outuncut
} # end sub mrna_trimvec(mrna)


sub hasBadGaps {
  my($fa)= @_;
  $fa =~ s/\s+//g; # do this once not each use?
  my $clen=length($fa);
  my $nbig= index($fa,$GAPSMAX); 
  my $n1= index($fa,'N'); 
  my $ne= rindex($fa,'N'); 
  my $hasNNN=HasGapNone; # 0
  $hasNNN |= HasGapEnd5 if($n1 >= 0 and $n1 < $ENDGAP);
  $hasNNN |= HasGapEnd3 if($ne >= $clen - $ENDGAP);
  $hasNNN |= HasGapMaxSpan if($nbig >= 0);
  return $hasNNN;
}

sub readVectab
{
  my($vtab)= @_;
  %vecscreen=(); # global now
  my $nvecid= 0;
  my($ok,$hin)= openRead($vtab); #open(F,$vtab) 
  if($ok) {
    while(<$hin>){ my($id,$b,$e,$vt)=split"\t"; $vecscreen{$id} .= "$b\t$e\t$vt\n"; } close($hin);
    ## compress overlaps here 
    foreach my $oid (keys %vecscreen) {
      my @vec=split"\n", $vecscreen{$oid};  
      if(@vec>1) { 
        @vec= sort { $a <=> $b } @vec; 
        my @vec2=(); my $ncut=0;
        for(my $i=@vec - 1; $i>0; $i--) { 
           my($ub,$ue,$vty)=split"\t",$vec[$i-1];
           my($xb,$xe,$xty)=split"\t",$vec[$i]; 
           next if($xb<1 or $xe < $xb or $xty =~ /Weak/i); # bad data?
           if($xb <= $ue) { $ub=$xb if($xb<$ub); $vec[$i-1]=join("\t",$ub,$xe,$vty); $ncut++; }  
           else { unshift @vec2, $vec[$i]; }
           }
        unshift @vec2, $vec[0];  $vecscreen{$oid}= join("\n",@vec2);
      }  
    }
  }
  # return %vecscreen;
  $nvecid= scalar(keys %vecscreen);
  return $nvecid;
}

 
sub putVecOrNNN { 
  my($houtf, $houtaa, $houtcds, $houtidlist,
     $oid,$hdr,$fain,$hasVecOrNNN)= @_;
  ## param of outhandles: $houtf, $houtaa, $houtcds  for OUT AAOUT CDSOUT
  #now: hasVecOrNNN & 1 == uvec; & 2 == end5gap; & 4 == end3gap; & 8 == biggap
  ##c{ HasGapNone => 0, HasGapVector => 1, HasGapEnd5 => 2, HasGapEnd3 => 4, HasGapMaxSpan => 8, }; # hasBadGaps

  my $retval=0;
  $fain =~ s/\s+//g;
  my $olen=length($fain);
  
  my $tblinfo= parse_evgheader($oid,$hdr,$olen);
  my $pubid= $tblinfo->{'pubid'};
  my $cdsoff= $tblinfo->{'cdsoff'}; 
  my ($cdsb,$cdse)= split/[-]/,$cdsoff; 
  my $oldaaq= $tblinfo->{'aaqual'};    
  my($aafull)= $oldaaq =~ m/(complete|partial\w*)/; 
  my $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1;
  my $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;
  my $oldorflen = 1 + (($cdsb>$cdse)? $cdsb - $cdse : $cdse - $cdsb);
	my $oldcdsnn=0;
  my $gnamed= $tblinfo->{'name'} || $genenames{$oid} || 0;
  my $namepct= $genenamepct{$oid} || 0; # is structured: 99%,123/345,678
  my $nameref= $genedbxref{$oid} || 0; #  
  my ($okname,$uniqname)= (0,0);
  if($gnamed) {  
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    $okname= (not $gnamed or $gnamed =~ /^CDD:|^hypothetical|^uncharacterized|^na$/i or $npct<50) ? 0 : 1;
    my $allids= $namedgenes{$gnamed}||"na"; 
    $uniqname=1 if($okname and $allids eq "$oid,");
  }
        
  
  my($ucut,$rue,$rub)= (0) x 9;
  my($fac,$facgap,$fav,$vectype)=("") x 10;
  
  if($vecscreen{$oid}) {
    # FIXME: need special case for vector cutting start-codon, insert new ATG (or part missing) in facgap
    ($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype)= cutVector($oid,$fain,$cdsb,$cdse); # if( $hasVecOrNNN & 1 > 0);
    if($ucut > 0.5*$olen and not $fac) { $hasVecOrNNN = HasGapVector; } #cut ALL ; dont test NNN ; dont change to fain
  }
  
  # FIXME: bug DROPS are not being dropped, or kept no change by mistake:
  ## ^^ Problem here, fac becomes blank, but that test below changes to fain : 
  #   sowhiteflyv8k61loc90865t1  UVector! really big uvcut,
  #   uvcut4b #BAD neworf uvtype=Strong match, sowhiteflyv8k61loc90865t1: loss orflen=0-1176; loss named=92%,293/320,391,CDD
  #   uvcut5f sowhiteflyv8k61loc90865t1 uvcut=Strong,1736,1-1736,cdsinend3; clen=1736/1736; << NO CUT

  # FIXME2: lots of ENDGAP in (partial) CDS, but have valid cds/prot before/after those endgaps
  # .. instead of trimming to end and chopping valid cds-bases, should squeeze these gaps down to 3+frame minimum
  ## vectype="endtrim" for non-UVector but endtrim; ..
  ## dang, do we need to trim both fac, facgap ?
  ## FIXME: which uvec vals need endtrim update? ucut? rub,rue?
  my $trimtype=""; my $ntrim= 0; my $nNcut=0; 
  
  if( $hasVecOrNNN > HasGapVector) {   
    my $fadegap= $fac || $fain; # NOT for fac == cutVector entirely
    my $gappy= hasBadGaps($fadegap); # check again vector may have cut
    unless($gappy == HasGapNone) # NO!! or ($GAPSOK && $gappy == HasGapMaxSpan)
    {
      my ($trimgap); 
      ## change this to not endTrim in CDS ! is chopping too many valid bases of partials, but need new cds after vectrim
      ## eg: litovavel3k55Loc3448t5,cdsoffs=2-400: 1-=TGAAATCGTCCGTNNNNNNNNNNNGCTTTGCTA trim>GCTTTGCTACA
      my @cdsvalid= ($fac)? () : ($cdsb,$cdse);
      ($fac,$trimtype,$ntrim,$nNcut,$oldcdsnn)= endTrimNNN($oid,$fadegap,@cdsvalid);
      ($facgap,$trimgap)= endTrimNNN($oid,$facgap) if($facgap);
      $vectype.=$trimtype if($trimtype);
      ## what when ntrim==0 && trimtype == "" ?? nocut
      # $ucut += $ntrim; # below
    }
  }
  ## FIXME: badcut should not count cds-nnn-squeezes, inframe reduction of nnn shouldnt affect prot aligns
  ## .. subtract $nNcut  from tests for bad changes = num NNN cut
  
  my $vecNotStrong = ($vectype =~ /Strong/i)?0:1;
  my $vecNotSqueeze= ($vectype eq "cdsns")?0:1; # special trimtype for cds gaps
  
  ## fixme .. old: rub,rue
  my $clen=length($fac); my $fl=""; 
  if($clen == $olen and not $vectype) { $vectype = "nocut"; } # vectype="nocut" # nochange? uncut?
  elsif(not $vectype) { $vectype = "errcut"; } # what? err?
  
  #BAD# $fac= $fac || $fain; # # NOT for fac == cutVector entirely
  unless($fac =~ /\w/) {  # skip to dropit?? $fac ||= $fain; is this ever right?
    $fl.="allcut"; $clen=0; # clen should == 0
    if( $ucut+$ntrim < 0.5*$olen) { } #problem
  }

  if($rub < $cdse and $rue > $cdsb) { if($rub <= $cdsb+2 and $rue<$cdse) { $fl.="cds5"; }
  elsif($rue >= $cdse-2 and $rub > $cdsb) { $fl.="cds3"; } else { $fl.="cdsin"; } }
  if($rue > $olen-9) { $fl.="end3"; } elsif($rub <= 9) { $fl.="end5"; }
  
  #FIXME2: need bestorf -partial option (always?) so that uvcut doesnt trigger further chop for complete-aa
  #FIXME: here? add protein qual check: cdna_bestorf($fac) : 
  ## if( aasize << cutsize/3? AND vectype = Moderate? and aahomol > minalign) 
  ##   then cancel cutvec; keep orig w/ note; check homol-align also?
  
  my $note=""; my $badcut=0;
  my $expect_neworflen= $oldorflen - $ucut; # not quite right for ucut in UTR
  my ($newaahdr,$newaa,$newcdshdr,$newcds) = ("") x 10;
  my ($newaalen,$newpcds,$newcompl,$neworflen)= (0) x 10; 
  
  if($clen > 0) {
  ($newaahdr,$newaa,$newcdshdr,$newcds, $newaalen,$newpcds,$newcompl,$neworflen)
      = getbestorf($oid,$hdr,$fac,$clen, $expect_neworflen);

  if($neworflen < $oldorflen) {
    my $cglen= length($facgap);
    # warn "# fagap.$oid=$cglen,$facgap\n" if $DEBUG;

    if($cglen > $clen) {
    my ($newaahdr1,$newaa1,$newcdshdr1,$newcds1,
       $newaalen1,$newpcds1,$newcompl1,$neworflen1)= getbestorf($oid,$hdr,$facgap,$cglen, $expect_neworflen);
    my $nga= $newaa1 =~ tr/X/X/; # my $ngc= $newcds1 =~ tr/N/N/;   
    my $newleng= $newaalen1 - $nga;
    # warn "# aacut=$newaahdr; lens=$newaalen1-$newaalen; aacgap.$newaahdr1\n" if $DEBUG;

    if( $newleng > $newaalen) { 
      $note .= "framefixlen=$newleng-$newaalen; ";
      ($newaahdr,$newaa,$newcdshdr,$newcds, $newaalen,$newpcds,$newcompl,$neworflen) = 
        ($newaahdr1,$newaa1,$newcdshdr1,$newcds1, $newaalen1,$newpcds1,$newcompl1,$neworflen1);
      ($fac,$clen)= ($facgap,$cglen);   
      } 
    }
  }
  ## FIXME-maybe: old aahdr has evgclass tags:  evgclass=main,okay,match:locust1sop4p4k23loc30t48,pct:99/96;
  ## .. copy to new??  in update_mrna_fileset ?
  } # fac/clen
  
  ## ?? here, do after recompute cds  
  my $CDShasTooManyXs=0; my $newcdsnn=0;
  if(1) {
    my($newcdsb,$newcdse)= $newaahdr =~ m/offs=(\d+).(\d+)/; 
    my $cdsfa= substr($fac,$newcdsb-1,1+$newcdse-$newcdsb); ## isnt this $newcds == formatted \n 
    my $cdsw= length($cdsfa);
 		$newcdsnn= $cdsfa =~ tr/N/N/; 
    if($cdsw > 1 and $newcdsnn > 0.48*$cdsw) { # ERR
      $CDShasTooManyXs= $newcdsnn; # flag it for below .. uniqname ?? check below
      $fl .= ",CDShasTooManyXs:$newcdsnn/$cdsw";
      $hasVecOrNNN |= HasGapTooManyXs;
    }
  }
  
  my $cdsnncut= ($oldcdsnn > 0) ? $oldcdsnn - $newcdsnn : 0; # use to adjust badcut tests
	my $neworflentest=$neworflen;  $neworflentest += $cdsnncut if($cdsnncut > 2);
  
  #FIXME here for no cut-able gaps?
  if($hasVecOrNNN == HasGapNone or ($GAPSOK && $hasVecOrNNN == HasGapMaxSpan)) { 
    $retval= -99;     # return and print asis
    return($retval);
  } 
  
  loggit(0, "getbestorf: $oid oldaa=$oldaaq; new.$newaahdr\n"); 
        
  if($neworflen > $oldorflen) { # DEBUG check these
    # my $namepct= $genenamepct{$oid} || 0; # is structured: 99%,123/345,678
    # my $nameref= $genedbxref{$oid} || 0; #  
    # my $named= $genenames{$oid} || 0; #  same as $gnamed
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    my $odiff= $neworflen-$oldorflen;
    loggit(0, "named:  $oid LONGER dorf=+$odiff, named=$namepct,$nameref,$gnamed;\n") if($npct>=50 or $nameref =~ /:/);  
  }
  
  if($neworflentest < $oldorflen - 9) { # check named align; also debug check names for neworflen >> oldorflen ?
    #  %genenames=%genenamepct=%genedbxref
    # my $namepct= $genenamepct{$oid} || 0; # is structured: 99%,123/345,678
    # my $nameref= $genedbxref{$oid} || 0; #  
    # my $named= $genenames{$oid} || 0; #  
    my $odiff= $neworflentest-$oldorflen;
    
    ## ?? ignore 'CDD:' as main name, it lacks related species homolog
    ## .. not sure ignore CDD: is right yet.
    # SEE ABOVE: my $okname= 1; 
    # my $okname= (not $gnamed or $gnamed =~ /^CDD:|^hypothetical|^uncharacterized|^na$/i) ? 0 : 1;
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    loggit(0, "named:  $oid dorf=$odiff, named=$namepct,$nameref,$gnamed;\n") if($npct>=50 or $nameref =~ /:/); 
    
    my($minpct,$minaln)= ($vecNotStrong) ? (66,50) : (90,90);
    ## for vecIsStrong, want to measure size of uvcut vs size of align? really need to know if uvcut is align-supported
    ## but dont yet have ref-align spans as inputs
    
    #FIXME: adjust badcut to accept Strong when also partial-cds at cut end?)
    # .. eg: shrimpt nbadcut=26 but 20 are uv=Strong, cut ~30 cds; probably accept all Strong uvcut
    # .. also have some perfect huge vector matches w/ complete prot: NUBP1 = Ecoli vec
    if($okname and not $vecNotStrong) {
      if($oldaaq =~ /partial/ and $neworflentest >= 0.95 * $expect_neworflen) { $okname=0; }
      elsif($oldaaq !~ /partial/  and $ucut > 99) { $okname=0; } # this is a long-Strong uvec, named likely vector protein
    }    
    $okname=0 if(not $vecNotSqueeze); # if($vecIsSqueeze); # not badcut for these, I hope..
     
    ## ?? need levels of npct test here, dont know unless npct ~100 if uvcut affects alignment
    
    if($okname and $npct >= $minpct and $naln >= $minaln and $neworflentest/3 < 0.99*$naln) { ##  ??
      $fl .= ",refalignlosscut:$odiff";
      $badcut++; $note.="loss named=$namepct,$nameref,$gnamed; "; # d=$odiff, 
    }
  }

  ## Disallow  $vectype =~ /Strong/ here?   ; skip this if refalignlosscut ?
  ## FIXME: add note for any largish change in orflen, smaller/bigger, whether uvec or trimNNN
  if(not $badcut and ($newaalen<1 or $neworflentest < 0.90 * $expect_neworflen)) {
    my $odiff= $neworflentest-$oldorflen; # neworflentest or neworflen ??
    $fl .= ",shortorfcut:$odiff"; # note even for Strong? yes. 
    if($vecNotStrong and $vecNotSqueeze) { $badcut++; $note.="loss orflen=$neworflentest-$oldorflen; "; }
  }
  
  $ucut += $ntrim; #?
  my $idflag="";  
  my $uvhdr="uvcut=$vectype,$ucut,$rub-$rue,$fl;";
  if($badcut) {  
    $idflag.="PROBLEMCUT,";
    # complain, maybe cancel uvcut; not badcut if $vectype=~/Strong/ always?
    loggit(1, "BAD neworf $uvhdr $oid: $note\n");
    $retval= -1;
  }
  
  chomp($hdr); 
  $hdr=~s,clen=,clen=$clen/,; 
  ## update mRNA hdr from $newaahdr : 
  
  if(length($fav)>45) { $fav=substr($fav,0,20)."..".substr($fav,-20); } ## dont stick LONG fav in header ?
  my $uvh2=$uvhdr; $uvh2.=" uvfa=$fav;" if($fav);
  $hdr=~s/ / $uvh2 /; 
  $fac =~ s/(.{60})/$1\n/g; $fac.="\n" unless($fac=~/\n$/);
  
  ## MINSIZE should change w/ named value, as per evgmrna2tsa.pl:putseq() 
  ## unique(name)/strong homol-name should keep shorter clen; but w/ annotation to support short len
  ## FIXME: dropit conflict with uniqnamed genes. from mrna2tsa

  my $dropit= ($clen<$MINSIZE or $CDShasTooManyXs)?1:0;  
  if($dropit and $uniqname) { 
    if($CDShasTooManyXs) { } # drop anyway, uniqname is suspect?
    elsif(($newaalen > 15 and $newcompl =~ /complete/) or $newaalen > 20) { $dropit=0; }
    my $dval= ($dropit)? "DROP.":"KEEP.";
    $dval .= ($CDShasTooManyXs) ? "gaps:$CDShasTooManyXs" : "short:$clen";
    $idflag.="PROBLEMCUT.uniquename:$gnamed,"; $idflag.="error:$dval," unless($dropit);
    loggit(1,"PROBLEM: unique name '$gnamed' but error:$dval"); 
  }
  
  if($dropit) { 
    my $dval= ($CDShasTooManyXs) ? "gaps:$CDShasTooManyXs" : "short:$clen";
    $idflag .= "DROPCUT.$dval,";
    loggit(1, "DROPCUT.$dval: $hdr\n"); 
    $retval= ($badcut)? -1 : 0;  # change from 0 to ?? -2
  } else { 
    $idflag.="OKCUT," unless($idflag =~ /CUT/);
    # $houtf, $houtaa, $houtcds, .. change to hash of handles?  $houts{mrna}, $houts{aa} ..
    print $houtf   $hdr,"\n",$fac;  
    print $houtaa  ">$oid $newaahdr; $uvhdr\n$newaa\n";    
    print $houtcds ">$oid $newcdshdr; $uvhdr\n$newcds\n";   
    $retval= ($badcut)? -1 : 1;
  }
  print $houtidlist  join("\t",$oid,$idflag,$uvhdr,$newaahdr)."\n"; # if($houtidlist);
  return $retval;
}    


=item UVCUT > ORFloss too big: bad bestorf
    
# ** Problem w/ getbestorf() only gets new starts at ATG/M, but cut over start can damage this
# .. need getbestorf() to look for partial5 after each stop codon.. patch here when cut thru ATG startcodon
# .. should add back that plus cds-inframe gap:  s/$uvectorwithstart/UtrnnnATGnnCds/; check 'cds5' flags for cdsloss > uvcut

grep BAD log.uvcut5f | perl -ne'($uvc)=m/uvcut=\w+,(\d+)/; ($oc,$ob)=m/orflen=(\d+).(\d+)/; $oloss
1=$ob-$oc; ($oloss)=m/cut:-(\d+);/; $oloss||=$oloss1; print if($oloss > 1.1*$uvc); ' | head
#BAD neworf uvcut=Moderate,39,144-182,cds5,refalignlosscut:-810; socatfishv1k95loc34193t1: loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 
#BAD neworf uvcut=Moderate,50,134-1779,cds3end3,shortorfcut:-219; socatfishv1k25loc6796t5: loss orflen=1356-1575; 

.. not too many cases ..
/bio/bio-grid/aabugs4/tsaevgc/tsaoutz

cacao3all7f/log.uvcut5f
#BAD neworf uvcut=Moderate,26,429-454,cdsin,shortorfcut:-120; cacao3vel14sc9k35Loc1726t4: loss orflen=489-609; 
  >> has cdsgaps, problem likely is those gaps, not cut but mangle new orf call, same prot start.
#BAD neworf uvcut=Moderate,26,695-720,cdsin,shortorfcut:-105; cacao3sopcsc2k29loc6700t1: loss orflen=522-627; 
  >> ditto, NNNN cdsgaps trail uvcut, mangle new orf call
  
catfish1all4cf/log.uvcut5f
#BAD neworf uvcut=Moderate,39,144-182,cds5,refalignlosscut:-810; socatfishv1k95loc34193t1: loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 
  below>> DEFINITELY cancel uvcut socatfishv1k95loc34193t1 (or replace ATG start in cutspan)
#BAD neworf uvcut=Moderate,50,134-1779,cds3end3,shortorfcut:-219; socatfishv1k25loc6796t5: loss orflen=1356-1575; 
  below>> SHOULD ok this uvcut; >> 2 uvcuts, damage cds, which also has NNN spans;

whitefly1evgcf/log.uvcut5f
#BAD neworf uvcut=Moderate,25,308-332,cds5,shortorfcut:-207; whitefly1vel6k39Loc14435t1: framefixlen=40-24; loss orflen=129-336; 
  >> uvcut is in startcodon-base3, replace 'G' corrects this overcut
  orig>  aalen=111,37%,complete-utrpoor; clen=906; strand=+; offs=306-641;
  G+cut> aalen=104,35%,complete-utrpoor; clen=885; strand=+; offs=306-620;  same prot
  G-cut> aalen=42,14%,complete-utrbad; clen=885; strand=+; offs=259-387;    diff prot
  
banana1all3cf/log.uvcut5f
litova1all3f/log.uvcut5f
locust1evgcf/log.uvcut5f
pogonus1all3cf/log.uvcut5f
shrimpt1evgf/log.uvcut5f
zticktr2acf/log.uvcut5f

=cut

sub squeezeNNN {
  my($fa,$nlower,$gapsmax,$keepnnn,$inmax)= @_;
  my $ncut=0;
  my $gapw= length( $gapsmax);
  # keepnnn = 0 or 3 for squeeze keep
  my $NNNs= 'NNNNNNNNNNN'; my $N1='N';
  if($nlower) { map{ $_= lc($_) } ($NNNs,$N1,$gapsmax); } # NOT USED here
  my $dorev=0; if($inmax < 0) { $fa=reverse($fa); $dorev=1; $inmax= -$inmax; }

  for (my $in= index($fa,$gapsmax); $in >= 0; $in=index($fa,$gapsmax)) {
    last if($inmax>0 and $in>=$inmax);
    my $w=length($fa); my $en=$in+$gapw; 
    $en++ while($en<$w and substr($fa,$en,1) eq $N1); 
    my $wn= $en-$in; 
    my $keep= $keepnnn + ($wn % 3); # PROBLEM: got +10 nnn in cases; keep shold not be that big
    $keep=3 if($keep==0);
    my $cut= $wn-$keep; $ncut+=$cut; 
    my $facut= substr($fa,0,$in).substr($NNNs,0,$keep).substr($fa,$en); 
    $fa=$facut; 
  } 
  if($dorev) { $fa=reverse($fa); }
  return($fa,$ncut);
}

#    ($fain2,$trimtype)= endTrimNNN($oid,$fain2);
sub endTrimNNN {
  my($oid,$fa,$cdsb,$cdse) = @_;
  my $trimtype="";
  my $ncut=0; 
  my $olen=length($fa);
  $fa =~ s/n/N/g;
  my $nNold= $fa =~ tr/N/N/; my $cdsNold= 0;
  
#  ## drop cds stuff? or not? dont want endtrim over valid cds bases
  my $validcds= (defined $cdse and $cdse>0)?1:0;
  my ($lcdsb,$lcdse)= ($validcds)?($cdsb,$cdse):(0,0);

## 5e:    
  if($validcds) {
    my $fixit=0;
    my $cdsfa= substr($fa,$cdsb-1,1+$cdse-$cdsb); 
    my $cdsw=length($cdsfa);
		$cdsNold= $cdsfa =~ tr/N/N/;
		my $fixone=0;
		
    my $ne= rindex($fa,'NN'); 
    if($ne >= $olen - $ENDGAP and $olen - $cdse < 3*$ENDGAP) {
      $ne= rindex($cdsfa,'NN'); 
      if($ne >=0 and $cdsw - $ne < 2*$ENDGAP) { 
        my $endc= substr($cdsfa,0,-$ENDGAP);
        $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
        $fixone |=2; $fixit |=2 if($zz >= 3);#?? is this zz test bad? bad if fixit 1 or 2 and other zz fails, need fix both
      }
    }

    my $n1= index($fa,'NN'); 
    if( $n1 >= 0 and $n1 <= $ENDGAP and $cdsb < 3*$ENDGAP) {
      $n1= index($cdsfa,'NN'); 
      if($n1 >=0 and $n1 < 2*$ENDGAP) { 
        my $endc= substr($cdsfa,0,$ENDGAP);
        $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
        $fixone |=1; $fixit |=1 if($zz >= 3); #?? is this zz test bad?
      }
    }
    ## zz test is to chop end NNN unless have enough non-N for squeeze.
    if($fixit and $fixone) { $fixit= $fixone; } # ignore zz test if one end passes
    
    if($fixit) { 
      my $cdsfac= $cdsfa;
      my $ccut=0;
      ## fixme: endgap only for 1..end5gap, need reverse(cdsfa) for end3gap
      ## maxgap 4 here, must be bigger than squeezemax of 3
      
      ## BAD here ?? got no SqueezeNNN .. long 'nnnnnnnn' spans in cdsns same as orig mRNA.untrim
      ## loggit(ERR) if ccut == 0
      
      if($fixit & 2 ) {## == 2
        ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,-2*$ENDGAP);
        $ncut += $ccut;
      }
      if($fixit & 1) { ##  == 1
        ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,2*$ENDGAP);
        $ncut += $ccut;
      }
      $cdsfac =~ s/^N+//; $cdsfac =~ s/N+$//;
      $trimtype.="cdsns";
      $cdsfac =~ s/N/n/g; # lower to prevent following trim
      $fa= substr($fa,0,$cdsb-1). $cdsfac . substr($fa,$cdse);
    }
    
  }
        
    ## not right; NN may start in UTR, end in CDS
    # my $incds= ( $n1 >= 0 and $n1 <= $ENDGAP and $n1 >= $cdsb)
    #  or ($ne >= $olen - $ENDGAP and $ne <= $cdse);

    # //no//FIXME4: No NN End trim in CDS for aastart,aastop: keep stop/start codon 
    # FIXME: cdsb,cdse adjust for inner gaps
    # fixme2: must adjust cdsb,e when cut BEFORE cdsb
    # YES: fixme3: this is a mess; better to a. cut NNN, b. rerun cdna_bestorf for new cds offset?
    # FIXME: single end N happens.. need 
    
  $fa=~s/(N+)$//;
  my $curlen= length($fa); 
  my $ne= rindex($fa,'NN');
  for(my $iter=0; $ne >= $curlen - $ENDGAP and $curlen > $ENDGAP and $iter<5; $iter++) {
    $fa= substr($fa,0,$ne); 
    if($fa=~s/(N+)$//) {  my $ncut=length($1); $ne-=$ncut; }
    $curlen= length($fa); 
    $ne= rindex($fa,'NN');  
  }
  $trimtype.="end3trim" if($curlen < $olen);
  
  ## FIXME: cds-phase/codon_start changes w/ mod 3 of n1   
  ## FIX2: see above  nnnAnnn < need to recurse chop endgaps ?
  my $ol2= length($fa); 
  $fa=~s/^(N+)//; 
  $curlen= length($fa);  
  my $n1= index($fa,'NN'); 
  for( my $iter=0; $n1 >= 0 and $n1 <= $ENDGAP and $iter<5; $iter++) {
    $n1++; $fa= substr($fa,$n1);  
    if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
    $curlen= length($fa); 
    $n1= index($fa,'NN');
  }
  $trimtype.="end5trim" if($curlen < $ol2);

  unless($GAPSOK) {
    my($fac,$cut)= squeezeNNN($fa,0,$GAPSMAX,3,0);
    $fa= $fac; $ncut+=$cut; 
  }
  
  my $nNnew= $fa =~ tr/N/N/;
	my $nNcut = $nNold - $nNnew;
  $ncut = $olen - length($fa);
  return($fa,$trimtype,$ncut,$nNcut,$cdsNold);
}


=item endTrim tests
    
#   if($validcds and ($cdsb < 3*$ENDGAP or $olen - $cdse < 3*$ENDGAP)) {
#     my $cdsfa= substr($fa,$cdsb-1,1+$cdse-$cdsb); my $cdsw=length($cdsfa);
#     my $n1= index($cdsfa,'NN'); 
#     my $ne= rindex($cdsfa,'NN'); 
#     
#     ## this isnt right yet... now other mistake: chopping too many valid cds bases
#     ## unfixit if %N > 50%? or >80%? at ends
#     if($n1 >=0 and $n1 < 2*$ENDGAP) { 
# ## 5d:    
#       my $endc= substr($cdsfa,0,$ENDGAP);
#       $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
#       $fixit |=1 if($zz >= 3);
# ## 5c:      
# #       my $endc= substr($cdsfa,0,2*$ENDGAP);
# #       my $nn= $endc =~ tr/N/N/; 
# #       $fixit |=1 unless($nn >= $ENDGAP);
#      }
#     if($ne >=0 and $cdsw - $ne < 2*$ENDGAP) {
# ## 5d:    
#       my $endc= substr($cdsfa,0,-$ENDGAP);
#       $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
#       $fixit |=2 if($zz >= 3);
# #       my $endc= substr($cdsfa,-2*$ENDGAP);
# #       my $nn= $endc =~ tr/N/N/; 
# #       $fixit |=2 unless($nn >= $ENDGAP);
#     }
#       
#     if($fixit) { 
#       my $cdsfac= $cdsfa;
#       my $ccut=0;
#       ## fixme: endgap only for 1..end5gap, need reverse(cdsfa) for end3gap
#       ## maxgap 4 here, must be bigger than squeezemax of 3
#       if($fixit & 2 == 2) {
#         ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,-2*$ENDGAP);
#         $ncut += $ccut;
#       }
#       if($fixit & 1 == 1) {
#         ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,2*$ENDGAP);
#         $ncut += $ccut;
#       }
#       $cdsfac =~ s/^N+//; $cdsfac =~ s/N+$//;
#       $trimtype.="cdsns";
#       $cdsfac =~ s/N/n/g; # lower to prevent following trim
#       $fa= substr($fa,0,$cdsb-1). $cdsfac . substr($fa,$cdse);
#     }
#   }

## problem 5c: chopping too many valid cds bases; should squeeze this instead.
#BAD neworf uvcut=end5trim,34,0-0,end5,refalignlosscut:-33; litovavel2rk35Loc20449t1: loss named=72%,144/200,148,CDD:201192,UniRef50_B1PT29,UniProt:B1PT29_ARTSF,,Cuticle protein; 
# uncut>litovavel2rk35Loc20449t1 type=cdna; aalen=148,99%,partial; clen=447;  strand=+; offs=2-445;
# CCCTGCTTNNNNNNNNNNNNNNNNNNNNNNNNNN GCCCCTGCTCCTGCTTACAAAGCCCC
# cut5c>litovavel2rk35Loc20449t1 uvcut=end5trim,34,0-0,end5,refalignlosscut:-33; type=cdna; aalen=148,99%,partial; clen=413/447;  strand=+; offs=2-445;
#  GCCCCTGCTCCTGCTTACAAAGCCCCTGAGCCTACCTACTCTGCCCCTTCCCCTAGCTAC

## problem 5b case: mostly NNN at cds end
##BAD neworf uvcut=cdsnsend5trim,26,0-0,end5,refalignlosscut:-27; litovavel1k25Loc19454t5: loss named=75%,308/411,311,CDD:201393,UniRef50_Q9GZS9,UniProt:CHST5_HUMAN,,Carbohydrate sulfotransferase 5; 
# uncut>litovavel1k25Loc19454t5 type=cdna; aalen=311,99%,partial; clen=937;  strand=+; offs=3-935;
# ANNNNNNNNNNNNNNNNNNNNNNNANNNNNN AGCAACGGCCGACGAAGGAGAGAAAGCCA
# cut5b>litovavel1k25Loc19454t5 uvcut=cdsnsend5trim,26,0-0,end5,refalignlosscut:-27; type=cdna; aalen=311,99%,partial; clen=911/937;  strand=+; offs=3-935;
# nAnnn AGCAACGGCCGACGAAGGAGAGAAAGCCAACAGCACCGAGAGCATCATCGCCTCC
#...
#BAD neworf uvcut=cdsnsend5trim,29,0-0,end5,refalignlosscut:-30; litovavel2rk21Loc21431t4: loss named=95%,346/365,342,CDD:173624,UniRef50_E9Q3W1,UniProt:E9Q3W1_MOUSE,,Casein kinase I; 
# uncut>litovavel2rk21Loc21431t4 type=cdna; aalen=342,99%,partial; clen=1029;  strand=+; offs=3-1028;
# CNNNNNNNNNNNNNNNNNNNTNNNNNNNNNNNNNNN CTTCCTCCTCTTCCCTCCTCTCCG
# cut5b>litovavel2rk21Loc21431t4 uvcut=cdsnsend5trim,29,0-0,end5,refalignlosscut:-30; type=cdna; aalen=342,99%,partial; clen=1000/1029;  strand=+; offs=3-1028;
# nnnTnnn CTTCCTCCTCTTCCCTCCTCTCCGCCAAGATGTCTTCGGGAATCATGGGGTGC

=cut
  
# sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
# sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub cutVector {
  my($oid,$fain,$cdsb,$cdse) = @_;
  my $olen=length($fain);
  my($fac,$facgap,$fav,$vectype,$lue,$lub,$ucut,$rue,$rub);
  $fac=$facgap=$fav=$vectype=""; $lue=$lub=$ucut=$rue=$rub=0; 
  return($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype) unless($vecscreen{$oid});
    
  my @vec=split"\n", $vecscreen{$oid}; my $nv=@vec;
  
## moved up....
#  ## fixmed: not wrong order, but overlapped;
#  #MISORDER catfishvel3ik45Loc262t5 uvec 1/2, 528-561 .. 549-567
#  #MISORDER catfishvel3ik45Loc513t4 uvec 1/2, 1-27 .. 15-34
  
  # FIXMEd: special case for vector cuts start-codon, keep ATG (or part missing) in facgap
  my $hasstart= ($cdsb>0 and uc(substr($fain,$cdsb-1,3)) eq 'ATG')?1:0;
  
  for(my $i=0; $i<$nv; $i++) { 
    my($ub,$ue,$vty)=split"\t",$vec[$i]; 
    
      # special case: retain ATG; need also handle uvcut in middle of codon..
      # BUT there are uvectors w/ proteins, drop-all is right answer, this way preserves only 3-base cds and is dropped
      #DROP too short clen=3: >sowhiteflyv8k61loc90865t1 uvcut=Strongc1,1733,1-1736,cdsinend3,shortorfcut:-1176; uvfa=TTCTTATCTCCTTTTGTAGT..TGGCGAGCTGGATGATGAGC; type=cdna; aalen=391,67%,complete; clen=3/1736;  strand=+; offs=132-1307;
      # n=12 cases in whitefly.
      
    my $overstart=($hasstart and $ue >= $cdsb and $ub <= $cdsb+2)?1:0; # special case: retain ATG in facgap
    $overstart=0 if($overstart and $ub <= $cdsb and $ue >= $cdse); # skip prot contained in vector
    if($overstart) {
      my $ue1= $cdsb-1; my $ub2= $cdsb+3; 
      my @add=();
      if($ub < $ue1) { push(@add, "$ub\t$ue1\t$vty"."c1"); }
      if($ub2 < $ue) { push(@add, "$ub2\t$ue\t$vty"."c1"); }
      if(@add) { 
        splice(@vec,$i,1,@add); # replace [$i] + add 1 
        $nv = @vec; 
        ($ub,$ue,$vty)=split"\t",$vec[$i]; # continue on w/ new vec[i]
      } else { next; }
     # now push these into @vec to replace overstart?
    }  
    
    my $uw=1+$ue-$ub; 
    $ucut += $uw;  ## add ucutInCDS += nb for ue - ub in ce - cb
    $vectype= $vty unless($vectype); # take 1st only?
    # $vectype .= $vty unless($vectype =~ /$vty/); # or this?
    if($i == 0) {  $rub=$ub; if($ub>1) { my $s=substr($fain,0,$ub - 1); $fac.=$s; $facgap.=$s; } }
    else { 
      my $wnu= $ub - 1 - $lue; if($wnu>0) { my $s= substr($fain,$lue,$wnu); $fac.=$s; $facgap.=$s; }
      loggit(1, "MISORDER $oid uvec $i/$nv, $lub-$lue .. $ub-$ue\n") if($ub < $lue);
      }
    
##?? add UVGAP only if in cds-span? if($ub < $cdse and $ue > $cdsb) 
## For ~3 cases, this is WRONG, leave out UVgap and get better new aa
## In 1 case, moderate uvec cuts cds5 makes much worse prot; 
## socatfishv1k95loc34193t1 (origaa=446,full and 100% match zfish gene; cutaa=176,utrbad;)
## vecscreen finds ~39 bp align w/ 2+ mismatch in cds5span
## w/o uvcut, full align to alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase [Danio rerio]
## .. use mrna.ann.txt homol-align vs vecsreen-medium to decide not to cut?
## FIXmaybe: some vector spans surrounded by NNN ; cut those?

    ## try both, use neworf to pick best : seems to work  # if($UVGAP)
    if(1) { 
      my $nug=3 + ($uw % 3); my $uvg= substr("nnnnnnnnn",0,$nug); $facgap .= $uvg; 
    }
    
    if($i == $nv-1) { $rue=$ue; if($ue < $olen) { my $s= substr($fain,$ue); $fac.=$s; $facgap.=$s; } }
    $fav.="n" if($i>0 and $ub>$lue+1); $fav .= substr($fain,$ub-1,$uw); 
    $lue=$ue; $lub=$ub; 
  }

  ##? do this here or wait for endtrimNNN
  $fac =~ s/[Nn]+$//; $fac =~ s/^[Nn]+//; #? yes
  $facgap =~ s/[Nn]+$//; $facgap =~ s/^[Nn]+//; #? yes
  $vectype =~ s/ match//g;  

  return($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype);
}
    

sub getbestorf {
  my($id, $hdr,$cdnain,$cdnasize, $oldorfcutlen)= @_;
  
  # my $cmd="$APPcdnabest -nostop -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq"; # -minaa=$MINAA 
  my $MINSIZE4CDS= 60; # what?
  
  if($cdnasize > $MINSIZE4CDS) {
  ## from cdna_bestorf.pl:cdna_bestorf()
  my $fullpart= "best.fwdstrand"; # bothstrand or only fwdstrand ? this is mrna, oriented
  ## ^^ need longorf vs fullorf option
  ## -- use if uvcut-long-partial > uvcut-best but < origaa (dont want > origaa)
  
  my( $bestorf)= getBestProt2( $fullpart, $cdnain);  
  # my( $orfprot, $orient)= orfParts($bestorf, [qw(protein orient)]);
  if(ref($bestorf)) {
    my $orient= undef;
    my $asfasta=1; # ($action =~ /fasta/)
    my($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)
          = proteindoc($bestorf,$cdnasize,$orient,$asfasta);
# warn "best.fwdstrand: $id $aalen aa, cdnasize=$cdnasize\n" if($DEBUG);
          
    if($orflen < 0.98*$oldorfcutlen and $compl =~ /complete/) {
      my($longorf)= getBestProt2( "long.fwdstrand", $cdnain);  # ask for longest partial orf
      my($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1)
           = proteindoc($longorf,$cdnasize,$orient,$asfasta);
# warn "long.fwdstrand: $id $aalen1 aa, cdnasize=$cdnasize\n" if($DEBUG);
           
     if($orflen1 > $orflen) {
      $bestorf= $longorf;
      ($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)=
        ($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1);
     }      
    }
          
    my $cdsfa= $bestorf->{sequence};  $cdsfa  =~ s/(.{60})/$1\n/g;
    (my $cdshdr=$orfprothdr) =~ s/ / type=cds; /;
    ## NOTE: >id is not in these hdr
    return($orfprothdr,$orfprotfa,$cdshdr,$cdsfa, $aalen,$pcds,$compl,$orflen);          
    }
  }
  return("","","","",0,0,0,0);          
}


## move from  evigene/scripts/evgmrna2tsa.pl
sub vecscreen
{
  my($cdnaseq,$vectab,$skiprun)=@_;
  $vectab= makename($cdnaseq,".$VECSUF") unless($vectab);
  unless(-f $vectab) { my($dt,$ft)= getFileset("trimset",$VECSUF); $vectab=$ft if($ft); }
  unless(-f $vectab) { my($dt,$ft)= getFileset('.',$VECSUF); $vectab=$ft if($ft); }
  return($vectab) if( -s $vectab or $skiprun); # or dryrun ..
  
  our($id,$vb,$ve,$ty,$vd,$outh,$inh);  
#   sub putv { our($id,$vb,$ve,$ty,$vd,$outh,$inh); 
#     print $outh join("\t",$id,$vb,$ve,$ty,$vd)."\n" if($id and $ty and not $ty=~/Weak/);  }

  sub putv2 { 
    my($outh,$id,$ty,$vd,@vbe)= @_; 
    return 0 unless($id and $ty and not $ty=~/Weak/i);
    my $no=0; foreach my $vbe (@vbe) { my($b,$e)= @$vbe;
      print $outh join("\t",$id,$b,$e,$ty,$vd)."\n"; $no++;
    } return $no;
  }
  
  ## FIXME: can we use ncbic++ instead? output not same as c-vecscreen... need to check curr ncbi source
  ## ncbic++ doesnt yet support vecscreen .. need own blastn -db UniVec parser to match.. but ncbic- seems obsolete
  
  ##  $ncbi/bin/vecscreen
  my $univecdb="";
  (my $ncbid=$APPvecscreen) =~ s,/vecscreen,/..,; ## want option for db UniVec path
  if( -f "$UniVecDB.nsq") {
    $univecdb= "$UniVecDB"; 
  } elsif( -f "$ncbid/data/$UniVecDB.nsq") {
    $univecdb= "$ncbid/data/$UniVecDB"; 
  } else {
    loggit(1,"ERR: $APPvecscreen missing ../data/$UniVecDB.nsq"); return; 
  }
  
  ## lots of this warn: [vecscreen] WARNING:  [000.000]  Blast: No valid letters to be indexed on context 0
  #old# my $ok= open($inh,"$APPvecscreen -i $cdnaseq -d $univecdb -f3 |");
  
  my ($cmddone,$err)=(0,0);
  my $cmd0="$APPvecscreen -f3 -d $univecdb"; # add -i cdna.split1.fa -o cdna.split1.vec 2> log
  my $vectmp= makename($cdnaseq,".vecscreen.tmp");
  my $veclog= makename($cdnaseq,".vecscreen.log");
  # fixme?? look in trimset/
  # unless(-f $vectmp) { my($ft)= grep /\.vecscreen.tmp/, @$trimd; $vectmp=$ft if($ft); }
  
  if(-s $vectmp) { $cmddone=1; } # skip runcmd if have old data
  elsif( $NCPU > 1 ) {
    my $ccount= facount($cdnaseq); # use this, not fasize
    if($ccount >= 50*$NCPU) {   
      ($err)= vecscreen_ncpu($NCPU,$cmd0,$ccount,$cdnaseq,$vectmp); # cat outparts.
      $cmddone=1; 
    }
  } 
  unless($cmddone) {
    $err= runcmd("$cmd0 -i $cdnaseq -o $vectmp 2> $veclog");
  }
 
  my($ok,$hin)= openRead($vectmp); # my $ok= open($inh,$vectmp); 
  $ok= open($outh,'>',$vectab) if($ok);
  unless($ok) { loggit(1,"ERR: $APPvecscreen -i $cdnaseq -d $univecdb TO $vectab"); return; }
  my($nvid,$so, @vbe); $nvid=0;
  while(<$hin>) {
    chomp; 
    if(/^>/) { 
      $nvid += putv2($outh,$id,$ty,$vd,@vbe) if($id);
      ($id)=m/>Vector (\S+)/; ($vd)=m/Database: (\S+)/; 
      $vb=$ve=$so=$ty=0; @vbe=();
    } elsif(/^No hits/) { 
      $id=0; 
    } elsif(/^\w/) {  
      if(/ match/) { 
        $ty=$_; 
      } elsif(/^Suspect origin/) { 
        $so=1;  # So follows all Match, need to append following b,e to appropriate @vbe offby1 of match range.
      } elsif(/^(\d+)\s+(\d+)/) { 
        my($b,$e)=($1,$2); 
        # if($so and $ve) { $vb=$b if($b<$vb); $ve=$e if($e>$ve);   ## THIS IS WRONG, dont take range, need list of cuts
        if($e < $b) { } # bug? skip
        elsif($so) {  
          for(my $i=0; $i<@vbe; $i++) {
            my($mb,$me)= @{$vbe[$i]};
            if($e == $mb-1) { $vbe[$i]= [$b,$me]; last; }
            elsif($b == $me+1) { $vbe[$i]= [$mb,$e]; last; } # ($b >= $me and $b <= $me+2) is range possible ?
          }
        } elsif($ty) { 
          push @vbe, [$b,$e]; #($vb,$ve)=($b,$e); 
        }  
      } 
    }
  } 
  $nvid += putv2($outh,$id,$ty,$vd,@vbe) if($id);
  close($outh); close($hin);
  
  loggit(0, "vectors found in ntr=",$nvid,$vectab); 
  return($vectab);
}

=item FIXME: bad vecscreen parsing

    ## FIXME: bad parse here for 2+ hits : getting max range, not vec spans
    ## bug in NNN trim: cut to 0 len : n=17 cut to zero for pogonus
    ## -- these are all BIG vectrim = all of cds; need to check that.
    # pogonusvel1pk45Loc1488t3	1	356	Strong match	UniVec_Core   << ** BAD parse of vecscreen.output
    # >Vector pogonusvel1pk45Loc1488t3 Screen: VecScreen Database: UniVec_Core
    # Strong match
    # 3	36
    # 320	356
    # Suspect origin
    # 1	2
    #............
    # pogonusvel1pk45Loc18227t1	1	631	Strong match	UniVec_Core
    ## eg name=Cytochrome c oxidase subunit IV ; PogonEG0004606t2	oid=pogonusvel1pk45Loc1488t3
    ##    oldcds=1-507,vectrim=509,
    ## CDShasTooManyXs:621 #PogonEG0011277t7	oid=pogonusvel1pk45Loc1897t  name=NEL-like	cutcds=339--1,oldcds=339-959,vectrim=961,
    # CDShasTooManyXs:497 #PogonEG0011080t2	oid=pogonusvel1pk45Loc18227t1	len=0; olen=631; nnn=0/631;	name=Glutathione S-transferase, putative	cutcds=134--1,oldcds=134-631,vectrim=631
    # #er2g: PROBLEM: keep unique name but CDShasTooManyXs:354 #PogonEG0004606t2      oid=pogonusvel1pk45Loc1488t3
    # len=0; olen=356; nnn=0/356;     name=Cytochrome c oxidase subunit IV    cutcds=2--1,oldcds=2-355,vectrim=356,
    #.........................
    
=item eg vecscreen 

  >Vector sobeetlepogo1ak39loc7758t1 Screen: VecScreen Database: UniVec (build 6.0)
  Strong match
  1387    1428
  Suspect origin
  1429    1431
  
  >Vector sobeetlepogo1ak31loc100502t1 Screen: VecScreen Database: UniVec (build 6.0)
  Weak match
  9       24    : r1
  563     578   : r2
  Suspect origin
  1       8     << add to 1st range
  579     588   << add to 2nd range

  >Vector sobeetlepogo1ak39loc88670t1 Screen: VecScreen Database: UniVec (build 6.0)
  Strong match
  1       32
  461     490
  Suspect origin
  491     491
  
  updated output vectable:
  sobeetlepogo1ak39loc88670t1     1       32      Strong match    UniVec
  sobeetlepogo1ak39loc88670t1     461     490     Strong match    UniVec << MISSED origin due to b == e

=cut
  


sub vecscreen_ncpu
{
  my($npart,$cmd0,$ccount,$cdnaseq,$vecout)=@_;
  
  # my $ccount= facount($cdnaseq); # use this, not fasize
  my $splcount= int(0.99 + $ccount/$npart);
  my $spldir= makename($cdnaseq,"_vecsplit/",'cdna|fasta|fsa|fa');  # use _tsasubmit/ instead?
  mkdir($spldir); # dryrun?
  
  # push @tmpfiles, $spldir;  ## tmpfiles or erasefiles ?
  
  ## my $err= runcmd("$APPvecscreen -f3 -d $univecdb -i $cdnaseq -o $vectmp 2> $veclog");
  my @splset= fasplitcount( $cdnaseq, $spldir, $npart, $splcount,"fa"); 
  my @vecset;
  my $npartgot= @splset;
  my $icpu= 0;   my $err=0;
  
  #NO: chdir($spldir); ## BAD!! chdir("../"); ?
  ## forkCMD= /home/ux455375/bio/ncbic11/bin/vecscreen -f3 -d /home/ux455375/bio/ncbic11/bin/../data/UniVec_Core\
  #  -i ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.fa \
  #  -o ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.vecout \
  #  2> ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.veclog
  ## errlog: ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split22.veclog : No such file or directory
  
  for(my $ip=0; $ip< $npartgot; $ip++) {
    my $cdna1= $splset[$ip];
    (my $veco1=$cdna1) =~ s/\.\w+$/.vecout/;
    (my $dlog1=$cdna1) =~ s/\.\w+$/.veclog/;
    push @vecset, $veco1;
    my $cmd1= $cmd0 . " -i $cdna1 -o $veco1 2> $dlog1";
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npartgot) { while (wait() != -1) { }; $icpu= 0; }
  }
  while (wait() != -1) { };
  
  # my $cmd= "cat ".join(' ',@vecset)." > $vecout"; runcmd($cmd);
  # preserve dlogs also?
  open(my $outh,'>',$vecout);
  foreach my $vf (@vecset) { 
    if(open(FT,$vf)) { while(<FT>) { print $outh $_; } close(FT);} 
    
    ## if($err||$DEBUG||$dryrun) { push @erasefiles, $vf; } else { unlink $vf if(-f $vf); }
    push @erasefiles, $vf;
    $vf=~s/vecout/veclog/; push @erasefiles, $vf;
    }
  close($outh);
  # rmdir($spldir) unless($err||$DEBUG||$dryrun); ## if($nerase == $npartgot);
  return($err); # $vecout
}


## in cdna_evigenesub
# sub parse_evgheader
# {
#   my($oid,$hdr,$trlen)= @_;
# 
#   ## drop some parts from mrna2tsa for now
#   my $pubid= $oid; #D $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?
# 
#   my $protid= $pubid;
#   $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
#   #D $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
#   #  protein_id      gnl|CacaoGD|Thecc1EG016762p1
# 
#   my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
#       aaqual => "na", trlen => $trlen, cdsoff => "", cdsor => 1, 
#       name => "", namepct => 0, dbxref => "na" ); 
# 
#   if( $genenames{$oid} ) { # $genenames and 
#     $tblinfo{'name'}= $genenames{$oid};
#     $tblinfo{'namepct'}=  $genenamepct{$oid} || 0;
#     $tblinfo{'dbxref'}=  $genedbxref{$oid}||"na";
#     $tblinfo{'cdd'}=  $cddnames{$oid}||"na";
#   }
#         
#   my($cdsb,$cdse,$aafull)=(0,0,0);
#   if($hdr =~ m/\boffs=([\d-]+)/) { my $cdsoff=$1; $tblinfo{'cdsoff'}= $cdsoff; 
#     ($cdsb,$cdse)= split/[-]/,$cdsoff;  } # do in putseq
#   if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo{'aaqual'}= $aq; 
#     ($aafull)= $aq =~ m/(complete|partial\w*)/; }
#   if($hdr =~ m/\bclen=(\d+)/ and not $trlen) { $trlen=$1; $tblinfo{'trlen'}= $trlen; } # skip? 
#   if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo{'cdsor'}= $1; } # expect all '+' from traa2mrna
#   if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct)=([^=\n;]+)/) { my $na=$1; 
#      $tblinfo{'name'}= $na unless($tblinfo{'name'}); }
# 
#   return \%tblinfo;
# }

## in cdna_evigenesub
# sub parse_genenames  # from evgmrna2tsa.pl
# {
#   my($genenames)= @_;
#   my($ngot,$nin)=(0,0);
#   # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
#   %genenames=%genenamepct=%genedbxref=%namedgenes=%cddnames=();
#   
#   open( my $inh, $genenames) or warn("ERR: parse_genenames reading $genenames\n");
#   while(<$inh>) { 
#     next unless(/^\w/ and /\t/);
#     chomp; $nin++;
#     my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name
#     my $xtra; ($name,$xtra)=split";",$name,2; 
#     $name =~ s/\s+$//;
#     
# #     if($pctalign =~/^\d/ and $pctalign < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
# #       if($pctalign >= $MIN_IDLIKE) { unless($name =~ /\blike|^Uncharacterized/) {
# #         $name =~ s/, putative//; 
# #         unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
# #         }
# #       } else { next; } ## should we preserve for ann.txt table ?
# #     }
#     
#     $refid =~ s/RefID://; 
#     $genedbxref{$id} .= "$refid," if($refid);
#     $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
#     $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
#     $cddnames{$id}= $name if($name =~ /CDD:/ and not $cddnames{$id});
#     unless($genenames{$id} and $name =~ /CDD:/) { # or refid =~ /CDD:/
#       $pctalign ||= 0; $refid ||= 0;
#       $genenames{$id}= $name;  $ngot++;
#       $genenamepct{$id}= $pctalign;
#     }
#   } close($inh);
#   
#   return($ngot,$nin);
# }


__END__

