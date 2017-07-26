#!/usr/bin/perl
# evgmrna2tsa.pl; was evgrna2genbanktsa.pl  

=item notes

  EvidentialGene evgmrna2tsa.pl
  process tr2aacds.pl outputs for ncbi tsa  submit
  
  -- main/alternate id table
  -- pubids version of main/alt ids
  -- vecscreen mrna.tr
  -- asmrna2ncbitsa.pl process vecscreen data, annot table?
  -- tbl2asn project.trclass ...

  parts from evigene/scripts/
   evigene2genbanktbl.pl
   asmrna2ncbitsa.pl
   bestgenes_update.pl 
   bestgenes_puban_kfish.pl ??

=item current steps

  1. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
  2. trclass2maintab() : tables of main-alt tr, pubids
  3. vecscreen(mRNA)  : tabulates vec-locs for later putseq that NNN's out vectors
  4. trprocess()/putseq() : link tr,names,pubids; make tbl2asn, annotation file set; trim seq NNN, 
  5. call tbl2asn (as desired)
  
=item fixme vecscreen

  vecscreen can be/is cutting major CDS parts now, with reason..
  should preceed most of this processing, with re-call of .aa, .cds for veccut.mrna

=item new pre-vecscreen steps

  0. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
    -- need genenames/stats for vecscreen counterbalance, ie some vecs are wrong!
    ## socatfishv1k95loc34193t1 (origaa=446,full and 100% align zfish gene; cutaa=176,utrbad;)
    ## vecscreen finds Moderate ~39 bp align w/ 2+ mismatch in cds5span
    ## w/o uvcut, full align to alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase [Danio rerio]
  1a. getmRNA() only   : for vecscreen()  
  1b. vecscreen(mRNA)  : cut out vectors, recall .aa,.cds for those, revise evgtrset 
     .. reprocess evgtrset, new mRNA fileset, ..
     -- add trim NNN ends here? so new mRNA set is clean of that
     
  2. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
  3. trclass2maintab() : tables of main-alt tr, pubids
  4. trprocess()/putseq() : link tr,names,pubids; make tbl2asn, annotation file set; trim seq NNN, 
      -- ? reannotate .mrna,.cds,.aa fileset with pubids, select .ann.txt (gene names, dbxref,..)?
      -- for public use outside ncbi tsa submit
  5. call tbl2asn (as desired)
  
=item old script

  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl

=item tbl2asn test

  need these input templates; generate some from other configs?
    evgr_tsamethods.cmt evgr_tsadesc.cmt evgr_tsasubmit.sbt
  
  pt=tsasub1 
  cp -p ../tsasubmit/evgr_*.{cmt,sbt} $pt/   
  # .. edit cmt per needs .. 
  # Should add PRJ of data source: PRJNA73443 .. where? in evgr_tsamethods.cmt ?
  
  org='[organism=Litopenaeus vannamei] [bioproject=PRJNA12345]'
  sra='[SRA=SRR346404]' 
  
  $ncbin/tbl2asn -p $pt/ -Z $pt/$pt.discrep.log \
     -w $pt/evgr_tsamethods.cmt -Y  $pt/evgr_tsadesc.cmt -t  $pt/evgr_tsasubmit.sbt \
     -a r10k  -l paired-ends -Vtb -Mt -XE \
     -j "[moltype=mRNA] [tech=TSA] $org $sra"

=cut

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
#maybe# use cdna_proteins;
#add # use cdna_evigenesub;

use constant VERSION => '2013.05.20'; # 07; # 05; # '04.20'; # .16'; # 03.20
   # 13.05.07: fix big vecscreen parsing bug .. redo all data from that
   
## evigene path = self path
my $EVIGENES="$FindBin::Bin"; #??
### new cdna_evigenesub globals:
# our $EVIGENES="$FindBin::Bin/..";  
# our $EGAPP='egtrimvec';  
# our $dryrun=0; ## $DRYRUN ?
# our $DEBUG= $ENV{debug}|| 0;

## Evigene tr2aacds subdirs: see tidyup
## add for trimvec, mrna2tsa:  trimset?  publicset? tsasubmit/submitset ?
## change trimvec,mrna2tsa output subdir: publicset? adding pubids, main2alt, ...
## separate tsasubmit subdir ..
our @evgdirs = qw(okayset dropset inputset tmpfiles erasefiles publicset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset); # tidyup file sets
our $outdir='publicset'; # pubset ?

my $debug= 0;
my $GENEINFO_VERS=2; # constant
my $MAXGAP=15; # NCBI 
my $ENDGAP=20; ## was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
  ## FIXME, NCBI tbl2asn changed this to 20 !
  ## [SEQ_INST.HighNContentStretch] Sequence has a stretch of at least 10 Ns within the first 20 bases BIOSEQ: lcl|MusacuEGm017348t1: delta, rna len= 509

my $MINSIZE=200; # NCBI
my $GAPSOK=1; # default on? new policy 2012Dec for TSA tbl2asn: -a r10u -l paired-ends
my $MINGENEIDENT=85; # for asm == gene identity
my $IDPREFIX= $ENV{idprefix} || "EVGmRNA"; ## "evgr"; #  opt
my $GDB_PREFIX='gnl|Evigene|';  #see below; 'gnl|Evigene|'; # use IDPREFIX ? or not, since ID has this also
my $pubidnum_start=0;
my $NCPU= 1; 
my $UniVecDB= $ENV{UniVec} || "UniVec_Core"; # Not UniVec

my $DOtbl2asn=0; 
my($organism,$sraids)=("Noname","SRR000000");  ## SRR 346404
my $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
    ## ^^ need config for these; generate some of this?
my $DATE=`date '+%Y%m%d'`; chomp($DATE); # default is today; use perl func?

## namegenes messed up, didn't screen out loqualnames:
##  1%,309/22971,716        Dumpy
##  0%,126/34350,248        Titin
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname

my ($vecscreenf,$trclass,$geneinfo,$genenames,$cdnaseq,$output,$logfile,$tblfile,$runsteps,$dryrun);

my %DEFAULT_SETTINGS= (
  IDPREFIX=>$IDPREFIX, TSADESC=>$TSADESC, DATE=>$DATE, 
  MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
  sraids => $sraids, organism => $organism,
  trclass => '', mrna => '', genenames=>'', vecscreen => '',
  );
  

# FIXME below: save Options to trname.info.stash : TSADESC, organism, ... 
my @saveopt= grep /^\-/, @ARGV;
my $optok= GetOptions(
  # "config=s", \$config,
  # "cadd=s", \@configadd,
  "mrna|cdna=s", \$cdnaseq,
  "class|trclass=s", \$trclass,
  "vectors|vecscreen=s", \$vecscreenf,  
  "names|genenames=s", \$genenames, ## ? allow for 2 files: myspecies.namerefids + allrefprot.names
  "organism|species=s", \$organism,   
  "sraids=s", \$sraids,   
  "TSADESC=s", \$TSADESC,   
  # "geneinfo=s", \$geneinfo,  # old
  "output:s",  \$output,
  "tblfile:s", \$tblfile,  
  "logfile:s", \$logfile,
  "idprefix=s", \$IDPREFIX,  
  "gdbprefix=s", \$GDB_PREFIX, #?? IDPREFIX default 
  "DATE=s", \$DATE,  
  # "proteins=s", \$proteins, # maybe
  # "version|MSRC=s", \$MySRC,
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  "NCPU=i", \$NCPU,## "MAXMEM=i", \$MAXMEM,  
  "runsteps=s", \$runsteps,   
  "runtbl2asn!", \$DOtbl2asn, 
  # "MINGENEIDENT=i", \$MINGENEIDENT,  ## not used?
  "GAPSOK!", \$GAPSOK, # default on always?
  "dryrun|n!", \$dryrun, 
  "debug!", \$debug, 
  );


# OPTION: add -outdir opt ; keep list of created files, move there, or work there..
# OPTIONed: -class evg_tr2aacds.trclass only requirement, gives paths to mrna, names.tab ?
# OPTIONed: here or caller? make mrna/cdnaseq  from evg traa2cds.pl okayset/name.{okay,okalt}.{tr,aa}  
# $evigene/scripts/prot/traa2cds.pl -trout -cdna $pt.tr.gz -aa $pt.aa.gz -out -log
##  -trout : cdna output, not cds, as name.mrna.tr
#?No# $cdnaseq= shift @ARGV unless($cdnaseq);
## caller script bug: #er2g: ERR: unused arguments: vannamei'
## from  opts="$opts -species='$species'"

die "EvidentialGene mrna2tsa VERSION ",VERSION,"
  make mRNA fasta and annotation table for TSA submit via tbl2asn
  using mRNA classes from Evigene tr2aacds, gene product names, vecscreen
Usage: evgmrna2tsa.pl -trclass mrna.trclass -mrna mrna.fasta ...
opts: -idprefix Thecc1EG -names mrna.names -log -dryrun -debug
    -runtbl2asn -organism=$organism -SRAids=$sraids -TSADESC='$TSADESC'
    -NCPU=$NCPU -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP\n"
  unless($optok and ($cdnaseq or $trclass));  
#  -out=out.fsa  -conf=evigene.conf  -vecscreen=infile 

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, };
my $logh= undef;
sub loggit{ my $dowarn=shift; 
  my $s= join(' ',@_); chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($logh){ print $logh "#er2g: $s\n"; } elsif($dowarn>0||$debug){ warn "#er2g: $s\n"; }
  if($dowarn == LOG_DIE) { die "#er2g: $s\n" ; }
}

if(not $logfile and defined $logfile) { # use output name
  $logfile= $cdnaseq || $trclass;
  $logfile= makename($logfile,".mrna2tsa.log");  
}
if($logfile) { open(LOG, ">>$logfile") or die $logfile; $logh= *LOG; } # open($logh,'>>',$logfile) or die;

my $GAPSMAX = ('N') x $MAXGAP;
# evigene_config($config, \@configadd); # always even if $config null
unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # fixme

 
# #*** defer IDPREFIX use till read sra_result.cvs ..
# #  get species/organism from that, make IDPREFIX from Spp.org. abbev. 
# my $ChangeDefaultIdPrefix= ($IDPREFIX eq "EVGmRNA")?1:0; 
# my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : 6;
# my $pubid_format = $IDPREFIX.'%0'.$nd.'d'; # $public_options{'publicid'} || "evgr000000";
# my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
# my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";
# #? No?# unless($GDB_PREFIX) { $GDB_PREFIX= "gnl|$IDPREFIX|"; } #?

loggit(1, "EvidentialGene mrna2tsa.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

my $APPvecscreen=  findapp("vecscreen"); #if($DOvecscreen); ?? ncbi/bin/... also need ncbi/data/UniVec
 ## find this also $UniVecDB="UniVec_Core"; # not UniVec

my $APPtbl2asn  =  findapp("tbl2asn") if($DOtbl2asn); # add this, needs configs + data doc files (.cmt, .sbt,..)
## .. tbl2asn can be real slow .. use fasplit(ncpu); forkcmd(tbl2asn,ncpu,..)
# my $APPmakeblastdb= findapp("makeblastdb");
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
# my $APPcdnabest= "$EVIGENES/cdna_bestorf.pl"; # allow ENV/path substitutions?
# my $APPaaqual=    "$EVIGENES/prot/aaqual.sh";

## FAIL at this point if any apps missing?
#-------------------------------------


my( %vecscreen, %gene, %genenames, %genedbxref, %genenamepct, %namedgenes, %cddnames, 
    %pubids, %geneinfo, $trpath,$trname,$sradatah, %settings);

do_settings("restore",$trclass||$cdnaseq,); # ("log|restore|save");

## add OPTION $steps : what to do/not here ??
$runsteps||="";
# use e.g. to make only trclass2maintab
my $skipmaintabrun=($runsteps and $runsteps !~ /main|ids|pub/) ? 1 : 0;
my $skipvecrun=($runsteps and $runsteps !~ /vec/) ? 1 : 0;
my $skiptrrun= ($runsteps and $runsteps !~ /tblout|process/) ? 1 : 0; # annot and ncbi input/output files
$DOtbl2asn= ($runsteps and $runsteps !~ /asn|ncbi|tbl2/) ? 0 : $DOtbl2asn;
$DOtbl2asn=0 if($skiptrrun);

($cdnaseq,$trpath,$trname,$sradatah)= get_evgtrset($trclass,$cdnaseq,);
  loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdsseqnr)
  loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

## FIXMEd below: trclass2maintab $pubids may exist ; use that for IDPREFIX if so ..
my($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX(); # default abbrev of $organism now

## FIXME: write some of the collected info to $trname.info.stash or such, so don't have to
##    refind it: ftp SRR ids, -opts like IDPREFIX, .. @saveopt

my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass);
loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 
if($pubids) { # require?
  open(F,$pubids) or loggit(1,"ERR: reading $pubids");
  my $first=0;
  while(<F>) { next unless(/^\w/); chomp; 
    my($id,$oid,$gid,$alti)=split"\t";     
    $pubids{$oid}= $id; $pubids{$id}="$oid\t$gid\t$alti";
    if(++$first == 1 and $nalltr == 0) { # reset to existing IDPREFIX
      ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($id); 
    }
  } close(F);
}

#/////////////////////////
## moved to asmrna_trimvec.pl
($vecscreenf)= vecscreen($cdnaseq,$vecscreenf||"",$skipvecrun);
if($vecscreenf) { 
  my %vids=();
  open(F,$vecscreenf) or loggit(1,"ERR: reading $vecscreenf");
  while(<F>) { chomp; my($vd,$vb,$ve,$vt)=split"\t"; $vecscreen{$vd}.="$vb\t$ve\t$vt\n"; $vids{$vd}++; } close(F);
  my $nvid=scalar(keys %vids); 
  loggit(0, "vectors found in ntr=",$nvid,$vecscreenf); 
}
#/////////////////////////

do_settings("log|save",$trclass||$cdnaseq,); # or after last call?  ("log|restore|save");

# REVISE trprocess to NOT do vecscreen or gaptrim .. now in asmrna_trimvec
my($outfa, $tblout, $annot,$notr,$nocds) = trprocess($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
loggit(0,"DONE output ntr=$notr, ncds=$nocds in files $maintab, $pubids, $outfa, $tblout, $annot"); 

if($DOtbl2asn) {
  my($npartgot,$spldir,$sqnoutlist)= tsa_tbl2asn($outfa,$tblout,$organism,$sraids);
  loggit(0,"DONE tsa_tbl2asn nparts=$npartgot, submitset=$sqnoutlist"); #?? $spldir  
}


#---------------------------------------------------  

sub do_settings {
  my($action,$pathname)= @_;
  ## write these to work dir; reread next go
  ## action == 'log|save|restore' ; restore called AFTER read new options, shouldnt replace
  my $PRES='m2t';
  my $trpname= makename($pathname,".mrna2tsa.info"); 
  $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; # ??
  
  my %mysettings= ( # current; global is %settings, and DEFAULT_SETTINGS before GetOptions
  IDPREFIX=>$IDPREFIX, TSADESC=>$TSADESC, DATE=>$DATE, 
  MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
  sraids => $sraids, organism => $organism,
  trclass => $trclass, mrna => $cdnaseq, genenames=>$genenames,
  vecscreen => $vecscreenf,
  );
  
  if($action =~ /restore|read/ and -f $trpname) {
    open(my $inh, $trpname); # or loggit(warn ..);
    while(<$inh>) { chomp; if(s/^$PRES.//) { 
      my($k,$v)=split /\s*=\s*/,$_,2; $v=~s/;*$//;
      my $ov= $mysettings{$k}; 
      unless($ov and $ov ne $DEFAULT_SETTINGS{$k}) { 
        ## fixme: need to reset global defaults
        $organism=  $v if($k eq 'organism');
        do { $v=~s/ +;/;/g; $v=~s/ +/;/g; $sraids=    $v; } if($k eq 'sraids');
        $IDPREFIX=  $v if($k eq 'IDPREFIX');
        $TSADESC=   $v if($k eq 'TSADESC');
        $DATE=      $v if($k eq 'DATE');
        $MINSIZE=   $v if($k eq 'MINSIZE'); # require int
        $MAXGAP=    $v if($k eq 'MAXGAP'); # require int
        $genenames=  $v if($k eq 'genenames');
        # $trclass=  $v if($k eq 'trclass');
        # $cdnaseq=  $v if($k eq 'mrna'); # fixme: key != varname
        # $vecscreen=  $v if($k eq 'vecscreen');
        $mysettings{$k}=$v; # after possible changes
        } 
      }
    } close($inh);
  }
  
  %settings = %mysettings; # make global now
  my $settings= join "\n", map{ "$PRES.$_=".$settings{$_} } sort keys %settings;
  if($action =~ /log/) { loggit(0, "mrna2tsa.info:\n$settings");  }
  if($action =~ /save/) { open(my $outh, '>', $trpname); print $outh $settings,"\n"; close($outh); }
}

=item get_evgtrset

  my($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,);

  add: parse sra_result.csv if exists, for species, sraids .. other metadata
  see new in asmrna_trimvec.pl
    
=cut

sub get_evgtrset {
  my($trclass,$cdnaseq)= @_;
  my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);
  
  if($cdnaseq and !$trclass) { $trclass= makename($cdnaseq,".trclass"); }
  if($trclass) {
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # BADDDDD
    else { $trname=$trpname; }
    $trpath ||= '.';  
    
    my $sracvs="$trpname.sra_result.csv";    
    $sracvs="$trpath/sra_result.csv" unless(-f $sracvs);  
    ($nsra,$sradatah)= parse_sra_result_cvs($sracvs) if(-f $sracvs);
    loggit(0,"sra_result from",$sracvs,"nresult=",$nsra);
    
     ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
    my($nupinfo,$tsamethf,$tsadescf,$tsasubf)= tsa_infotemplates($trpath, $trname, $sradatah);
    loggit(0,"info updated $nupinfo TSADESC=",$TSADESC);
   
    unless($genenames) { #? put in get_evgtrset
      my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
      loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
    }

    my $okpath="$trpath/okayset"; # fixme change to publicset for output, okayset for input
    ($cdnaseq)= getmRNA($okpath,$trname) if(!$cdnaseq and -d $okpath);

  }
  
  return($cdnaseq,$trpath,$trname,$sradatah);
}

# see new in asmrna_trimvec.pl
sub getmRNA
{
  my($okpath,$trname)= @_;
  my($cdnaseq)=(""); # == mrna (oriented), not cdna.parts.tr
  
  #? FIXME? suffix .tr may not be used: .cdna? .fa? .fasta? ...
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna  
  
  opendir(D,$okpath); 
  my @okd= map{ "$okpath/$_" } readdir(D); 
  closedir(D);

  ## Ugh! bad call: #er2g: get_evgtrset= ./okayset/pogonus1all3.mrna0.ann.txt . pogonus1all3
  ## need grep /\.mrna$|\.mrna.gz$|.mrna.fa$/ ??? 
  my ($trf)= grep /\.mrna$|\.mrna.gz$/, @okd;
  if($trf and -s $trf) { $cdnaseq= $trf; }
  else { 
    my $okall=0;
    my $cdnatmp="$okpath/$trname.mrna";
    loggit(0,"Make mRNA $cdnatmp from okayset transcripts");
    my($oktr) = grep /.okay\.($TRSUFFIX)/, @okd;  
    my($alttr)= grep /.okalt\.($TRSUFFIX)/, @okd; 
    my $okaa= makename($oktr,'.aa');   ##mess (my $okaa=$oktr) =~ s/\.gz//; $okaa=~s/\.\w+$/.aa/;
    my $altaa=makename($alttr,'.aa');  ##(my $altaa=$alttr)=~ s/\.gz//; $altaa=~ s/\.\w+$/.aa/;
    #.. problem here .aa.gz not .aa
    $okaa= "$okaa.gz" if( -f "$okaa.gz");
    $altaa= "$altaa.gz" if( -f "$altaa.gz");
    if($oktr and -f $oktr and -f $okaa) { 
      my $err= runcmd("$APPtraa2cds -trout -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    if($alttr and -f $alttr and -f $altaa) { 
      my $err= runcmd("$APPtraa2cds -trout -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    $cdnaseq= $cdnatmp if(-s $cdnatmp);
    loggit(1,"ERR: No-make mRNA $cdnatmp from $oktr,$okaa + $alttr,$altaa") unless($okall>1);
  }
  return($cdnaseq);      
}

# my($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX();
sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  #?? defer IDPREFIX use till read sra_result.cvs ..
  #  get species/organism from that, make IDPREFIX from Spp.org. abbev. 
  
  # my $ChangeDefaultIdPrefix= ($IDPREFIX eq "EVGmRNA")?1:0; 
  my $ChangeDefaultIdPrefix= ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
    $ChangeDefaultIdPrefix=0;
  }
  
  if($ChangeDefaultIdPrefix and not($organism eq "Noname" or $organism !~ m/\w\w/)) {
    my $prefix="";
    my($gen,$spp)=split" ",$organism;
    if($spp and length($gen)>1) { $prefix= ucfirst(substr($gen,0,3)) . lc(substr($spp,0,3)) . "EGm"; }
    else { $prefix= ucfirst(substr($organism,0,6)) . "EGm"; }
    ##else { $prefix=$IDPREFIX; } # make default?
    $IDPREFIX= $prefix;
  }
  my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  my $pubid_format = $IDPREFIX.'%0'.$nd.'d'; # $public_options{'publicid'} || "evgr000000";
  my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
  my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";
  #? No?# unless($GDB_PREFIX) { $GDB_PREFIX= "gnl|$IDPREFIX|"; } #?
  return($pubid_format,$altid_format,$GBPROID);
}



=item tsa_infotemplates

  - Make these if needed; update for trasm name as needed

evgr_tsamethods.cmt:
  StructuredCommentPrefix ##Assembly-Data-START##
  Assembly Method EvidentialGene v2013.03.02; << from VERSION
     Velvet/Oases v1.2.03/o0.2.06 (2012.02); SOAPdenovo-Trans v2011.12.22; Trinity v2012.03.17  << USER inputs
  Assembly Name   evigeneR13     << from local myspecies.trclass name
  Sequencing Technology   Illumina   << from sra_result

evgr_tsadesc.cmt : insert parts of sra_result
  RNA-Seq data of [[this species]] are
  assembled de-novo with [[various RNA assemblers]], using multiple options
  for kmer fragmenting, insert sizes, read coverage, quality and
  abundance filtering. 
  EvidentialGene tr2aacds pipeline software is used
  to process the [[several million]] resulting assemblies by coding
  sequences, translate to proteins, score gene evidence including
  CDS/UTR quality, homology, and classify/reduce into a biologically
  informative transcriptome of primary and alternate transcripts.

evgr_tsasubmit.sbt
Submit-block ::= {
  contact { contact { name name { ... } }
  ..
  Seqdesc ::= user { type str "DBLink",
  data { { label str "BioProject", num 1, data strs { "PRJNA99999" } } } << need PRJNA here..
  }
  

=cut

sub tsa_infotemplates  # make if not found as files; edit w/ sra_result info, other
{
  my($trpath, $trname, $sradata)= @_; ## add sraresult hash info
  my($methtxt,$desctxt,$subtxt)=("") x 3;
  my($tsamethf,$tsadescf,$tsasubf)=("") x 3;
  my $update=0; my $nupdate=0;
  
  # FIXME: must have \tabs in StructuredComment
  $tsamethf= "$trpath/$trname.tsamethods.cmt";  $update=0;
  unless( -s $tsamethf) {  $tsamethf= "$trpath/evgr_tsamethods.cmt"; $update=1; }
  if( -s $tsamethf) {
    open(F, $tsamethf); $methtxt= join "", <F>; close(F);
  } else {
    $update=1; $methtxt=""; # TEMPlATE
    my $Illumina= $sradata->{"Instrument"} || "Illumina"; # from sradata
    my $asmname="$trname.evigene"; # from $trname
    my $egvers= VERSION;
    $methtxt=<<"EOT";
StructuredCommentPrefix\t##Assembly-Data-START##
Assembly Method\tEvidentialGene $egvers
Assembly Name\t$asmname
Sequencing Technology\t$Illumina
EOT
  }
  if($update) {
    $tsamethf= "$trpath/$trname.tsamethods.cmt"; $nupdate++;
    open(O,'>',$tsamethf); print O $methtxt,"\n"; close(O);
  }
  
  $tsadescf= "$trpath/$trname.tsadesc.cmt"; $update=0;
  unless( -s $tsadescf) {  $tsadescf= "$trpath/evgr_tsadesc.cmt"; $update=1; }
  if( -s $tsadescf) {
    open(F, $tsadescf); $desctxt= join "", <F>; close(F);
  } else {
    $update=1; $desctxt=""; # TEMPlATE
    my $asmsoft="various RNA assemblers"; #?
    my $datasize= $sradata->{"Total Size, Mb"} || ""; ## "Total Size, Mb"
    my $nspots= $sradata->{"Total Spots"} || "";
    my @ds=split";",$datasize;  my @ns=split";",$nspots; 
    if(@ds>1) { my $ds=0; my $ns=0; 
      for(my $i=0; $i<@ds; $i++) { my $n;
        ($n)= $ds[$i] =~ m/(\d+)/; $ds+=$n; 
        ($n)= $ns[$i] =~ m/(\d+)/; $ns+=$n; 
        } 
      $datasize=$ds; $nspots= $ns; }
    $datasize.=" Mb in $nspots read-pairs";
    my $asmcount="many"; # from where? tr2aacds.log has counts
    $desctxt=<<"EOT";
RNA-Seq data, $datasize, of $organism are assembled with
$asmsoft, using multiple options. EvidentialGene tr2aacds pipeline
software is used to process the $asmcount resulting assemblies by
coding sequences, translate to proteins, score gene evidence and
classify/reduce to a biologically informative transcriptome of primary
and alternate transcripts.
EOT
  }
  if($update) {
    $tsadescf= "$trpath/$trname.tsadesc.cmt"; $nupdate++;
    open(O,'>',$tsadescf); print O $desctxt,"\n"; close(O);
  }

  $tsasubf= "$trpath/$trname.tsasubmit.sbt"; $update=0;
  unless( -s $tsasubf) {  $tsasubf= "$trpath/evgr_tsasubmit.sbt"; $update=1; }
  if( -s $tsasubf) {
    open(F, $tsasubf); $subtxt= join "", <F>; close(F);
  } else {
    $update=1; $subtxt=""; # TEMPlATE cant really guess this.
    ## need "BioProject",  "PRJNA99999" from somewhere ..
  }
  if($update) {
    $tsasubf= "$trpath/$trname.tsasubmit.sbt"; $nupdate++;
    open(O,'>',$tsasubf); print O $subtxt,"\n"; close(O);
  }

  ## edit TSADESC if($nupdate or paths not in $TSADESC)
  ## my $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
  unless($TSADESC =~ m/$tsamethf/ and $TSADESC =~ m/$tsadescf/ and $TSADESC =~ m/$tsasubf/ ) {
    $TSADESC="-w $tsamethf -Y $tsadescf -t $tsasubf"; 
    #see.above# loggit(0,"info updated $nupdate TSADESC=",$TSADESC);
  }
  
  return($nupdate,$tsamethf,$tsadescf,$tsasubf);
}



=item sra_ftp2srr from FTPpath

  curl  'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  dr-xr-xr-x 1073741824 ftp      anonymous        0 Aug 12  2012 SRR424344
  
  >> best:
  curl -s -l 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX160070/'
  SRR521835
  SRR521944
  
  wget -A 'sra' -r -l 2 -nv -nd -nH -np --spider \
   'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  14:25:43 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/SRR424344.sra [0] -> "SRR424344.sra" [1]

=cut

sub sra_ftp2srr {
  my(@ftps)= @_;  return () unless(@ftps);
  my $APPcurl= findapp('curl'); return () if($APPcurl =~ /MISSING/);
  my @srrs=(); foreach my $ftppath (@ftps) {
    my $srrs= `$APPcurl -s -l $ftppath/` if($ftppath =~ m,^ftp://,); # or http: ??
    push @srrs, grep /SRR/, map{ m/(SRR\w+)/; $1; } split " ",$srrs;
    }
  return @srrs;
}


=item parse_sra_result_cvs

  add: parse sra_result.csv if exists, for species, sraids .. other metadata
  evgr2tsa/litova1all3f/
    whiteshrimp_sra_result.csv # rename litova1all3.sra_result.csv

  "Experiment Accession","Experiment Title","Organism Name","Instrument",
    "Submitter","Study Accession","Study Title","Sample Accession","Sample Title",
    "Total Size, Mb","Total RUNs","Total Spots","Total Bases","FTP Path to Experiment",
    "Library Name","Library Strategy","Library Source","Library Selection"
    
  "SRX098246","Litopenaeus vannamei  transcriptome","Litopenaeus vannamei","Illumina HiSeq 2000","BGI",
    "SRP008317","BGI Litopenaeus vannamei transcriptome sequencing","SRS265043","Pacific white shrimp",
    "1515.37","1","13697473","2465545140","ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX098/SRX098246",
    "Ex-zai-2_l1","RNA-Seq","TRANSCRIPTOMIC","cDNA"
  
=cut

sub parse_sra_result_cvs
{
  my($sracvs)= @_;
  my($ngot,$nin)=(0,0);
  my (%sradata, @hd);
  
  open( my $inh, $sracvs) or loggit(1,"ERR: parse_sra_result_cvs reading $sracvs");
  while(<$inh>) { 
    next unless(/,/); chomp; 
    my @cols= map{ s/^"//; s/"$//; $_; } split /\",\s*\"/, $_;
    if($cols[0] =~ /Experiment Accession/) { @hd= @cols; }
    else { 
      $ngot++;
      for(my $i=0; $i<@cols; $i++) { 
        my $hd=$hd[$i]||$i; my $v=$cols[$i]; 
        $sradata{$hd}.="$v;" unless($sradata{$hd} and $sradata{$hd} eq "$v;"); 
        }
    }
  } close($inh);
  
  foreach my $k (keys %sradata) { $sradata{$k}=~s/;$//; }
  
  #reset defaults:  
  # my($deforg,$defsra)=("Noname","SRR000000"); 
  my($deforg,$defsra)= ($DEFAULT_SETTINGS{'organism'},$DEFAULT_SETTINGS{'sraids'});
  my($sorg,$sids)    = ($sradata{"Organism Name"},$sradata{"Experiment Accession"});
  $organism= $sorg if($sorg and ($organism eq $deforg or $organism !~ m/\w\w/));
  $sraids=   $sids if($sids and ($sraids eq $defsra or $sraids !~ m/SRR/));

  my @SRAK=("Experiment Accession","Organism Name","Instrument",
            "Submitter", "Total Size, Mb","Total Spots");
  if($ngot>0) { my @v= map{ $_ .'='. $sradata{$_}.';'; } @SRAK; loggit(0,"sra_result:",@v); }
  
  ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
  ## rewrite template .cmt, .sbt unless updated already.

  ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
  my $NEEDsrrids=0; 
  $NEEDsrrids= ($sraids =~ /SRR/ and $sraids ne $defsra)?0:1; ## Need to save to info/config file
  if($NEEDsrrids) {
    my @ftps=  grep /ftp/, split/;/, $sradata{"FTP Path to Experiment"};  
    my @srr= sra_ftp2srr(@ftps);
    if(@srr>0) {
      $sraids= join(";",@srr);
      $sradata{'SRAids'}= $sraids;
      loggit(1,"sra_id=",$sraids);
      }
  }
   
  return($ngot, \%sradata);
}

## in cdna_evigenesub.pm
sub parse_genenames
{
  my($genenames)= @_;
  my($ngot,$nin)=(0,0);
  # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
  
  ## FIXME2: ** use uniq names vs ERR: too short to keep valid tr, e.g. 
  #er2g: ERR: too short:183 #LitvaEG0018688t4     oid=litovavel2k35Loc15824t1     len=183 name=CDD: clpS, ATP-dependent Clp protease ada..-like    
  # grep  'CDD: clpS,' *.names = 1 only = litovavel2k35Loc15824t1
  
  # FIXME: need better reader; 2+ rows/id; pick best .. format may change..
  # names.tab ==  id, name, pctalign, refid, repid  : now
  #  trid1  C-ets-2 protein, putative       89%,103/116,197 RefID:UniRef50_E0VFI2   RepID:E0VFI2_PEDHC
  #  trid2  DBH-like monooxygenase protein 1        73%,445/613,516 RefID:UniRef50_Q6UVY6   RepID:MOXD1_HUMAN
  
  open( my $inh, $genenames) or loggit(1,"ERR: parse_genenames reading $genenames");
  while(<$inh>) { 
    next unless(/^\w/ and /\t/);
    chomp; $nin++;
    my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name
    my $xtra; ($name,$xtra)=split";",$name,2; 
    $name =~ s/\s+$//;
    
    # FIXME: 2 names/id maybe: CDD: xxx and gene xxx; keep both in ann.txt ? and pctalign?
    ## pctalign == 100%,450/450,446 : pct,naln/nref,ntrg
    ## old geneinfo:
    #  ($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
    ## usage below in putseq
    # my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
    #  aaqual => "na", trlen => 0, cdsoff => "", cdsor => 1, 
    #  name => $genenames{$oid}||"", dbxref => $genedbxref{$oid}|| "na"); 

    if($pctalign =~/^\d/ and $pctalign < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
      ## bad: Uncharacterized protein-like ; Nuclease HARBI1, putative-like; Protein-like
      if($pctalign >= $MIN_IDLIKE) { unless($name =~ /\blike|^Uncharacterized/) {
        $name =~ s/, putative//; 
        unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
        }
      } else { next; } ## should we preserve for ann.txt table ?
    }
    
    ## fixme: CDD:206692,cd04107,RefID:UniRef50_Q9NX57,UniProt:RAB20_HUMAN, 
    ## drop  RefID:; drop? cd04107
    $refid =~ s/RefID://; 
    $genedbxref{$id} .= "$refid," if($refid);
    $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
    $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
    
    ## FIXME: keep CDD names for .ann.txt, maybe .tbl submit as 2nd note
    $cddnames{$id}= $name if($name =~ /CDD:/ and not $cddnames{$id});
    
    unless($genenames{$id} and $name =~ /CDD:/) { # or refid =~ /CDD:/
      $pctalign ||= 0; $refid ||= 0;
      $genenames{$id}= $name;  $ngot++;
      $genenamepct{$id}= $pctalign;
      ## $genedbxref{$id}= ($repid) ? $repid : $refid;  # do list here for all dbxref
       # repid Yes or No? this is by default RefID=UniRef50_xxxx and RepID=UniProt:xxxx_HUMAN
    }
  } close($inh);
  
  return($ngot,$nin);
}


## in cdna_evigenesub.pm
sub parse_evgheader
{
  my($oid,$hdr,$trlen)= @_;

  my $pubid= $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?

  my $protid= $pubid;
  $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
  $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
  #  protein_id      gnl|CacaoGD|Thecc1EG016762p1

  my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      aaqual => "na", trlen => $trlen, cdsoff => "", cdsor => 1, 
      name => "", namepct => 0, dbxref => "na" ); 

  if( $genenames and $genenames{$oid} ) {
    $tblinfo{'name'}= $genenames{$oid};
    $tblinfo{'namepct'}=  $genenamepct{$oid} || 0;
    $tblinfo{'dbxref'}=  $genedbxref{$oid}||"na";
    $tblinfo{'cdd'}=  $cddnames{$oid}||"na";
  }
        
  my($cdsb,$cdse,$aafull)=(0,0,0);
  if($hdr =~ m/\boffs=([\d-]+)/) { my $cdsoff=$1; $tblinfo{'cdsoff'}= $cdsoff; 
    ($cdsb,$cdse)= split/[-]/,$cdsoff;  } # do in putseq
  if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo{'aaqual'}= $aq; 
    ($aafull)= $aq =~ m/(complete|partial\w*)/; }
  if($hdr =~ m/\bclen=(\d+)/ and not $trlen) { $trlen=$1; $tblinfo{'trlen'}= $trlen; } # skip? 
  if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo{'cdsor'}= $1; } # expect all '+' from traa2mrna
  if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct)=([^=\n;]+)/) { my $na=$1; 
     $tblinfo{'name'}= $na unless($tblinfo{'name'}); }

  return \%tblinfo;
}

  # my($outfa,$tblout,$annot,$notr,$nocds) = trprocess($cdnaseq,$trclass);    
sub trprocess
{
  my($cdnaseq,$trclass,$skiprun)=@_;
  my($notr,$nocds)=(0,0);
  my $outfa = ($output)  ? $output  : makename($cdnaseq,".fsa"); # was .fna ; tbl2asn wants fsa
  my $tblout= ($tblfile) ? $tblfile : makename($cdnaseq,".tbl");
  my $annot =  makename($cdnaseq,".ann.txt"); ## FIXME: change suffix: was .annotab
  return($outfa,$tblout,$annot,$notr,$nocds) if( -s $outfa or $skiprun); # or dryrun ..

  ## maybe option: $REGEN_AA
  (my $outaa=$outfa) = s/\.fsa/.faa/;
  
  my($nnamed,$namin)= parse_genenames($genenames) if($genenames);
  
  my ($inh,$outh,$tblh,$annoth,$hd,$oid,$fa,$ok);
  ## cdnaseq eq stdin also?
  if($cdnaseq =~ /\.gz$/) { $ok= open($inh,"gunzip -c $cdnaseq|"); }
  else { $ok= open($inh,$cdnaseq); }  
  $ok= open($outh,'>',$outfa) if($ok);
  $ok= open($tblh,'>',$tblout) if($ok);
  $ok= open($annoth,'>',$annot) if($ok);
  # $ok= open($outaah,'>',$outaa) if($outaa and $ok);
  unless($ok) { loggit(1,"ERR: trprocess $cdnaseq TO $outfa,$tblout"); return; }
  
# REVISE trprocess/putseq to NOT do vecscreen or gaptrim .. now in asmrna_trimvec
  # replace >oid >pubid; but need to keep oid for other input tables.
  # add gene names?
  my($itr,$otr,$ocds,$oerr)= (0) x 10;
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; 
      ## add: $outaah
      ($otr,$ocds,$oerr)= putseq($outh,$tblh,$annoth,$oid,$hd,$fa,$itr) if($fa); 
      $notr+= $otr; $nocds+= $ocds;
      $oid=$d; $hd=$_; chomp($hd); $fa=""; $itr++;
      }
    elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
  ## add: $outaah
  ($otr,$ocds,$oerr)= putseq($outh,$tblh,$annoth,$oid,$hd,$fa,$itr) if($fa); # last
  $notr+= $otr; $nocds+= $ocds;
  close($inh); close($outh);
  return($outfa,$tblout,$annot,$notr,$nocds); 
}


=item tr2main

 aabugs4/aabugs4qual/tsaevg/daphmag5_evgt2c_2013mar7.txt
 ## .. this main{md} isnt right either, as alt-2ndid can be other alt.

  grep okay $pt.trclass | perl -ne'($td,$ok,$cl,$md,$piad,$aq,$fl)=split; $cl=~s/a2$//; \
  ($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd); if($cl =~ /main|noclass/) { $main{$td}=$cl; } \
  else { $alt{$md}{$td}= $cl; $balt{$td}=$md; } $n++; \
  END { @amain= grep { not $main{$_} } sort keys %alt; \
  for $am (@amain) { $md= $balt{$am}; \
  if($main{$md}) { @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_} } @at; } \
  elsif($md) {  $main{$am}="NOMAIN";  } } \
  foreach $md (sort keys %main) { @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}};\
  $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; $mc=$main{$md}; print join("\t",$md,$mc,$ad)."\n"; } }'\
    > okayset/$pt.mainalt.tab

=cut

sub trclass2maintab
{
  my($trclass)=@_;
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  my $ntr=0;
  my $mainindex= $pubidnum_start;
  return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..
  
  my(%main,%alt,%balt,%drop,$outh,$outpubidh,$inh);
  ## my $ok= open($inh,"grep okay $trclass |");
  my $ok= open($inh,$trclass);
  $ok= open($outh,'>',$maintab) if($ok);
  $ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  ## FIXME: only althi are reliably locus alternates; altmid .. are more likely paralogs
  while(<$inh>) {
    my($td,$ok,$cl,$md,$piad,$aq,$fl)=split;
    unless($ok eq 'okay') { $drop{$td}=1; next; } # OPTION: include drops?
    $cl=~s/a2$//;  #$n++; 
    my($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd); 
    if($cl =~ /^main|^noclass/) { $main{$td}=$cl;  $balt{$td}=$td; } 
    else { $alt{$md}{$td}= $cl; $balt{$td}=$md; }  
  }

  ## Fix MISSING main links, from alt to other alts ..  
  ## FIXME2: some of these NOMAIN are drops *** dont retain;
  ## FIXME3: adding drop{xxx} has screwed up alt-links somewhere; ** STILL MESSED UP
  # ... use drop{xx} only on output?
  
  my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; 
  foreach my $am (@amain) { 
    # next if($drop{$am});
    my $md= $balt{$am} || $am; ## $md=$am if($drop{$md});
    # next if($drop{$am} or $drop{$md}); #??
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) {  $main{$am}="NOMAIN";  } # "nomain" ?
  }
  foreach my $td (keys %balt) {
    # next if($drop{$td});
    my $md= $balt{$td} || $td; 
    # $md=$td if($drop{$md}); # what ??
    $main{$md}="NOMAIN" unless($main{$md});
  }
     
  ## add headers to these:
  #originalID     MainClass  Alternates
  #Public_mRNA_ID         originalID      PublicGeneID    AltNum
  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum))."\n"
    if($outpubidh);

  #above# my $mainindex= $pubidnum_start;
  foreach my $md (sort keys %main) { 
    my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}}; # not here?  grep { ! $drop{$_} } 
    my $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc=$main{$md}; 
    
    print $outh join("\t",$md,$mc,$ad)."\n"; 
    if($outpubidh) { # should be required ??
      $mainindex++;
      my $ialt= 0;
      unless($drop{$md}) {
      my ($pubmrnaid,$pubgeneid)= make_pubid($md, $mainindex, ++$ialt);
      print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$ialt)."\n"; $ntr++;
      }
      foreach my $ad (@ad) {
        unless($drop{$ad}) {
        my ($altmrnaid,$altgeneid)= make_pubid($ad, $mainindex, ++$ialt);
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$ialt)."\n"; $ntr++;
        }
      }
    }
  }
  close($inh); close($outh);
  return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}

sub make_pubid
{
  my($oid, $mainindex, $altnum)= @_;

  ## use/check oid? keep global hash for pubid <=> oid ?
  ## my $pubid_format = $IDPREFIX.'%06d'; # $public_options{'publicid'} || "evgr000000";
  ## my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";

  # $pubidnum_start++ if($altnum == 1); ## or $mainindex == last index?
  # my $pubidnum= $pubidnum_start;  # ONLY if altnum == 1 ? or if not seen this case..
  my $pubidnum= $mainindex;
  $pubidnum_start= $pubidnum; #?
  
  my $pubgene = sprintf( $pubid_format, $pubidnum); 
  my $pubid   = $pubgene . sprintf( $altid_format, $altnum);
  return($pubid,$pubgene);
}


=item putseq

  ## putseq() from asmrna2ncbitsa.pl
  Note: evg mrna.tr may have basic cds annot info
  okayset/dmag5xau13c2011_okmrna.tr.gz
  >sodmag4nalk25loc591t24 type=cdna; aalen=906,76%,complete; clen=3562;  strand=+; offs=93-2813;
  >sodmag4nalk25loc2322t3 type=cdna; aalen=812,86%,complete; clen=2827;  strand=+; offs=203-2641;

=cut

  ## see evigene2genbanktbl.pl:putCDSloc()
sub putCDSloc 
{
  my($cdsoff,$partial,$cdsphase)= @_;
  
  ## is cdsoff == "<123->456" allowed here?
  my($start,$stop)= split/[-]/,$cdsoff; # or $cdsoff =~ m/(\d+)-(\d+)/; # 
  
  my($p5,$p3,$codonstart)= ("<",">",0);
  if($partial =~ /complete/) { }
  else {
    unless($partial =~ /partial[53]/) { $partial.="53"; } # both
    if($partial =~ /3/) { $stop="$p3$stop"; }
    if($partial =~ /5/) { $start="$p5$start";  
      $codonstart=$cdsphase+1; # is this right?
    }
  }
  my $tbl= join("\t",$start,$stop,"CDS\n");
  $tbl .= "\t\t\tcodon_start\t$codonstart\n" if($codonstart>0);
  return $tbl;
}

sub putseq 
{
  ## add param: $outaah
  my ($outh, $tblh, $annoth, $oid, $hdr, $fa, $itr)=@_; 
  my $tlog="";  
  
  # tblinfo: my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)= @$rinfo;
  # maybe write all tblinfo as separate table, more info than .tbl .. put in logfile instead of below crap?

  my $ol= length($fa); 
  my($pubid,$def,$tblinfo);
  my($cdsb,$cdse,$cdsphase,$aafull,$aastart,$aastop,$ntrout,$ncdsout)=(0) x 10;

  $tblinfo= parse_evgheader($oid,$hdr,length($fa));
  
  $pubid= $tblinfo->{'pubid'};
  $def= $pubid; # only this for ncbisubmit.fsa ?
  
  my $cdsoff= $tblinfo->{'cdsoff'}; ($cdsb,$cdse)= split/[-]/,$cdsoff; 
  my $aq= $tblinfo->{'aaqual'};    ($aafull)= $aq =~ m/(complete|partial\w*)/; 
  $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1;
  $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;
  
#  #..........................
# REVISE trprocess/putseq to NOT do vecscreen or gaptrim .. now in asmrna_trimvec
  
  ## FIXME: vectrim in CDS should be disallowed *@#($ causing loss of stop codons in valid/hi-homology cds
  ## FIXME2: use vtype vt, skip all but Strong if in cds span, regardless of aastart/stop
  ## Ooops : cdse is last base-3 of stopcodon, not first
  my($vectrimw)=(0,0);
  if(my $vec=$vecscreen{$oid}) { 
    # FIXMEd: multi locs per vd; sep by \n
    my @vec=split"\n",$vec; my $nv= @vec;
    foreach $vec (@vec) {
      my($vb,$ve,$vtype)=split"\t",$vec; 
      my $overcds= ($vb <= $cdse and $ve >= $cdsb)?1:0;
      if($overcds) {
        next unless($vtype =~ /Strong/i);
        if($aastop  and $ve >= $cdse and $vb <= $cdse+1) { $vb=1+$cdse; } # dont trim stop codon !!!
        if($aastart and $vb <= $cdsb+2 and $ve >= $cdsb) { $ve=$cdsb-1; } # dont trim start codon !!!
      }
      next if($vb >= $ol or $ve<=$vb);
      my $trimw=1+$ve-$vb; $vectrimw += $trimw;
      substr($fa,$vb-1,$trimw)= ("N") x $trimw; 
      }
    } 

  # geneinfo for ncbi tsa submit: move this to feature.tbl not cdna.fsa; add cds offsets
  # >Feature Thecc1RA_L_g13025t00001
  #  bg  eg  gene  locus_tag TCM_000002
  #  bc  ec  CDS   product Cystathionine beta-synthase (CBS) family protein isoform 1
  #                protein_id  gnl|CacaoGD|Thecc1EG000002p1
  
  #....... GAP Cleaning is a messy case .............
  my $nn= $fa =~ tr/Nn/Nn/; 
  my $ncut=0;
  my $CDShasTooManyXs=0;
  if($nn>0) { 
    my ($lcdsb,$lcdse)= ($cdsb,$cdse);
    $fa =~ s/n/N/g;
    
    # FIXME5: SEQ_FEAT.CDShasTooManyXs : need to count NNN inside CDS and disallow a few soaptr gappy cds; more gaps than cds   
    my $cdsfa= substr($fa,$cdsb-1,1+$cdse-$cdsb); ## cdsb was OFF-By-1
    my $cdsnn= $cdsfa =~ tr/N/N/; 
    if($cdsnn > 0.48*length($cdsfa)) { # ERR
      $CDShasTooManyXs= $cdsnn; # flag it for below .. uniqname ?? check below
    }

    # FIXME4: No NN End trim in CDS for aastart,aastop: keep stop/start codon 
    # FIXME: cdsb,cdse adjust for inner gaps
    # fixme2: must adjust cdsb,e when cut BEFORE cdsb
    # fixme3: this is a mess; better to a. cut NNN, b. rerun cdna_bestorf for new cds offset?
    ## SEQ_INST.HighNContentStretch: stretch of at least 5 Ns within the last 10 bases << CHANGED to 20 bases

    ##  [SEQ_INST.LeadingX] Sequence starts with leading X BIOSEQ: MusacuEGm045505p18 == bananavel2k29Loc2240t7
    # >bananavel2k29Loc2240t7 type=cdna; aalen=185,94%,partial5; clen=590;  strand=+; offs=3-560;
    # ANNNNNNNNNNNNNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    
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
    
    my $ne= rindex($fa,'N'); 
    my $curlen= $ol; my $iter=0;
    while($ne >= $curlen - $ENDGAP and $curlen > $ENDGAP and $iter<5) {
      $iter++;
      if($ne >= $curlen - $ENDGAP and $ne <= $cdse and $aastop) { $ne= index($fa,'N',$cdse+1); }
      if($ne >= $curlen - $ENDGAP) { 
        $fa=substr($fa,0,$ne); 
        if($fa=~s/(N+)$//) {  my $ncut=length($1); $ne-=$ncut; }
        if($ne < $cdse) { $cdse = $ne; } #??
        $curlen= length($fa); 
        $ne= rindex($fa,'N'); 
        }
    }
      
    ## FIXME: cds-phase/codon_start changes w/ mod 3 of n1   
    ## FIX2: see above  nnnAnnn < need to recurse chop endgaps ?
    my $n1= index($fa,'N');  $iter=0;
    while( $n1 >= 0 and $n1 <= $ENDGAP and $iter<5) {
      $iter++;
      if($n1>=0 and $n1 <= $ENDGAP and $n1 >= $cdsb and $aastart) { $n1= rindex($fa,'N',$cdsb-1); } # is this right? partial/5 check?
      if($n1>=0 and $n1 <= $ENDGAP) { 
        $n1++; $fa= substr($fa,$n1);  
        if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
        
        ## FIXME3: w/ recurse, cdsb < 0 often here... leads to this crap
        ##    225 ERROR: SEQ_FEAT.InternalStop; SEQ_FEAT.NoStop
        ##    236 REJECT: SEQ_FEAT.Range
        #bad.now# if($cdsb>0) ; use lcdsb or cdsoff
        if($lcdsb>0) { $cdsb -= $n1; $cdse -= $n1;  } # $cdsphase = $n1 % 3; ?? ** PROBLEMS now w/ phase change
        $n1= index($fa,'N'); 
        }
    }
      
    $ncut=0; my $gapw= length( $GAPSMAX); #== MAXGAP
    unless($GAPSOK) {
      for (my $in= index($fa,$GAPSMAX); $in >= 0; ) {
        my $w=length($fa); my $en=$in+$gapw; 
        $en++ while($en<$w and substr($fa,$en,1) eq "N"); 
        my $wn= $en-$in; my $keep= 3 + ($wn % 3); my $cut= $wn-$keep; $ncut+=$cut; 
        my $facut= substr($fa,0,$in).substr("NNNNNN",0,$keep).substr($fa,$en); 
        $fa=$facut; 
        if($cdse>0) {
          if($en < $cdsb) { $cdsb -= $cut; $cdse -= $cut; } ##  $cdsphase = $cut % 3;
          elsif($in < $cdse and $en > $cdsb) {
            if($in <= $cdsb) { $cdsb -= $cut; $cdse -= $cut; } #??  $cdsphase = $n1 % 3; ??
            else { $cdse -= $cut; }
          }
        }
        $in=index($fa,$GAPSMAX); 
      } 
    }

    ## FIXME: cant have neg cdsb ..
    #er2g: >dmag5xevgr001990t1      oid=dmag4vel4xbxk75Loc4074t1    len=607; olen=615; cut=0; nnn=12/13;    cutcds=-5-501,oldcds=3-509,
    unless($cdse == $lcdse and $cdsb == $lcdsb) {
      if($cdsb <= 0) { 
        $cdsphase = (2 + $cdsb) % 3; # ?? this seems right, dont know why !!
        $cdsb=1;  my $aafull0= $aafull;
        $aafull="partial5" if($aafull=~/complete/); # 
        $aafull="partial"  if($aafull=~/partial3/); # 
        $tblinfo->{'aaqual'} =~ s/$aafull0/$aafull/ if($aafull0 ne $aafull);
        #FIX: $aaqual =~ s/complete|partial3/partial/; 
        }
      $tblinfo->{'cdsold'}= $tblinfo->{'cdsoff'}; # for annotab?
      $tblinfo->{'cdsoff'}= "$cdsb-$cdse";
      $tlog.="cutcds=$cdsb-$cdse,oldcds=$lcdsb-$lcdse,"; 
    }
  } 
  
  my $nl= length($fa); 
  $tblinfo->{'trlen'}= $nl;
  my $nn1= $fa=~tr/N/N/; 
  my $lendelta= ($nl==$ol and $nn==0) ? $nl : "$nl; olen=$ol; nnn=$nn1/$nn;";
  $lendelta .= " cut=$ncut;" if($ncut>0); 
  $tlog.="vectrim=$vectrimw," if($vectrimw); 
  ##er2g: >PogonEG0028883t1        oid=sobeetlepogo1ak31loc13234t4 
  ## len=413; olen=436; cut=0; nnn=18/41;    cutcds=1-314,oldcds=2-337,vectrim=23,
  ## for annotab:  cdsoff = cutcds; oldcds,olen not needed;  add: nnn=xxx,vectrim=xxx
  my $annogaps= ($nl==$ol and $nn==0) ? "0" : "gaps=$nn1/$nn";
  $annogaps.= ",oldcds=".$tblinfo->{'cdsold'} if($tblinfo->{'cdsold'});
  $annogaps.= ",cut=$ncut" if($ncut>0); 
  $annogaps.= ",vectrim=$vectrimw" if($vectrimw); 
  $annogaps = "TOOSHORT=$nl,$annogaps,cdsnn$CDShasTooManyXs" if($nl<$MINSIZE or $CDShasTooManyXs);
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct)};
  my $cddname= $tblinfo->{'cdd'} || "na";

  ## FIXME: always annotate CDD name if exists.
  if($annoth) { 
    # usable output table ; FIXMEd: add header at top
    # NOTE: add vectrim info, other?; will need to regen .aa, .cds from .fsa output to be accurate..
    # .. user choice: ncbi submit restrictions may not be desired.
    # FIXmaybe: recode dbxref tags: as for .tbl ?? UniProt => SwissProt, xxxx:ENS => ENSEMBL:ENS ..
    
    my $aname= $gname || "hypothetical protein"; # || "na";  # which ??
    print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name))."\n" if($itr==1);
    print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$aname,$cddname)."\n";
  }
  
  if($nl < $MINSIZE or $CDShasTooManyXs) { #? do AFTER annoth so table has full record of skips?
    my $errtype= ($CDShasTooManyXs) ? "CDShasTooManyXs:$CDShasTooManyXs": "too short:$nl";
    my $uniqname=0;  ## retain shorties w/ uniq name **
    if($gname) { 
      my $allids= $namedgenes{$gname}; 
      $uniqname=1 if($allids eq "$oid,");
    }
    my $logt="#$def\toid=$oid\tlen=$lendelta";
    $logt.= "\tname=$gname" if($gname); $logt.= "\t$tlog" if($tlog);
    if($uniqname) {
      loggit(0,"PROBLEM: keep unique name but $errtype",$logt);
    } else {
      loggit(0,"ERR: skip $errtype",$logt);
      return (0,0,"ERR: $errtype");
    }
  }
        
  map{ $_="" if($_ eq "na"); } ($lotag,$gname,$cddname,$dbxref,$aaqual);  

  $aafull ||= $aaqual;
  my $hascddname= ($cddname and $cddname ne $gname)?1:0;
  $gname .= " fragment" if($aafull =~ /partial/ and $gname);

  if($cdsoff =~ m/\d+\-/) { # allow for <123->456
    $ncdsout++; #? always same count as ntrout ?
    print $tblh ">Features\t$pubid\t$GBPROID\n\n"; #?? is this used now
    # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827
    
    my $cdsline= putCDSloc($cdsoff,$aafull,$cdsphase);
    print $tblh $cdsline;
    # $cdsoff =~ s/\-/\t/; print $tblh "$cdsoff\tCDS\n";
    
    print $tblh "\t\t\tprotein_id\t$protid\n" if($protid);
    print $tblh "\t\t\tlocus_tag\t$lotag\n" if($lotag);
    (my $npct=$namepct) =~ s/,.*//;
    print $tblh "\t\t\tnote\talignment:blastp is $npct\n" if($gname and $npct>0);
    print $tblh "\t\t\tnote\toriginal_id:$oid\n" if($oid);
    my $pname= $gname || "hypothetical protein";
    print $tblh "\t\t\tproduct\t$pname\n";
#         # ^^ FIXME: SEQ_FEAT.MissingCDSproduct; need some name
    print $tblh "\t\t\tnote\tconserved domain $cddname\n" if($hascddname);

    if($dbxref =~ /\w/) {
      foreach my $dx (split",",$dbxref) { 
        # DBXREF_RECODE fix for ncbi's  dbx: restrictions
        $dx =~ s/^UniRef/SwissProt:UniRef/; ## UniProtKB:UniRef/; # or /TrEMBL:UniRef/;  SwissProt
        $dx =~ s/UniProt:/SwissProt:/; ## UniProtKB: # need list to check/replace as per 
        # is CDD:id ok? accepted by tbl2asn .. YES
        
        ## ==> catfish1all4cf/okayset/catfish1all4.mrna.tbl2asn.sumval <==
        ## 60800 WARNING: SEQ_FEAT.IllegalDbXref
        ##  db_xref DRERI:ENSDARG00000019365 ## >> ENSEMBL:  FIXME in where? catfish/zfish.names
        ## anything:ENS\w+\d\d+ should become ENSEMBL:
        if($dx =~ /:ENS/ and $dx !~ /^ENSEMBL:/) { $dx =~ s/^\w+:(ENS\w+\d\d+)$/ENSEMBL:$1/; }
        ## FIX2: arath:AT5G60040.1 << should be TAIR: dammit ; fix .names instead of here..
        
        # evigene2genbanktbl.pl sub reformatval()
        #  my($d)= m/^(\w+):/; if( $DBXREF_RECODE{$d} ) { $d= $DBXREF_RECODE{$d}; }
        
        print $tblh "\t\t\tdb_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); 
      }
    }
    print $tblh "\n";
  }
  
  $fa =~ s/(.{60})/$1\n/g; 
  print $outh ">$def\n$fa\n";  $ntrout++;
  loggit(0,">$def\toid=$oid\tlen=$lendelta\t$tlog"); # \tgene=$gn
  
  ## Maybe add option regen .aa at least from updated fa, cdsoff
  # if($outaah) { # == $REGEN_AA 
  #   my($XXaalen,$XXpcds,$XXcompl,$XXorflen,$fahead,$faprot) = cdna_proteins:translate1($fa,$cdsb,$cdse,AS_FASTA=1,...);
  #   print $outaah ">",$def," $fahead\n",$faprot,"\n";
  # }
  
  return ($ntrout,$ncdsout,"OK"); # no, return ($ntrout,$ncdsout)
}





## move to rnaseq/asmrna_trimvec.pl
sub vecscreen
{
  my($cdnaseq,$vectab,$skiprun)=@_;
  $vectab= makename($cdnaseq,".vector.tab") unless($vectab);
  return($vectab) if( -s $vectab or $skiprun); # or dryrun ..
  
  our($id,$vb,$ve,$ty,$vd,$outh,$inh);  
  sub putv { our($id,$vb,$ve,$ty,$vd,$outh,$inh); 
    print $outh join("\t",$id,$vb,$ve,$ty,$vd)."\n" if($id and $ty and not $ty=~/Weak/); 
  }

  sub putv2 { 
    my($outh,$id,$ty,$vd,@vbe)= @_; 
    return unless($id and $ty and not $ty=~/Weak/i);
    foreach my $vbe (@vbe) { my($b,$e)= @$vbe;
      print $outh join("\t",$id,$b,$e,$ty,$vd)."\n"; 
    }
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
  
  my $ok= open($inh,$vectmp); 
  $ok= open($outh,'>',$vectab) if($ok);
  unless($ok) { loggit(1,"ERR: $APPvecscreen -i $cdnaseq -d $univecdb TO $vectab"); return; }
  my($so, @vbe);
  while(<$inh>) {
    chomp; 
    if(/^>/) { 
      # putv() if($id); 
      putv2($outh,$id,$ty,$vd,@vbe) if($id);
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
  # putv() if($id); 
  putv2($outh,$id,$ty,$vd,@vbe) if($id);
  close($outh); close($inh);
  return($vectab);
}


sub vecscreen_ncpu
{
  my($npart,$cmd0,$ccount,$cdnaseq,$vecout)=@_;
  
  # my $ccount= facount($cdnaseq); # use this, not fasize
  my $splcount= int(0.99 + $ccount/$npart);
  my $spldir= makename($cdnaseq,"_vecsplit/",'cdna|fasta|fsa|fa');  # use _tsasubmit/ instead?
  mkdir($spldir); # dryrun?
  
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
  open(my $outh,'>',$vecout);
  foreach my $vf (@vecset) { if(open(FT,$vf)) { while(<FT>) { print $outh $_; } close(FT);} }
  close($outh);
  return($err); # $vecout
}


## move to cdna_evigenesub.pm
# sub revcomp { # UNUSED here
#   my ($seq) = @_;
#   my $reversed_seq = reverse ($seq);
#   $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
#   return ($reversed_seq);
# }

sub facount {
  my($infile)=@_;
  my $nrec=0;
  if(-s $infile) { 
    if($infile=~/\.gz$/) { $nrec=`gunzip -c $infile | grep -c '^>'`; }
    else { $nrec=`grep -c '^>' $infile`; }
    chomp($nrec); }
  return $nrec; # what? "nrec=".
}

# sub fasize { # UNUSED here
#   my $fa=shift; my $b=0; my $ok;
#   #if($fa=~/\.gz$/) { $b=`gunzip -c $fa | grep -v '^>' | wc -c | sed 's/ .*//'`; } 
#   #else { $b=`grep -v '^>' $fa | wc -c | sed 's/ .*//'`; }
#   if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
#   else { $ok= open(F,$fa); }
#   while(<F>) { $b += length($_) unless(/^>/); } close(F);
#   chomp($b); return $b; 
# }

# sub fasplit { # UNUSED here; this splits by size; need alternate split by count: facount/ncpu
#   my($fa, $spldir, $npart, $splsize, $fasuf)=@_;
#   $fasuf ||= "fa";
#   my @splist= (); my $ok= 0;
#   my($atsize,$atfile)=(0,0);
#   my $fabase= basename($fa);
#   if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
#   else { $ok= open(F,$fa); }
#   unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splsize"); return @splist; }
#   mkdir($spldir) unless(-d $spldir);
#   while(<F>) {
#     if($atfile==0 or ($atsize > $splsize && /^>/)) {
#       $atfile++; if($atfile>$npart) { } # finishup???
#       close(SPL) if($atfile>1);
#       my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
#       $ok= open(SPL,'>',$spname);  $atsize= 0;
#       unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
#       push @splist, $spname;
#       }
#     print SPL;
#     $atsize+= length($_) unless(/^>/);
#   } 
#   close(SPL); close(F);
#   return @splist;
# }

sub fasplitcount { # this splits by size; need alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splcount, $fasuf)=@_;
#  $cdnaseq, $spldir, $npart, $splcount,"fsa"
  $fasuf ||= "fa";
  my @splist= ();
  my($ok, $irec,$atfile)=(0,0, 0);
  my $fabase= basename($fa); #? is this bad
  if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
  else { $ok= open(F,$fa); }
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splcount"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  
  while(<F>) {
    if($atfile==0 or ($irec >= $splcount && /^>/)) {
      $atfile++;  $irec=0;
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname); 
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $irec++ if(/^>/);
    ## $atsize+= length($_) unless(/^>/);
  } 
  close(SPL); close(F);
  return @splist;
}

sub tsa_tbl2asn
{
  my($cdnaseq,$cdnatbl,$organism,$sraids)=@_;
  return unless(-s $cdnaseq);
  
  ##.. this file part input works..
  # pt=litova1all3.mrnatop2000
  # $ncbin/tbl2asn -i tsa.$pt.fsa -f tsa.$pt.tbl -Z discrep.$pt.log \
  # -a r10k -l paired-ends -Vtb -Mt -XE \
  # -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt \
  # -j "[moltype=mRNA] [tech=TSA] $org $sra"
  ##  -r  Path for Results [String]  Optional
  
  # my $wdir= dirname($cdnaseq);
  # my $cdnabase= basename($cdnaseq);
  
  my $spldir= makename($cdnaseq,"_tsasubmit/",'cdna|fasta|fsa|fa');  # ok for _split dir
  
  ## -Vb = gen genbank gbf, not needed, option
  my $tsaopts="-a r10k -l paired-ends -Vt -Mt -XE"; # read from config.
  my $tsadesc= $TSADESC; # "-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
    ## ^^ need config for these; generate some of this?
    ## FIXME: tbl2asn dies silently if missing files -w .. -Y .. or -t ..
    ## UGH: qsub gave this via env sra='xxx':
    ##   -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']'
    ## egmrna2tsa.19433.err sh: -c: line 0: unexpected EOF while looking for matching `''
#er2g: forkCMD= /home/ux455375/bio/ncbi2227/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE \
#  -w ./pogonus1all3.tsamethods.cmt -Y ./pogonus1all3.tsadesc.cmt -t ./pogonus1all3.tsasubmit.sbt \
# -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']' \
# -i ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.fsa \
# -f ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.tbl \
# -Z ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.discrep

  # map{ s/^['"]+//; s/['"]+$//; } ($organism, $sraids); ## $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; ??
  map{ s/^\W+//; s/\W+$//; } ($organism, $sraids); ## $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; ??
  my $tsaqual="-j \'[moltype=mRNA] [tech=TSA] [organism=$organism] [SRA=$sraids]\'";
  my $cmd0="$APPtbl2asn $tsaopts $tsadesc $tsaqual "; # add in loop: -i $pt.fsa -f $pt.tbl -Z discrep.$pt.log
  
  my ($cmddone, $nout, $sqnout)=(0) x 10;
  my $runok= ($APPtbl2asn =~ /MISSING/) ? 0 : 1;

  #?? check for tbl2asn file set?  
  # $sqnout=$cdnaseq;  $sqnout=~s/\.fsa/.sqn/; 
  $sqnout= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa');   
  return(1,"",$sqnout) if( -s $sqnout ); ## all sqnout have full path ?
  if( -d $spldir ) { ## ! -s $sqnout and 
    opendir(D,$spldir); 
    my @sqnout= map{ "$spldir/$_" } grep /\.sqn/, readdir(D); 
    closedir(D);
    $sqnout= join ", ", @sqnout; $nout= @sqnout;
    return($nout,$spldir,$sqnout) if($nout>0); 
  }
  
  $sqnout="";
  mkdir($spldir); # dryrun? only ncpu or use for all cases?
  if($runok and $NCPU > 1) {
    my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
    if($ccount >= 50*$NCPU) {
      ($nout,$sqnout)= tbl2asn_ncpu($NCPU,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl);
      $cmddone=1;
    }
  } 
  
    # put into $spldir ?
  if($runok and ! $cmddone) {
    my $dlog1=$cdnatbl; $dlog1 =~ s/\.\w+$//; $dlog1.=".discrep";
    my $cmd1= $cmd0 . " -i $cdnaseq -f $cdnatbl -Z $dlog1"; # want -r $spldir or not ?
    my $err= runcmd($cmd1);    
    my $sqnt= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa'); # add $spldir ??
    if( -s $sqnt ) { $nout=1; $sqnout=$sqnt; }
    # $nout= ($err==0 and $runok)?1:0;
    # if($nout>0) { $sqnout=$cdnaseq;  $sqnout=~s/\.fsa/.sqn/; }# makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa'); 
  }
  return($nout,$spldir,$sqnout); # list all?
}


sub tbl2asn_ncpu
{
  my($npart,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl)=@_;
  # $NCPU,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl
  
  # my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
  my $splcount= int(0.99 + $ccount/$npart);
  # my $spldir= makename($cdnaseq,"_split/");  # use _tsasubmit/ instead?
  # mkdir($spldir); # dryrun?
  
  ## NOTE: fasplit adds spldir to set path
  my @splset= fasplitcount( $cdnaseq, $spldir, $npart, $splcount,"fsa"); 
  my @tblset= fasplitcount( $cdnatbl, $spldir, $npart, $splcount,"tbl");  # need to split into same parts.
  my $npartgot= @splset;
  my $icpu= 0;  
  # NO: chdir($spldir); # bad!, need files in starting path ... use -i spldir/in -r spldir/
  for(my $ip=0; $ip< $npartgot; $ip++) {
    my $cdna1= $splset[$ip];
    my $tbl1 = $tblset[$ip];
    (my $dlog1=$tbl1) =~ s/\.\w+$/.discrep/; # $dlog1= makename($cdna1,".discrep");
    my $cmd1= $cmd0 . " -i $cdna1 -f $tbl1 -Z $dlog1"; 
      ## DONT need -r spldir ? w/o put parts in same dir as -i, -f
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npartgot) { while (wait() != -1) { }; $icpu= 0; }
  }
  while (wait() != -1) { };
  
  ## One part problem, errorsummary.val parts overwrite .. rebuild from $name.val ? or
  # my @valout= # my @gbfout= 
  # my @sqnout= map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset;
  
  ## deal with $name.split*.fixedproducts also?
  ## deal with $name.split*.discrep also?
  
  ## remake errorsummary.val from parts  .. maybe for 1cpu also? count ERROR: and loggit ?
  my %et;
  my $esumf= makename($cdnaseq,".tbl2asn.sumval"); # or cdnaseq.errorsummary.val ?
  foreach my $sf (@splset) { 
    (my $valf=$sf) =~ s/\.fsa/.val/;
    if(open(VF,$valf)) { 
      while(<VF>) { my($err,$typ)= m/^(\S+)[\w\s]+.([\w\.]+)/; $et{$err}{$typ}++ if($typ); } close(VF); 
    }
  }
  open(my $esumh,'>',$esumf);
  foreach my $e (sort keys %et) { my @t=sort keys %{$et{$e}}; 
    foreach my $t (@t) { my $c=$et{$e}{$t}; printf $esumh "%6d %s %s\n",$c,$e,$t; } 
  } close($esumh); 

  ## check existance -s sqn     # or readdir(D) as above  
  my @sqnout= grep { -s $_ } map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset; 
  my $sqnout= join ", ", @sqnout;
  my $nsqn= @sqnout;
  return($nsqn,$sqnout); # list all?
}

=item tbl2asn log

#er2g: app=tbl2asn, path=/home/gilbertd/bin/tbl2asn
#er2g: evigeneapp=prot/traa2cds.pl, path=/bio/bio-grid/mb/evigene/scripts/prot/traa2cds.pl
#er2g: get_evgtrset= ./okayset/allstrimpt1.mrna . allstrimpt1
#er2g: trclass2maintab primary n=  allntr=  allstrimpt1.pubids
#er2g: vectors found in ntr= 962 ./okayset/allstrimpt1.mrna.vector.tab
#er2g: DONE output ntr=0, ncds=0 in files allstrimpt1.mainalt.tab, allstrimpt1.pubids, ./okayset/allstrimpt1.mrna.fsa, ./okayset/allstrimpt1.mrna.tbl, ./okayset/allstrimpt1.mrna.ann.txt
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE \
 -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt \
 -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]' \
 -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.fsa \
 -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.tbl \
 -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.discrep
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]'  -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.fsa -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.tbl -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.discrep
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]'  -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.fsa -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.tbl -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.discrep
#er2g: DONE tsa_tbl2asn nparts=3, submitset=./okayset/allstrimpt1.mrna_tsasubmit/  ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.sqn, ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.sqn, ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.sqn

remake errorsummary.val

cat okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split?.val | perl -ne\
'($err,$typ)= m/^(\S+)[\w\s]+.([\w\.]+)/; $et{$err}{$typ}++; END{ foreach $e (sort keys %et) { @t=sort keys %{$et{$e}}; foreach $t (@t) { $c=$et{$e}{$t}; printf "%6d %s %s\n",$c,$e,$t; } } }' | head

     1 ERROR: SEQ_FEAT.NoStop
     2 WARNING: SEQ_FEAT.CDShasTooManyXs
     2 WARNING: SEQ_FEAT.EcNumberProblem
 33777 WARNING: SEQ_FEAT.PartialProblem
     8 WARNING: SEQ_INST.HighNContentStretch


=item tbl2asn work

# ## tbl2asn outputs:  .sqn == asn1, can it be catted ??? NO.. need submit each part file?
# ##  .val == errors list, can be catted
# ##  .gbf if made can be catted
# # my $cmd;
# # $cmd= "cat ".join(' ',@aaset)." > $aaseq";   runcmd($cmd);
# # $cmd= "cat ".join(' ',@cdsset)." > $cdsseq"; runcmd($cmd);
## keep as part files for tsa submit for now
#   if($debug||$dryrun) {
#     push @erasefiles, @splset, @aaset, @cdsset; # which?
#   } else {
#     foreach my $fn (@splset, @aaset, @cdsset) {  unlink($fn) if(-f $fn); } 
#     rmdir($spldir); 
#   }
  
=cut


## in cdna_evigenesub.pm
sub findapp
{
  my($aname, $nodie)=@_;
  my $app=""; $nodie ||= 0;
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=`which $aname` unless($app); 
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$debug||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "app=$aname, path=$app");
  return($app);
}

## in cdna_evigenesub.pm
sub findevigeneapp
{
  my($aname, $nodie)=@_;
  my $app= $aname; $nodie ||= 0;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$debug||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}



## in cdna_evigenesub.pm
sub runcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"CMD=",@cmd);  
  ## fail if $cmd[0] =~ /MISSING_/
  my $err= ($dryrun) ? 0 : system(@cmd);  
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

## in cdna_evigenesub.pm
sub forkcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"forkCMD=",@cmd);  
  unless($dryrun) {
    my $pid= fork();
    if($pid) { # parent
      return $pid;
    } else { # child
      my $err= system(@cmd);
      exit($err);
    }
  }
  # if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
}

## in cdna_evigenesub.pm
sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';  ## fixme need insuf: tr|fasta|fa
  my $outfile= $infile; $outfile =~ s/\.gz$//;
  $outfile =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}


__END__
 
 
=item tsasub11 cacao

  evigene=/Users/gilbertd/Desktop/dspp-work/genes2/evigene
  ncbin=/Users/gilbertd/Desktop/dspp-work/genomesoft//ncbi/bin
  nwbsra='[SRA=SRR531454; SRR531455; SRA058778; SRA058779]'
  vel3sra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782; SRR531454; SRR531455; SRA058778; SRA058779]'
  rnasra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782]'
    ## redo asmrna2ncbitsa.pl -GAPSOK
    
  cd tsasub11
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
  
  if [ -f $pt/TCM01.tsa_rasm.$pt.tbl ]; then continue; fi
  echo asmrna2ncbitsa $pt/TCM01.tsa_rasm.$pt
  
  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl
  
  } done
  
  
    # was -a r10u; unknown gap len; YES, r10k got rid of warning
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
   if [ $pt == cacao3nwb ]; then sra=$nwbsra; 
   elif [ $pt == cacao3vel ]; then sra=$vel3sra; 
   else sra=$rnasra; fi
   
   if [ -f $pt/TCM01.tsa_rasm.$pt.sqn ]; then continue; fi
  
   echo $pt : $sra
   
   $ncbin/tbl2asn -p $pt/ -Z $pt/$pt-discrep.log  -w tsamethods.$pt.cmt \
   -Y tsadesc.$pt.cmt -t cacao3i_tsasubmit.sbt \
   -a r10k  -l paired-ends -Vtb -Mt -XE \
   -j "$sra [bioproject=PRJNA51633] [moltype=mRNA] [tech=TSA] [organism=Theobroma cacao]"
  
  } done
  
=cut

=item run_evgmrna2tsa.sh

  #! /bin/bash
  ### env idprefix=MysppEGm trclass=myspecies.trclass datad=`pwd`  qsub -q shared run_evgmrna2tsa.sh
  #PBS -N evgmrna2tsa 
  #PBS -l nodes=1:ppn=8,walltime=5:55:00
  #PBS -o egmrna2tsa.$$.out
  #PBS -e egmrna2tsa.$$.err
  #PBS -V
  
  ncpu=8; # most 100K trsets run on 8cpu in 10m-20min; 300K set can take 1hr.

  evigene=$HOME/bio/evigene/scripts
  export vecscreen=$HOME/bio/ncbic11/bin/vecscreen
  export tbl2asn=$HOME/bio/ncbi2227/bin/tbl2asn
  
  ## FIXME: -TSADESC=tbl2asn flags for path/to/tsa.cmt files
  ## default $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt";
  
  opts="-debug -runtbl2asn -NCPU $ncpu"
  if [ "X" = "X$datad" ]; then echo "missing env datad=path/to/data"; exit -1; fi
  if [ "X" = "X$trclass" ]; then "echo env trclass=path/to/name.trclass"; exit -1; fi
  ## .. these are now read via sra_result.csv, species => idprefix
  if [ "X" != "X$species" ]; then opts="$opts -species=\'$species\'"; fi
  if [ "X" != "X$sra" ]; then opts="$opts -sraids=\'$sra\'"; fi
  if [ "X" != "X$idprefix" ]; then opts="$opts -idprefix $idprefix"; fi
  
  cd $datad/
  ## add -outdir opt or : mkdir tsasubmit; cd tsasubmit
  ##FIXME: template files ; need these  or generate defaults?
  if [ ! -f evgr_tsasubmit.sbt ]; then
    if [ -d ../tsasubmit ]; then cp ../tsasubmit/*.{cmt,sbt} ./; fi
  fi
    
  echo $evigene/evgmrna2tsa.pl  $opts -log -class $trclass
  $evigene/evgmrna2tsa.pl  $opts -log -class $trclass

=cut 