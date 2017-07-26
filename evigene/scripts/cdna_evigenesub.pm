# cdna_evigenesub.pm

# package cdna_evigenesub;
package main;

use strict;
use warnings;
use FindBin;
use File::Basename qw(basename dirname fileparse);

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (xxxx);

use vars qw ( $EVIGENES $EGAPP $GDB_PREFIX $DEBUG $dryrun
  $MIN_NAMEIDENT $MIN_IDLIKE  
  %genenames %genedbxref %genenamepct %namedgenes %cddnames 
  %pubids $APPtraa2cds
  );

our $DEBUG=0;
our $dryrun=0; ## $DRYRUN ?
our $EVIGENES="$FindBin::Bin"; #??
our $EGAPP='mrna2tsa'; # FIXME?
our $EGLOG='egr';
our $GDB_PREFIX='gnl|Evigene|';  #see below;  use IDPREFIX? No, ID has
our $APPtraa2cds= undef; #findevigeneapp("prot/traa2cds.pl"); # move to cdna_evigenesub for get_mRNA

## name stuff 
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, };
our $logh= undef; # package local?
sub loggit{ 
	# my $dowarn=shift; my $s= join(' ',@_); # dang warn @_ empty join for 1st call here, from where ??
	my($dowarn,@msg)= @_; return unless($dowarn or @msg);
  my $s= join(' ',@msg); chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($logh){ print $logh "#$EGLOG: $s\n"; } elsif($dowarn>0||$DEBUG){ warn "#$EGLOG: $s\n"; }
  if($dowarn == LOG_DIE) { die "#$EGLOG: $s\n" ; }
}

sub openloggit {
  my($logfile,$trname)= @_;
  if(not $logfile and defined $logfile) { # use output name
    $logfile= $trname || $EGLOG;
    $logfile= makename($logfile,".$EGAPP.log");  # need program suffix??
  }
  if($logfile) { open($logh, '>>', $logfile) or die $logfile; } 
}

sub openRead { # add to cdna_evigenesub.pm
  my($fna)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") : open($hin,$fna);  
  loggit(1,"ERR: openRead $fna") unless($ok);
  return ($ok,$hin);
}

## note these are in cdna_protein also; need more package local privacy.
sub _min1 { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max1 { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }



sub parse_genenames
{
  my($genenames)= @_;
  my($ngot,$nin)=(0,0);
  # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
  %genenames=%genenamepct=%genedbxref=%namedgenes=%cddnames=();
  return($ngot,$nin) unless($genenames and -f $genenames);
  
  ## FIXME2: ** use uniq names vs ERR: too short to keep valid tr, e.g. 
  #er2g: ERR: too short:183 #LitvaEG0018688t4     oid=litovavel2k35Loc15824t1     len=183 name=CDD: clpS, ATP-dependent Clp protease ada..-like    
  # grep  'CDD: clpS,' *.names = 1 only = litovavel2k35Loc15824t1
  
  # FIXME: need better reader; 2+ rows/id; pick best .. format may change..
  # names.tab ==  id, name, pctalign, refid, repid  : now
  #  trid1  C-ets-2 protein, putative       89%,103/116,197 RefID:UniRef50_E0VFI2   RepID:E0VFI2_PEDHC
  #  trid2  DBH-like monooxygenase protein 1        73%,445/613,516 RefID:UniRef50_Q6UVY6   RepID:MOXD1_HUMAN
  
  # open( my $inh, $genenames) or loggit(1,"ERR: parse_genenames reading $genenames");
  my($ok,$inh)= openRead($genenames);
  unless($ok) { loggit(1,"ERR: parse_genenames reading $genenames"); return; }
  
  while(<$inh>) { 
    next unless(/^\w/ and /\t/);
    chomp; $nin++;
    my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name
    my $xtra; ($name,$xtra)=split";",$name,2; 
    $name =~ s/\s+$//;

    ## BUG in data: missing pctalign fields ; dont know why.
    ## output of evigene/scripts/prot/namegenes.pl 
    ## whitefly1vel5k45Loc9888t2       Synaptojanin-1-like protein     RefID:UniRef50_B4E1Z3   UniProt:B4E1Z3_HUMAN
    if($refid and not defined $repid and $pctalign and $pctalign =~ /^\D/) {
    	$repid=$refid; $refid=$pctalign; $pctalign="";
    }
    $refid||=""; $repid||="";
    
    # FIXME: 2 names/id maybe: CDD: xxx and gene xxx; keep both in ann.txt ? and pctalign?
    ## pctalign == 100%,450/450,446 : pct,naln/nref,ntrg
    ## old geneinfo:
    #  ($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
    ## usage below in putseq
    # my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
    #  aaqual => "na", trlen => 0, cdsoff => "", cdsor => 1, 
    #  name => $genenames{$oid}||"", dbxref => $genedbxref{$oid}|| "na"); 

    ## keep this? option? 
    my($pcta)= $pctalign=~m/(\d+)/; # warn Argument "87%,80/92,86" isn't numeric 
    if($pctalign =~/^\d/ and $pcta < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
      ## bad: Uncharacterized protein-like ; Nuclease HARBI1, putative-like; Protein-like
      if($pcta >= $MIN_IDLIKE) { unless($name =~ /\blike|^Uncharacterized/) {
        $name =~ s/, putative//; 
        unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
        }
      } else { next; } ## should we preserve for ann.txt table ?
    }

        # DBXREF_RECODE fix for ncbi's  dbx: restrictions : evgmrna2tsa2.pl:putTblFsa()
        # FIXME999: more problems w/ gene.names table having odd/local DBprefix:ID
        #   .. fix where? should have here list of valid NCBI/ISxxx db prefixes. from where?
    
    ## fixme: CDD:206692,cd04107,RefID:UniRef50_Q9NX57,UniProt:RAB20_HUMAN, 
    ## drop  RefID:; drop? cd04107
    $refid =~ s/^RefID://;  $repid =~ s/^RepID://;  ## RepID: also 
    ## ?? try here add right DbPrefix: ? Uniprot/Uniref easy, others a mess.
    map { if(/:/) { } # asis
    	if(/^UniRef/i or /^[A-Z0-9]+_[A-Z][A-Z]+$/) { $_="TrEMBL:$_"; } 
    	elsif(/^ENS\w+\d\d+$/) { $_="ENSEMBL:$_"; }
      } ($refid,$repid);
    
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


sub parse_evgheader
{
  my($oid,$hdr,$trlen)= @_;

    ## this becomes param or not?
  my $pubid= $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?

  my $protid= $pubid;
  $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
  $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
  #  protein_id      gnl|CacaoGD|Thecc1EG016762p1

  my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      aaqual => "na", trlen => $trlen, cdsoff => "", cdsor => 1, 
      name => "", namepct => 0, dbxref => "na", cdd => "na" );  

  if( $genenames{$oid} ) {
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


sub getOkFileset
{
  my($okpath,$SUFFIX,$okfiles)= @_;
  my($oktr,$alttr)=("","");
  $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  my @okd=(); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  ($oktr) = grep /.okay\.($SUFFIX)/, @okd;  
  ($alttr)= grep /.okalt\.($SUFFIX)/, @okd; 
  return($oktr,$alttr,\@okd); ## change to (\@okd,$oktr,$alttr) : want @okd mostly
}

sub getFileset
{
  my($okpath,$SUFFIX,$okfiles)= @_;
  $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  my @okd=(); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  my @files = grep /\.($SUFFIX)/, @okd; my $nok=@okd;
  #ok here# warn "#getFileset($okpath,$SUFFIX,$okfiles)= dir:$nok, suf:@files \n" if($DEBUG); 
  return(\@okd, @files); 
}

sub tidyupFileset { 
	my($tod,@td)= @_;  
	return unless($tod and @td);
	my @tdlist;
	mkdir($tod) unless(-d $tod);
	foreach my $fn (@td) { 
  	if(-f $fn and not ($fn =~ m,$tod/,)) { 
  		(my $tfn=$fn)=~s,^[\w\.]+/,,;   ## ASSUMES tod is subdir in current path !!
  		# my $err= runcmd("mv $fn $tod/$tfn"); # system mv or perl 
  		rename($fn,"$tod/$tfn");
  		push @tdlist, "$tod/$tfn";
  	} 
  } 
	my $n=@tdlist; my $n1=_min1($n-1,4); loggit(0,"tidy: n=",$n, @tdlist[0..$n1]); 
	#loggit(0,"tidyup: n=",scalar(@tdlist), (grep/\w/,@tdlist[0..4])); 
}



sub getmRNA   ## move to cdnasubs.pm ?
{
  my($okpath,$trname,$pubdir,$ADDutrorf)= @_;
  my($cdnaseq)=(""); # == mrna (oriented), not cdna.parts.tr
  
	use constant ALSOMAKE_AACDS => 1;
  $ADDutrorf= 1; # what? always check for okayset/*.utrorf.mrna ?
  
  #? FIXME? suffix .tr may not be used: .cdna? .fa? .fasta? ...
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna|cds|aa,  
   ## ?? add .fa ?? BUT conflict in now publicset/*.{mrna,cds,aa}_pub.fa **

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
  # unless($trf) { ($trf)= grep /\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not?
  ## FIXME? add \.mrna_pub.fa[.gz]
  unless($trf) { ($trf)= grep /\.mrna_pub\.fa|\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not?
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
    
  ## FIXME: fail here if missing oktr ?
		$APPtraa2cds= findevigeneapp("prot/traa2cds.pl") unless($APPtraa2cds); #  
		my $tropts="-trout "; $tropts .= "-nomiss " if($ADDutrorf);

 		
    if($oktr and -f $oktr and -f $okaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    if($alttr and -f $alttr and -f $altaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    
    if($ADDutrorf and $okall > 0) {
 			my($okin) = grep /.utrorf.mrna$|.utrorf.mrna.gz$/, @okd;  
   		if($okin) { runcmd("cat $okin >> $cdnatmp"); loggit(0,"add $okin to $cdnatmp"); } # err check? loggit?
   		else { $ADDutrorf=0; } # dont do .aa,cds
    }
    $cdnaseq= $cdnatmp if(-s $cdnatmp);
    # loggit(1,"ERR: No-make mRNA $cdnatmp from $oktr,$okaa + $alttr,$altaa") unless($okall>1);
    loggit(LOG_DIE,"FATAL: make mRNA $cdnatmp: missing files oktr=$oktr,okaa=$okaa,alttr=$alttr,altaa=$altaa,
    	from (oktr,alttr)=getOkFileset(path=$okpath,suffix=$TRSUFFIX)") unless($okall>1);
    
    # FIXmaybe: also make pubdir/.aa,.cds along with .mrna ? see hassle in update_mrna_fileset 
    ## FIXME2: ok*.aa and utrorf.aa have dup entries, use utrorf if there.
    if(ALSOMAKE_AACDS and $pubdir) {
  		my($ok,$hin,$hout,$fout,$okin,$altin,$utrin,); my %ids=();
  		foreach my $suf (".aa",".cds") {
				$fout= makename($cdnatmp,$suf);  
				($okin) = grep /.okay$suf$|.okay$suf.gz$/, @okd;  
				($altin)= grep /.okalt$suf$|.okalt$suf.gz$/, @okd; 
				($utrin)= ($ADDutrorf) ? grep(/.utrorf$suf$|.utrorf$suf.gz$/, @okd) : (); 
				if($okin and $altin and not -f $fout) {
					$ok= open($hout,'>',$fout);	my(%did,$did);
					if($utrin) { ($ok,$hin)= openRead($utrin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin); }
					($ok,$hin)= openRead($okin);  
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin);
					($ok,$hin)= openRead($altin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin);
					close($hout);
					map{ $ids{$_}{$suf}=$did{$_} } keys %did;
# 					($ok,$hin)= openRead($okin);  while(<$hin>){ print $hout $_; } close($hin);
# 					($ok,$hin)= openRead($altin); while(<$hin>){ print $hout $_; } close($hin);
# 					if($utrin) { ($ok,$hin)= openRead($utrin); while(<$hin>){ print $hout $_; } close($hin); }
					}
				}
				
  	## should check for ID agreement among output .mrna, .aa, .cds  
			my $nidok=0;
			($ok,$hin)= openRead($cdnaseq); while(<$hin>){ if(/^>(\S+)/) { $ids{$1}{'.mrna'}++; } } close($hin);
			foreach my $id (sort keys %ids) {
				my @suf= sort keys %{$ids{$id}}; 
				if(@suf==3) { $nidok++; } else { loggit(1,"ERR:getmRNA-misspart $id:",@suf); }
			}	
			loggit(0,"getmRNA $pubdir/mrna,aa,cds nid=",$nidok);
    }
  }
  
  return($cdnaseq);  # , $aaseq, $cdsseq    
}


## from asmrna_trimvec : generalize mrna_update_fileset ..
## see tr2aacds.pl:asmdupclass_fileset
sub update_mrna_fileset
{
  my($trpath, $inmrna, $updatename, $trimids, @trimfiles)= @_; 
  my $upstatus=0;
  # outputs now should go to pubdir, but use inmrna path.
  my($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles; # hash instead? for below %fset
	## return hash of valid ids?
	my %okids=();
	
	my $ntrim=0;
	if($trimids and ref($trimids)) {
		$ntrim= scalar(keys %$trimids);

	} elsif($trimidfile and -f $trimidfile) { # hash of ids in trimfiles, regenerate from trimidfile
   	$trimids={};  
   	my($ok,$hin)= openRead($trimidfile); 
   	while(<$hin>) { if(/^\w+/) { 
   		my($id,$action,@info)=split; # parse more of trimidfile? expect: id, action=DROP/OKCUT/PROBLEM, trimnotes, aanewhdr 
   		$trimids->{$id}=$action; $ntrim++; } 
   	} close($hin);

	} else { # ERROR unless -f flagtrimvec
 		#below# loggit(1, "ERR update_fileset  empty trim ids= $trimidfile"); 
	}
	
## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??

## ........ **#@&@!*% this file name wrangling is a big waste of time ..........
## ........ use fewer naming choices ???  no okdir for pubdir set
  
	my $flagtrimvec= makename($inmrna,'.'.$updatename); 
  my $outaa = makename($inmrna,".aa");  
  my $outcds= makename($inmrna,".cds"); 

  my($pubdir);
  (my $mrnapath= $inmrna) =~ s,/[^/]+$,,; ## THIS MAY BE WRONG NOW .. publicset/ vs okayset/
  # (undef,undef,$pubdir)= getOkFileset($mrnapath);
  ($pubdir)= getFileset($mrnapath);
  
  ## drop okalt here? see above getmRNA/ALSOMAKE_AACDS, expect/require pubdir/aa,cds?
  my $aapatt= basename($outaa);
  my($inaaseq) = grep /$aapatt$|$aapatt.gz$/, (@$pubdir);  
  if($inaaseq) { $outaa=$inaaseq;  } else { } # fail?
  
  my $cdspatt=basename($outcds);
  my($incdsseq) = grep /$cdspatt$|$cdspatt.gz$/, (@$pubdir);  
  if($incdsseq) { $outcds=$incdsseq; } else { } # fail?
  
  ## also return \%okids ??
  # return ( 0, $inmrna, $outaa, $outcds, \%okids) if( -f $flagtrimvec); #check for flag-file made below  
  ## FIXME:for -novectrim set updatename == 
  if( $updatename =~ /SKIP/ or -f $flagtrimvec) {
  	my($ok,$hin)= openRead($inmrna); while(<$hin>) { if(/^>(\S+)/) { $okids{$1}=1; } } close($hin);
  	# what if outaa,outcds missing?
  	return ( 0, $inmrna, $outaa, $outcds, \%okids);
  }
  if($ntrim < 1 or ! -s $trimmrna) { # fixme for -novectrim
		loggit(1, "ERR update_fileset  empty trim files= $trimmrna, $trimidfile"); 
		return (-1, $inmrna, $outaa, $outcds, \%okids);
  	}

  my $upmrna  = makename($inmrna,".mrna_upd"); 
  my $upaaseq = makename($inmrna,".aa_upd"); # failed empty outaa from above
  my $upcdsseq= makename($inmrna,".cds_upd");  # failed empty outcds from above
  
  $outaa =~ s/.gz$//; $outcds =~ s/.gz$//; # output fupname
  (my $outmrna= $inmrna)=~ s/.gz$//;
  my %fset= (
            # $nup,$nsame,$fin,$ftrim,$fup,$fupname
    mrna => [ 0, 0, $inmrna, $trimmrna, $upmrna, $outmrna ],
    aa   => [ 0, 0, $inaaseq, $trimaa, $upaaseq, $outaa],  #?? fix here for .okay.aa + .okalt.aa ?
    cds  => [ 0, 0, $incdsseq, $trimcds, $upcdsseq, $outcds],
  );
  
  foreach my $suf (sort keys %fset) { # is inmrna always .mrna ?
    my($ok,$hin,$hup,%keptids);
    my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};
      
    if(-s $fin and -s $ftrim) {
      %keptids=();   
      $ok= open($hup,'>',$fup); # unless($ok)...
      ($ok,$hin)= openRead($fin);  
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { my $d=$1; $ok=($trimids->{$d})?0:1; if($ok){ $keptids{$d}++; $nsame++; } } 
      	print $hup $_ if($ok); }
      close($hin);
      
    	## pull trimset/uvcut.$suf, check? for trimids{id} and/or collect above kept ids
      ($ok,$hin)= openRead($ftrim); 
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { my $d=$1; $ok=($trimids->{$d} and not $keptids{$d})?1:0; 
      		if($ok){ $keptids{$d}++; $nup++; } }
      	print $hup $_ if($ok); } 
      close($hin);
      close($hup);
      $fset{$suf}->[0]= $nup; $fset{$suf}->[1]= $nsame;  # $fset{$suf}= [$nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname];
      $upstatus++ if($nup>0);
    	%okids= %keptids if($suf eq 'mrna');
    } else {
      ## error, warn/loggit, skip?
      loggit(1, "ERR update_fileset.$suf empty $fin or $ftrim"); 
    } 
  }
  ## end merge loop

  my (@outfiles,@tmpfiles);
  if($upstatus == 3)  { # $upstatus == 3 or > 0?
    ## rename input files to input.old, better: input.untrim .notrim? .pretrim? .old?
    ## rename newmerge files to input
    foreach my $suf (sort keys %fset) { 
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      if(-s $fupname) { rename($fupname,"$fupname.untrim"); push @tmpfiles, "$fupname.untrim"; }
      rename($fup,$fupname); push @outfiles, $fupname;
      ## push @tmpfiles, $ftrim; # from @trimfiles. dont call tmpfile? 
      loggit(0, "update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim > $fupname"); 
    }
    runcmd("touch $flagtrimvec"); ## touch flag-file that new input has uvtrim results ..
  } else {
    loggit(1, "ERR update_fileset missing status=$upstatus/3");  # list fupnames?
    foreach my $suf (sort keys %fset) { # sort: aa,cds,mrna
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      push @outfiles, $fup;
      loggit(1, "ERR update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim >$fup"); 
    }
  }
  
  ## also return \%okids
  return ($upstatus, \@outfiles, \@tmpfiles, \%okids); # $upaaseq, $upcdsseq,$upmrna, 
}


# sub revcomp {  # same as in cdna_proteins.pm
#   my ($seq) = @_;
#   my $reversed_seq = reverse ($seq);
#   $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
#   return ($reversed_seq);
# }

sub facount {
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa); # openRead()
  if($ok) { while(<$hin>) { $n++ if(/^>/); } close($hin); }
  return $n;  
}

sub fasize {  
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa);# openRead()
  while(<$hin>) { chomp; $n += length($_) unless(/^>/); } close($hin);
  return $n; 
}

sub fasizes_nogap {  
  my($fa,$okayc,$gapc)= @_;
  $okayc ||= 'ACGTagct'; $gapc ||= 'Nn';
  if($okayc =~ /^aa|amino|prot/) { $okayc='A-WYZa-wyz'; $gapc='Xx\*'; }
  my ($nokay,$ngap,$nt,$id)=(0,0,0,0); # NOT total, per record
  my %fasizes=();
  my($ok,$hin)= openRead($fa);
  while(<$hin>) { 
    if(/^>(\S+)/) {
      if($id) { $fasizes{$id}= join"\t",$nokay,$nt,$ngap; }
      $id=$1; ($nokay,$ngap,$nt)=(0,0,0);
    } else {
      chomp; s/\*$//;
      $nt += length($_);
      $nokay += tr/$okayc/$okayc/; # tr/A-WYZa-wyz/A-WYZa-wyz/; ## aa chars; na chars=/ACGTagct/ gaps=/Nn/
      $ngap  += tr/$gapc/$gapc/;  # tr/Xx\*/Xx\*/;
    }
  } close($hin);
  if($id) { $fasizes{$id}= join"\t",$nokay,$nt,$ngap; }
  return \%fasizes; # ($nokay,$nt,$ngap);
}


sub fasplit { #  this splits by size; need alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splsize, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= (); my ($ok,$hin)= (0,undef);
  my($atsize,$atfile)=(0,0);
  my $fabase= basename($fa);
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splsize"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  while(<$hin>) {
    if($atfile==0 or ($atsize > $splsize && /^>/)) {
      $atfile++; if($atfile>$npart) { } # finishup???
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname);  $atsize= 0;
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $atsize+= length($_) unless(/^>/);
  } 
  close(SPL); close($hin);
  return @splist;
}

sub fasplitcount { # this splits by size; need alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splcount, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= ();
  my($ok,$hin, $irec,$atfile)=(0) x 10; 
  my $fabase= basename($fa); #? is this bad
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splcount"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  
  while(<$hin>) {
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
  close(SPL); close($hin);
  return @splist;
}

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
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
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
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}

sub runcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"CMD=",@cmd);  
  ## fail if $cmd[0] =~ /MISSING_/
  my $err= ($dryrun) ? 0 : system(@cmd);  ## ?? add option run: @result=`$cmd`;
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

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

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';  ## fixme need insuf: tr|fasta|fa
  unless($infile) { warn "cdna_evigenesub:makename MISSING infile"; $infile="Noname$$"; } # or die / buggy soft
  my $outfile= $infile; $outfile =~ s/\.gz$//; # bad for infile empty/undef .. do what?
  $outfile =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}

=item evigene_config

	GetOptions(.. "config=s", \$configfile, "cadd=s", \@configadd,);
	$configh= evigene_config($configfile, \@configadd); # always even if $config null

	$DBID= $configh->{general}->{dbid} || "MyDBID"; 
	$LOCUSTAG= $configh->{general}->{locus_tag} || $DBID; 
	$IDPrefix= $configh->{pubopt}->{publicid} || "Evigene"; 
	
 	from  evigene/scripts/evigene2genbanktbl.pl
	but drop global special hashes: %evidence; %public_options; %geneset; %programs; 
	change to: return \%config # replace other config hashes w/ this 2-level: config{part}{key} = val

=cut

sub evigene_config {
  my($cfile, $addoptions)= @_;
  my %config=(); my $ctype=0;

  use constant{ kGENERAL => 'general', kEVFILE => 'evidence', kEVPROG => 'programs', kPUBOPT => 'pubopt',  };
  	#old: kEVOPT => 'evoption', kANOPT => 'anoption',  kEVGENES => 'geneset',
  	
  if($cfile) { #  and -f $cfile
    open(F,$cfile) or die "ERROR reading config: $cfile";
  
    my @CONFIG= <F>;
    push @CONFIG, @$addoptions if(ref $addoptions);
    
    my ($lastkey, $lastval);
    foreach (@CONFIG)
    {  # key => value
      s/^\s+//;  s/\#.*$//; s/\s+$//; # end of line comments 

    ## need now to handle continuation lines, end with \
    ## or prefix with "+" << use this

      my($key,$val);
      if(/^[+\.]/ and $lastkey) { $key= $lastkey; s/^.//; $val= $lastval.$_; }
      elsif($lastval =~ m,\\$, and $lastkey) { $key= $lastkey; $lastval=~s,\\$,,; $val= $lastval.$_; }
      else { ($key,$val)= split(/[=\s]+/, $_, 2); }
      
      next unless($key =~ /^\w/); 
      $val =~ s/\s*$//;  $val =~ s/\\n/\n/g; $val =~ s/\\t/\t/g; # dequote tabs,newlines
      # allow for val == regex in match/m or subs/s;  
      #  names:
      #   cutdbx = s/\((InterPro|TAIR):[\w\.-]+\)//g ?
      #   isgeneid = m/^(Os|At|AT)\d{1,2}g\d{3,}/

# revise to parse '^section:' into separate hash
      if($key =~ s/:$//) { $ctype=$key; }
      
      if($key =~ /^evidence/) { $ctype= kEVFILE; } # old/fixed set
      elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
      elsif($key =~ /^program/) { $ctype= kEVPROG; }
      elsif($key =~ /^end$/) { $ctype= 0; }
#       elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
#       elsif($key =~ /^anoption/) { $ctype= kANOPT; }
#       elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
#.. drop special config hashes...
#       elsif($ctype eq kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
#       elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
# #       elsif($ctype eq kEVOPT) { $evaluate_options{$key}= $val; } 
# #       elsif($ctype eq kANOPT ) { $annotate_options{$key}= $val; } 
# #      elsif($ctype eq kEVGENES) { $geneset{$key}= $val; }
#       elsif($ctype eq kPUBOPT ) { $public_options{$key}= $val; } 
#       elsif($ctype eq kEVPROG) { $programs{$key}= $val; }

      # generic keys: name date genome .. other?
      if($key =~ /\w/ and $val =~ /\w/) { 
        my $ogroup= $ctype || "general";
        $config{$ogroup}->{$key}= $val; # this one #which?  $config{$ogroup}{$key}= $val;
      }
      
      # also for now : overlap, other progs
      ($lastkey, $lastval)=($key, $val);
    } close(F);
  }
  return \%config;
}


1;
__END__
