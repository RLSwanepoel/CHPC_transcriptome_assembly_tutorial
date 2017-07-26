#!/usr/bin/perl
# splign2gff.pl : convert ncbi splign (mrna x genome) exon alignment table to gene.gff

use strict;
use warnings;
use Getopt::Long;

=item about

  convert ncbi splign to gff
  http://www.ncbi.nlm.nih.gov/sutils/splign/splign.cgi

=item FIXME

  -- need to filter out low qual 2ndary aligns,
    use splog table ** for each trid, order by -score, remove 2nd scores below 99%? 95%? of top
  ** splog score NOT GOOD for complex cases
    -- need also score canonical/standard splice bases; often have fwd-fwd good-splice, rev-rev poor-splice, higher score
    -- AND for split-genes, need align-span/scaffold checks: must keep all align parts that cover mRNA span
      
=cut


my $debug=0;
my $IDPREFIX="evgs"; 
my $SRC="splign";
my $MRNA_ATTR='clen|offs|oid|aalen|Name';
my $DUPMIN=  0.90; # 0.99; # score percent of best align : change this, splog score is unreliable
my($infile,$output,$inlog,$dryrun,$trinfo,$insorted)= (0) x 10;
$output= undef;

my @saveopt= grep /^\-/, @ARGV;
my $optok= GetOptions(
  "input|tblfile=s", \$infile,
  "logfile=s", \$inlog,
  "trinfo=s", \$trinfo, # mrna headers from mrna.fasta ?
  # "mrna|cdna=s", \$cdnaseq,
  # "class|trclass=s", \$trclass,
  "output|gff:s",  \$output,
  #unused# "idprefix=s", \$IDPREFIX,  
  "SRC|source=s", \$SRC,  
  # "dryrun|n!", \$dryrun, 
  "sorted!", \$insorted, 
  "debug!", \$debug, 
  );

die "convert ncbi splign mRNA alignments to gff gene locations.
usage: splign2gff.pl -in mymrna.splign  -log mymrna.splign.log [ -trinfo mymrna.headers ] [-out my.gff ]
  -debug -source xyz -sorted
  -trinfo is mRNA fasta headers with protein info as from Evigene: $MRNA_ATTR
" unless($optok);

my(%trinfo);
readTrinfo($trinfo) if($trinfo);

my($ntrkeep,$ntrdup,%trkeep,%trhasdup); $ntrkeep=0;
$ntrkeep= readSplignLog($inlog) if($inlog);

my $INH= *STDIN;  # require file for sort ?
# if($infile) { open($INH,$infile) or die "reading $infile"; }  # add gzip open
## sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n  input.splign  # << by trid / ipart / scaff / sbeg / send
if($insorted) {
	if($infile and $infile !~ /stdin|^-/) { open($INH,$infile) or die "reading $infile"; }
} else {
	die "ERR: need splign file as input to sort" unless($infile and -f $infile);
  open($INH,"sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n $infile |") or die "sort/read $infile";
}

my $OUTH= *STDOUT;
if(not $output and defined $output and $infile) {
  $output= $infile; $output =~ s/\.\w+$//; $output.=".gff";
}
if($output) { open($OUTH,'>',$output) or die "writing $output"; }  

my($ntrout,$ltid,$ltidi,$lidir,$gaps)= (0) x 10; 
my(@xon,%trdup,%trprintdup,%trexons);
my $lgeneid=""; my @geneset=();


## FIXME here? sort input.splign by trID? get all idir/trid, but preserve exon order per idir?
## sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n  input.splign  # << by trid / ipart / scaff / sbeg / send
#
# +875	FunhetEGm056896t1	Scaffold348	0.933	360	604	963	79954	80313	GG<exon>  	M16RM18RM28RM43RM16RM79RMRM12RM16RM2RM3R2MRM4RM2RM2RM20RM8RM5RM2RM5RM2RM12RM22RM17
#	# .. rev map is slightly different: 77082 .. 77306,77313 
#	-875	FunhetEGm056896t1	Scaffold348	0.932	353	963	611	80313	79961	  <exon>GT	M17RM22RM12RM2RM5RM2RM5RM8RM20RM2RM2RM4RMR2M3RM2RM16RM12RMRM79RM16RM43RM28RM18RM9

=item FIXME: need more parts scoring

	splign.log score and trkeep not good enougn
 ** splog score NOT GOOD for complex cases
    -- need also score canonical/standard splice bases; often have fwd-fwd good-splice, rev-rev poor-splice, higher score
    -- AND for split-genes, need align-span/scaffold checks: must keep all align parts that cover mRNA span

  ** need 2-pass thru splign table, per trancript (are all parts of trID together?)
    -- accumulate all of trID idir parts: +1, -1, +2, -2 ... then decide which to keep
    -- count canon splices, fwd: AG<exon>GT  rev: AC<exon>CT
    -- tabulate trid:span per idir, 
    
  ** But there are those 1-exon things that map 100s of places, slight diff in align score, no splice score.
  
=cut
		
sub MAINhere {}

while(<$INH>) {
	next unless(/^[+-]\d/); #? all ok?
	chomp;
	my($idir, $tid, $gref, $pident, $alen, $tb, $te, $gb, $ge, $xsplice, $mismats) = split"\t";

## change this trkeep/skip: dont skip from splog score of trkeep just yet; calc gene score from exons, canon splices, ..
## then decide to keep or drop.  For fwd/rev at same locus, need to look at both first here to decide?
## .. OK now, trkeep has up to MAXTOP=10 dups, regardless of log score, but this filter removes some trash w/ 100s align
  my $keepit= ( $ntrkeep > 0 ) ? $trkeep{$tid}{$idir} || 0 : 1;
  next unless($keepit);
    
	my $tidi= $tid.$idir;  # Trid+100  Trid-100 ..
	## DEFER  putgff() : put into trid array, then decide..
	if($tidi ne $ltidi) { 
		if(@xon) {
			my($geneid,$genescore,$genomespan,$mrnaspan,$genearr)= 
				putgff($ltid,$lidir,$gaps,\@xon);
			# returns 0 for skipit...
			if($geneid) {		
			if($geneid ne $lgeneid) {
				processGeneSet($lgeneid, @geneset ) if($lgeneid && @geneset); @geneset=();
				}	
			push @geneset, [$geneid,$genescore,$genomespan,$mrnaspan,$genearr]; #??
			$lgeneid= $geneid;
			}
		}
			
		$gaps=""; @xon=(); 
		}

  ## fixme: orient calc; NOTE tb,te,gb,ge are '-' for Gaps, other non-exon,non-locations
	# idir: +/- orient, item number
  my($dor)= $idir =~ /^([+-])/; # arghhhh .. - means mRNA is reversed (and tb > te ??)
	## ^^ irrelevant for needs here, use tb>te and gb>ge
 	my($strand,$tor,$gor,$hasloc)= ('.',0,0,0);
	$tor=1; if($tb > $te) { $tor=-1; ($tb,$te)= ($te,$tb); } # dont need rev trspan?
  ## FIXME: need to keep gor,tor and pass on to putgff splice test, or do splice flip here?
  if($gb =~ /\d/ and $ge =~ /\d/) {
  	$hasloc=1;
		$gor=1; if($gb > $ge) { $gor=-1; ($gb,$ge)= ($ge,$gb); } # should match idir; orient?
		$strand= ($gor * $tor < 0) ? '-' : '+';
	}
	
	# xsplice:  lbp <exon|gap|..> rbp; gap: M-Gap, L-Gap, what else?
	## if gap, then missing vals: pident, gb, ge, mismats
	my $xtype="dunno"; my $splchar="nnnn"; 
	my ($splscore,$rsplscore)= (0,0);
	if($xsplice =~ /exon/) {
  	$xtype="exon"; 
  	# my($sl,$sr)= $xsplice =~ m/(..).exon.(..)]/;  # AG<exon>GT .. 
  	my $sl = ($xsplice=~m/(..).exon/)?$1:'nn';
  	my $sr = ($xsplice=~m/exon.(..)/)?$1:'nn';
  	map{ $_='nn' if($_ eq '  ') }($sl,$sr);
  	
  	## ?? do splice flip here if gor<0 and +splice, or gor>0 and -splice?
  	## maybe not here..
#   	if($gor < 0 and $tor < 0) {
#   		$sl='AC' if($sl eq 'AG');
#   		$sr='CT' if($sr eq 'GT');
#   	} elsif($gor > 0 and $tor > 0) {
#    		$sl='AG' if($sl eq 'AC');
#   		$sr='GT' if($sr eq 'CT');
#   	}
  	
  	my $splor=($gor < 0)?'-': ($gor > 0)?'+' : ''; ## buggers gor not enough? gor*tor ?
  	# my $splor=$strand; ## buggers gor not enough? gor*tor ?
		$splchar=$sl.$sr.$splor;
		##$splchar.='-' if($gor < 0); # DOES THIS DO IT RIGHT?
		
		## this gor> gor< is probably wrong need tor also? ie splice eq '-' or '+'
		## eg. FunhetEGm000081t3, -strand, 20 exons, all 'AGGT' splices .. lots like this, and reverse:
		##   +strand, spice=ACCT
		## sum: splice=4,AGGT = 10258, splice=0,AGGT'  = 8933
		##     splice=4,ACCT = 56, 'splice=0,ACCT' = 29 << this is not good.
		# Scaffold200     splign  mRNA    215515  223375  100     -       .       ID=FunhetEGm000081t3;
		# cov=100%,3976/3991;nexon=20;splice=0;Target=FunhetEGm000081t3 1 3991;gescore=20;
		# splinfo=+14924,nd0,sc176.854;aalen=1329,99%,partial;clen=3991;offs=3-3989;oid=kfish2eg6velvk55Loc85687t1
		## gmap has Scaffold200:215515-223375:-,cov=100% == Titin.
		
		# right?  MOVE this to putgff scoring all exons, keep splice=splchar
		if($sl eq 'AG') { $splscore+=2; } elsif($sl eq 'AC') { $rsplscore+=2; } 
		if($sr eq 'GT') { $splscore+=2; } elsif($sr eq 'CT') { $rsplscore+=2; } 
		#wrong# if($gor > 0) { $splscore+=2 if($sl eq 'AG'); $splscore+=2 if($sr eq 'GT'); }
		#wrong# elsif($gor < 0) { $splscore+=2 if($sl eq 'AC'); $splscore+=2 if($sr eq 'CT'); }
	} else { 
		($xtype=$xsplice) =~ s/[<>-]//g; $gaps.="$xtype:$tb-$te,"; 
	}	
	
	
=item exon rows

	record 'AG<exon>GT' splice bases around exon for other uses,
	canon splice: fwd: AG<exon>GT  rev: AC<exon>CT
	spaces '  ' for no spl char/gap
	
   my($spl,$spr)= $xsplice =~ m/(..).exon.(..)]/;  
		what are variants on pattern?  ??<exon>?? missing/blank chars?  
     'CC<exon>  ' : blank splice before/after Gap
   need to revcomp xsplices, if or>0, but tor < 0?
   ALSO, can use canonical/std splice bases as qual score; some of splog scores are higher for rev-rev case,
     but those had many fewer std splices .. shifted exon breaks?

=item xsplice frequency

  according to augustus, GT/AG or GC/AG are two valid splice patts. AG-GT,AG-GC here
  547720 AG<exon>GT     canon. splice
   88468 AC<exon>CT      rev canon : when does this occur?
        ^ is this valid when mrna-rev + genome-rev ?
          splign is making some backass map orients, is this the case?
          
  350265 <L-Gap>        gaps
  346391 <M-Gap>
  333549 <R-Gap>
  42532   <exon>  
  
  62924   <exon>GT      half-splices
  61712 AG<exon>  
  22567 AC<exon>  
  17501   <exon>CT
  
  44958 AC<exon>CC      which of these are valid alt spicing?
  26913 TA<exon>CT
  23511 AC<exon>AC
  20438 AG<exon>GC     << count as valid? but occur often w/ many mismatches
  19190 CA<exon>CT
  18062 TA<exon>CC
  17852 AG<exon>AA
  17385 AG<exon>AT
  16846 AG<exon>CT
  16775 AC<exon>GT

=item AC-CT cases

>> this is useful case: various mappers have diff exon set, drop-outs as gap where good exon expression shown.
   .. introns give alternate paths that are not captured by tr-mapping.
   >> take a look at all of cufflinks direct location alt-transcripts, not tr-gmap. does it give same exon drop-outs?
       kf2bcuf2_Gsc74g70000t1, or  fungrk2cuf2_Gsc413g51437t1,
       kf2bcuf2_Gsc74g70000t1 kf2bcuf2_Gsc74g70000t2 kf2bcuf2_Gsc74g70000t3 kf2bcuf2_Gsc74g70000t4
       ^^ these all have same mappings as velv,soap,trin trasm, missing some of valid alt exon/intron cases
       .. splign gets other exon mappings of same tr than gmap or cufl.
       .. cant tell if this is trasm or genomeasm mixup, both?
       well-known gene: human:UniRef50_O75165, Name DnaJ subfamily C member 13; 2200aa : size requires extra exons

-323    FunhetEGm027001t1       Scaffold74      1       185     4619    4435    539792  539608  AC<exon>CT      M185
-323    FunhetEGm027001t1       Scaffold74      1       116     4219    4104    539076  538961  AC<exon>CT      M116
-323    FunhetEGm027001t1       Scaffold74      1       105     3257    3153    533407  533303  AC<exon>CT      M105
-323    FunhetEGm027001t1       Scaffold74      1       102     3152    3051    532137  532036  AC<exon>CT      M102
-323    FunhetEGm027001t1       Scaffold74      1       63      3050    2988    531940  531878  AC<exon>CT      M63
-323    FunhetEGm027001t1       Scaffold74      1       160     2987    2828    531768  531609  AC<exon>CT      M160
-323    FunhetEGm027001t1       Scaffold74      1       83      2569    2487    530868  530786  AC<exon>CT      M83
-323    FunhetEGm027001t1       Scaffold74      1       18      2486    2469    528828  528811  AC<exon>CT      M18
-323    FunhetEGm027001t1       Scaffold74      1       201     710     510     518887  518687  AC<exon>CT      M201
-323    FunhetEGm027001t1       Scaffold74      1       42      509     468     517005  516964  AC<exon>CT      M42
-323    FunhetEGm027001t1       Scaffold74      1       150     467     318     516648  516499  AC<exon>CT      M150
-323    FunhetEGm027001t1       Scaffold74      1       76      317     242     515179  515104  AC<exon>CT      M76

-7557   FunhetEGm021830t1       Scaffold10163   1       106     1592    1487    603990  603885  AC<exon>CT      M106
-7557   FunhetEGm021830t1       Scaffold10163   1       96      1486    1391    602760  602665  AC<exon>CT      M96
-7557   FunhetEGm021830t1       Scaffold10163   1       202     1390    1189    602558  602357  AC<exon>CT      M202
-7557   FunhetEGm021830t1       Scaffold10163   1       128     1188    1061    602156  602029  AC<exon>CT      M128
-7557   FunhetEGm021830t1       Scaffold10163   1       125     1060    936     601905  601781  AC<exon>CT      M125
-7557   FunhetEGm021830t1       Scaffold10163   1       112     935     824     601660  601549  AC<exon>CT      M112
-7557   FunhetEGm021830t1       Scaffold10163   1       72      584     513     598861  598790  AC<exon>CT      M72
-7557   FunhetEGm021830t1       Scaffold10163   1       29      512     484     598669  598641  AC<exon>CT      M29

>> utrbad : gene join?
-7704   FunhetEGm027172t1       Scaffold2       1       132     3124    2993    2562226 2562357 AC<exon>CT      M132
-7704   FunhetEGm027172t1       Scaffold2       1       52      2992    2941    2562482 2562533 AC<exon>CT      M52
-7704   FunhetEGm027172t1       Scaffold2       1       151     2712    2562    2564094 2564244 AC<exon>CT      M151
-7704   FunhetEGm027172t1       Scaffold2       1       71      2561    2491    2564343 2564413 AC<exon>CT      M71
-7704   FunhetEGm027172t1       Scaffold2       1       194     2328    2135    2579844 2580037 AC<exon>CT      M194
-7704   FunhetEGm027172t1       Scaffold2       1       55      2134    2080    2580400 2580454 AC<exon>CT      M55
	
=item non-exon rows

4248 <L-Gap>  == left-end (5p?) gap
4476 <M-Gap>  == middle-gap ?
3795 <R-Gap>  == right-gap 
 284 <poly-A> == what?
 165 <poly-T> == what?

=cut

      ## FIXME PROBLEM : skip dups if they cover same/nearly genomic locs, in reverse... need best prot strand
      ## save exons? check $trkeep{$tid}{$d} dup count?
      ## PROBLEM 2 	: split-genes, cant use single top score, need to look at align spans,
      ##   when align parts are disjunct on sep Scaffolds, must keep all.      
        
	if($hasloc) {  # only exons have location
 		# insert comment or what for no-location rows? 
   
    ## count R, D, I, other? r=replace? d=del, i=ins;  pident maybe == mmsum/alen
    my $mmsum=0; while($mismats =~ m/M(\d*)/g) { my $m=$1; $m=1 unless($m); $mmsum+=$m; }
    my $indel=0; while($mismats =~ m/[ID](\d*)/g) { my $m=$1; $m=1 unless($m); $indel+=$m; }

# 		if($rsplscore > $splscore) { 
# 			# want flipped case? need genesum of rsplscore > splscore
# 			# splice=-$rsplscore ? or splice=$splscore,$rsplscore,$splchar ?
# 		}

    #o# my $attr="Parent=$tid;Target=$tid $tb $te;splice=$splscore,$splchar;align=$mmsum/$alen";
    my $attr="Parent=$tid;Target=$tid $tb $te;splice=$splchar;align=$mmsum/$alen";
    $attr .=";indel=$indel" if($indel);
    #NOT here# $attr .=";gor=$gor" if($gor < 0); #?? need only for putgff splice fixup?
    push @xon, [$gref,$SRC,$xtype,$gb,$ge,$pident,$strand,".",$attr];
	}
	
	$ltid=$tid; $lidir=$idir; $ltidi= $tidi;
}

if(@xon) {
	my($geneid,$genescore,$genomespan,$mrnaspan,$genearr)= 
		putgff($ltid,$lidir,$gaps,\@xon);
	if($geneid) {
		push @geneset, [$geneid,$genescore,$genomespan,$mrnaspan,$genearr]; #??
		$lgeneid= $geneid;
	}
}
processGeneSet($lgeneid, @geneset ) if($lgeneid && @geneset); @geneset=();

close($OUTH); close($INH);
warn "# splign.gff: ntrout=$ntrout, \n" if $debug;

#-------------------------------

sub dupLoc
{
	my($trid,$exons)= @_;
	
	if($trdup{$trid}++) {
      ## PROBLEM : skip dups if they cover same/nearly genomic locs, in reverse... need best prot strand
      ## this trexons test for overlap is bad; use mrna span?
	  my $hasloc= 0;
	  foreach my $x (@$exons) {
	  	my $xbe=  join "", (split"\t",$x)[0,3,4];
	    # my $xbe= $x->[0] . $x->[3] . $x->[4]; # exact match or overlap test?
	    $hasloc++ if($trexons{$trid}{$xbe});	    
	  }
	  if($hasloc) { $trdup{$trid}--; return $hasloc; } #?? is return mistake? need score of dup before deciding
	}
	return 0;
}

sub  processGeneSet
{
	my( $geneid, @geneset )= @_;
	# geneset row == [$geneid,$genescore,$genomespan,$mrnaspan,$genearr];

	use constant kOVERSLOP => 19; # what? need to check, splign does some align adjusting to spread true tr-align
	use constant DROPDUPGENES => 0;
	
	if(@geneset == 1) {
		printgene(@{$geneset[0]});
	
	} elsif(@geneset>1) {
		@geneset= sort { $b->[1] <=> $a->[1] or $a->[3] cmp $b->[3] } @geneset; # sort gscore,mspan

if(DROPDUPGENES) { # this filt should be option, parts check keeps out dups...
		my @gsnodup=(); 
		for my $gn (@geneset) {  # remove dupl calls at same locus, usual is +dir,-dir for splign
			my @exons= grep /\texon/, @{$gn->[4]}; 
			my $hasdup= (@exons) ? dupLoc($geneid,\@exons) : 0;
			push @gsnodup, $gn unless($hasdup);
		}
		@geneset= @gsnodup;
}
			
		#wait: printgene(@{$geneset[0]}); #  look for more @parts before print; update attr if found; qnd answer
		
		my @parts=();  ## should re-sort geneset by mrnaspan. so dont get more overlaps.
		## need to check thru sgene for mrnaspan split parts
		my @gff1= @{$geneset[0]->[4]}; # too messy a structure.. drop all but @gff genearr?
		my($cov)= $gff1[0] =~ m/cov=(\d+)/; # or align= ??
		if($cov>1 and $cov<95) { # look for more parts?
			my($md1,$mb1,$me1)= split/[:-]/, $geneset[0]->[3]; # or gff1 Target=xxx b e
			foreach my $i (1 .. $#geneset) {
				my($md,$mb,$me)= split/[:-]/, $geneset[$i]->[3];
				if($mb > $me1 - kOVERSLOP) { push @parts, $geneset[$i];  $me1=$me; }
				elsif($me < $mb1 + kOVERSLOP) { push @parts, $geneset[$i];  $mb1=$mb; }
# 				if($mb > $me1 - kOVERSLOP or $me < $mb1 + kOVERSLOP) {
# 					push @parts, $geneset[$i]; # look for how many?
# 					if($mb > $me1 - kOVERSLOP) { $me1=$me; } # join.. # was: if($mb >$me1-9  and $mb < $me1+90)
# 					elsif($me < $mb1 + kOVERSLOP) { $mb1=$mb; }
# 				}	
			}
			#wait# foreach my $p (@parts) { printgene(@$p); } # and call it a chimera ! id=geneid_C2..
		}

		if(@parts) {
			unshift @parts, $geneset[0];
			my @order=();
			my $npart=@parts;
			
			foreach my $i (0..$#parts) { 	
				my($md,$mb,$me)= split/[:-]/, $parts[$i]->[3]; # mrna span
				push @order, [$i,$mb];
			}
			@order= sort { $a->[1] <=> $b->[1] } @order;
			
			# my @spart=(); foreach my $o (@order) { push @spart, $parts[$o->[0]]; }
			# foreach my $i (0..$#spart) 
			my $pn= 0;
			foreach my $o (@order) { 
				my $ip= $o->[0];
				my $part= $parts[$ip];
				$pn++; #my $pn=$i+1;
				my $geneid= $part->[0]; # $spart[$i]->[0];
				my $newid= $geneid."_C$pn"; 
				## ugh: FunhetEGm000040t3_C0, FunhetEGm000040t3_G2_C1
				## have to replace also $spart[$i]->[0]
				$part->[0]= $newid; # $spart[$i]->[0]= $newid;
				foreach (@{$part->[4]}) {  # @{$spart[$i]->[4]}
					 s/(ID|Parent)=$geneid/$1=$newid/; 
					 if(/\tmRNA/) { s,;,;Split=$pn/$npart;,; } # oops, pn is wrong after reorder
				}
				printgene(@$part); # @{$spart[$i]}
			}
			# foreach my $p (@spart) { printgene(@$p); } # and call it a chimera ! id=geneid_C2..
		
		} else {
			printgene(@{$geneset[0]}); # wait, look for more @parts before print; update attr if found; qnd answer
		}
			
		if($debug) {
			foreach my $i (1 .. $#geneset) {
				my($did,$dgenescore,$dgenomespan,$dmrnaspan,$dgenearr)= @{$geneset[$i]};
				print $OUTH "#d$i.",$dgenearr->[0] if($dgenearr); # show mRNA 2..n
			}
		}
	} 
	# else { printgene(@{$geneset[0]}); }
}

# $trinfo= readTrinfo($trinfo) if($trinfo);
# >FunhetEGm033615t1 type=mRNA; aalen=645,23%,complete-utrbad; clen=8340; offs=11-1948; oid=kfish2qf7soapk21loc3t1; organism=Fundulus heteroclitus;
# >FunhetEGm033619t1 type=mRNA; aalen=1366,37%,complete-utrpoor; clen=10794; offs=3242-7342; oid=kfish2qf7soapk21loc3t3; organism=Fundulus heteroclitus;

sub readTrinfo {
	my($trinfo)= @_;
  return {} unless($trinfo);
	## my $MRNA_ATTR='clen|offs|oid|aalen|Name';
  ## look for aalen=  offs=  clen= 
  open(my $INH,$trinfo) or die "reading $trinfo"; 
  while(<$INH>) {
    if(/^>(\S+)/) { my $tid=$1;
   		my $at= join ";", m/((?:$MRNA_ATTR)=[^;\s]+)/g;
   		$trinfo{$tid}= $at || "";
    }
  } close($INH);
	return \%trinfo; #?? global now
}

=item readSplignLog

# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
# +1	FunhetEGm046329t1	Scaffold0	Ok	-8.248		: has only 1 exon partly aligned, thus low -score
# +2	FunhetEGm046569t1	Scaffold0	Ok	-17.252
# ..
# +8	FunhetEGm057611t1	Scaffold0	Ok	-333.571
# +9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: + = forward
# -9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: - = reverse
#
## problem ?? from combining split-compute set: duplicate idir numbers, wont uniqely say which items are linked
## no, each trid is limited to one split-part.  

=cut

sub readSplignLog {
	my($inlog)= @_;
	%trkeep= %trhasdup= ();  # global now
  return 0 unless($inlog);
  my(%tscore,@dkeep);
  open(my $LOGH,$inlog) or die "reading $inlog"; 
  while(<$LOGH>) {
    next unless(/^[+-]\d/);
    chomp; my($idir, $tid, $gref, $okay, $score) = split"\t";
    $score = -$score; # largest -score is best.. NOT Always; have +/- often, some hiscore - are wrong dir
    $tscore{$tid}{$idir}= $score;
  } close($LOGH);
  
  my $MAXTOP=10;
  my @tid= sort keys %tscore;
  foreach my $tid (@tid) {
  	## FIXME: score from splign.log may not be best to use alone..
  	## 	also count canonical splices for same tr align fwd/rev at same locus..
  	## change to keep top 10? scores, so split-tr parts are kept (? 10 enough, only care for max 3,4 split parts)
  	
    my ($dtop,@d) = sort{ $tscore{$tid}{$b} <=> $tscore{$tid}{$a} } keys %{$tscore{$tid}};
    my $scoretop= $tscore{$tid}{$dtop};
    my $scut= $DUPMIN * $scoretop;
    $trkeep{$tid}{$dtop}=$scoretop; $trhasdup{$tid}=0;
    # push @dkeep, $dtop;
    my $ntop=1;
    foreach my $d (@d) {
      my $s=$tscore{$tid}{$d};
      if($s >= $scut or $ntop < $MAXTOP) { 
        $trkeep{$tid}{$d}=$s; # push @dkeep, $d;
        $trhasdup{$tid}++; $ntrdup++; $ntop++;
      } 
      ## PROBLEM : skip dups if they cover same/nearly exons, in reverse... need best prot strand
      } 
    $ntrkeep++;
    # $trkeep{$tid}= \@dkeep; #??
  }
  
  warn "# splign.log: ntr=$ntrkeep, ntrdup=$ntrdup\n" if $debug;
  return $ntrkeep;
}

sub printgene {
	my($trid,$genescore,$genomespan,$mrnaspan,$refgene)= @_;

	  ## FIXME: here: defered ID change til here.. from putgff; need to clear %trdup / use new hash
	my $id= $trid; my $changeid=0;
	## this also will be used for split-gene cases, where want same 1 ID but add attribs: Split=1/2,2/2 ..
	if($trprintdup{$trid}++) {  
	  my $idup= $trprintdup{$trid}; $changeid=1;
	  $id= $trid ."_G". $idup; # what ID format for 2ndary matches?  _G2..9 ?
	}
	
	foreach my $gff (@$refgene) {
		$gff =~ s/(ID|Parent)=$trid/$1=$id/ if($changeid);
		print $OUTH $gff; # has \n newline for now
	}
}

sub putgff {
	my($trid,$idir,$gaps,$exons)= @_;

	my $id= $trid;
	my @xon= sort{ $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] } @$exons;
  my $trattr= $trinfo{$trid} || ""; # source mrna attributes
  # clen= offs= oid= aalen= .. Name= ??

	## dup align check:
	## FIXME: need genescore of all such dups first, to decide which to keep.
	## NOT HERE, see sub dupLoc()	
# 	if($trdup{$trid}++) {
#       ## PROBLEM : skip dups if they cover same/nearly genomic locs, in reverse... need best prot strand
#       ## this trexons test for overlap is bad; use mrna span?
# 	  my $hasloc= 0;
# 	  foreach my $x (@xon) {
# 	    my $xbe= $x->[0] . $x->[3] . $x->[4]; # exact match or overlap test?
# 	    $hasloc++ if($trexons{$trid}{$xbe});	    
# 	  }
# 	  if($hasloc) { $trdup{$trid}--; return 0; } #?? is return mistake? need score of dup before deciding
# 	  
# 	  if(0) { # NOT NOW
# 	  my $idup= $trdup{$trid};
# 	  $id= $trid ."_G". $idup; # what ID format for 2ndary matches?  _G2..9 ?
# 	  }
# 	  ## FIXME:^ defer ID change now, till processGeneSet / printgene
# 	}
	
	my $tscore= $trkeep{$trid}{$idir} || 0;  $tscore =~ s/^\+//;
	my $tsdup=  $trhasdup{$trid}||0; #debug ??
	
	# make mRNA from exons, print
	## FIXME: mrna length needs gap spans added..
	my($mref,$mgb,$mge,$mor,$maln,$mlen,$mrnab,$mrnae,
		 $mnxon,$msplscore,$rsplscore,$gapspan,$ngaps)=(0) x 20;
		 
	if($gaps) {
	  my @gs= split",",$gaps; $ngaps=@gs;
	  foreach (@gs) { my($gt,$gb,$ge)= split/[:-]/,$_; $gapspan += 1+$ge-$gb; } # does this include polyA,polyT ?	  
	  $mlen += $gapspan;
	}
	
	# recalc cov= align + mismatchs, not indels or gaps ?
	$mnxon= @xon; my $gor=0;
	foreach my $x (@xon) {
		my($gref,$src,$xtype,$gb,$ge,$pident,$strand,$phx,$attr)= @$x;
		unless($mref) { $mref=$gref; $mor=$strand; }
		$mgb=$gb if($mgb==0 or $gb<$mgb); $mge=$ge if($ge>$mge);
		
		#o# if( my($splscore)= $attr=~m/splice=(\d+)/ ) { $msplscore+=$splscore; }  # score is # canon-splice-bases both sides of exon: 4,2,0.  div by 4? adjust ends? 2-max
		## add rev-splice score sum? splice=9,-12,... add n-gap count to score? or just gapspan?
		## right?  MOVE this to putgff scoring all exons, keep splice=splchar
		## use rsplscore > msplscore to flip strand, flag as sense=-1 as per gmap.gff
		# if($attr=~m/gor=-1/) { $gor--; } else { $gor++; } #splice= fix or add to splice= flag?
		if( my($sl,$sr)= $attr=~m/splice=(..)(..)/ ) {
			my $gorx= ( $attr=~m/splice=$sl$sr([+-])/ ) ? $1 : 0;
			if($gorx eq '-') { $gor--; } elsif($gorx eq '+') { $gor++; }
			if($sl eq 'AG') { $msplscore+=2; } elsif($sl eq 'AC') { $rsplscore+=2; } 
			if($sr eq 'GT') { $msplscore+=2; } elsif($sr eq 'CT') { $rsplscore+=2; } 
		}
		
		my($tgd,$tgb,$tge)= $attr=~m/Target=(\S+) (\d+) (\d+)/; # my $len= 1+$tge-$tgb;
		$mrnab=$tgb if($mrnab==0 or $tgb<$mrnab);
		$mrnae=$tge if($mrnae==0 or $tge>$mrnae);
		
		my($aln,$len)= $attr=~m/align=(\d+).(\d+)/;
	  $maln+=$aln; $mlen+=$len;
		#?? maybe drop align=100/100 from exon outputs, only mismatches?
		
	  my $xbe= $gref.$gb.$ge;# my $xbe= $x->[0] . $x->[3] . $x->[4]; # exact match or overlap test?
    $trexons{$trid}{$xbe}++; # for trdup test, only need if($trhasdup{$trid})
	}
	
	## FIXME3: is mis-orient problem in part due to bad ORF-call, mRNA is not right?  check for utrbad+misorient
	
	## FIXME2: need to flipor for both ways, test is strand vs fwdspl,revspl:
	##	main: $strand= ($gor * $tor < 0) ? '-' : '+';
	##  need genome-or only here, if fwdspl, expect gor > 0, flip if gor<0; vs revspl && gor < 0
	##  flipor=1 if(($strand eq '-' and fwdspl) or ($strand eq '+' and revspl))
	
  ## FIXME4: if both +score and -score are high, report in addat confusion == rev gene join likely

	my $flipor=0; my $addat="";
  my $msplscoresave= $msplscore;
	if($msplscoresave > 10 and $rsplscore > 10) {
	  $addat.=";splicemix=$msplscoresave,$rsplscore"; #? or append to splice=$msplscore, ?
	}
	if($rsplscore > 2+$msplscore and $rsplscore > 2) {
		$msplscore= $rsplscore;
		# if($gor > 0 and $mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; } ## and $tor > 0 ??
		if( $gor > 0 ) { 
		  if($mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; }
		  elsif($mor eq '-') { $addat.=";sense=-1";  }
		} elsif($gor < 0 and $mor eq '+') {
		  $mor='-'; $flipor=1; $addat.=";sense=-1"; 
		}
	} elsif( $msplscore > 2 ) { # messier, sense=-1 when mrna-flipped (tor <0); do flip both ways
		## if($gor < 0 and $mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; } ## and $tor < 0 ??
		if( $gor < 0 ) {
		  if($mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; } 
		} elsif($gor > 0 and $mor eq '-') {
		  $mor='+'; $flipor=1; $addat.=";sense=-1"; 
		}
	}
	
	if($trattr =~ /clen=(\d+)/) { $mlen=$1;	}
	my $pctalign= int(0.5 + 100 * $maln/$mlen);
	
## change to match gmap.gff attrib??
## aaln=479;cov=17.8;indels=43/0;nexon=10;pid=94.5;qlen=8340;
## qlen == clen; pid*coverage == malign/mlen ;  indels?  nexon?
## change exon score to percent from .proportion?
	
	##? genescore: multiply proportions for all score parts; max = 1; BUT not .90 aln x 0 splice .. need fudge
	my $genescore = $pctalign/100; 
	
	## fix for no splice at tr ends: +4 splice score: $msplscore+=4 if($msplscore>0)
	## fix for 1-exon spct= (2 + 4)/(5*2) = .60 ?? ;  (3 + 4)/(5*3) = 0.47
	
	use constant ONEX_SPLSCORE => 0.70; ## .50;
	if($mnxon>1) { my $spct= ($mnxon+$msplscore+4)/(5*$mnxon); $genescore *= $spct; } 
	else { $genescore *= ONEX_SPLSCORE; } # single-exon adjust? not good enough for dup choice w/ other has 2+

=item genescore problems
	
	.. problem deciding what part to weight more : align cover is most important, but
		 more valid splices w/ slightly lower cover is better (??)
		 
## bad genescore, splice=0 above splice=8 .. nxon=1 for bad choice, nxon=3 for good miss.
# Scaffold1810    splkf2eg7m1     mRNA    22562   23156   88      -       .       ID=Funhe2Eq7m000962t1;cov=88%,566/642;nex
#   on=1;splice=0;Target=Funhe2Eq7m000962t1 49 642;gaps=48,LGap:1-48,;gescore=79;splinfo=+8643,nd9,sc44.929;aalen=164,76%,par
#   tial3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1
# #d1.Scaffold474 splkf2eg7m1     mRNA    170935  187263  89      +       .       ID=Funhe2Eq7m000962t1;cov=89%,569/642;nex
#   on=3;splice=8;Target=Funhe2Eq7m000962t1 2 642;gaps=1,LGap:1-1,;gescore=65;splinfo=+1647,nd9,sc54.389;aalen=164,76%,partia
#   l3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1
# #d2.Scaffold271 splkf2eg7m1     mRNA    267767  268201  62      -       .       ID=Funhe2Eq7m000962t1;cov=62%,395/642;nex
#   on=1;splice=2;Target=Funhe2Eq7m000962t1 207 642;gaps=206,LGap:1-206,;gescore=56;splinfo=+7228,nd9,sc123.48;aalen=164,76%,
#   partial3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1

## more genescore quandry : is spliced low cov better/worse than unspliced hicov?
## >> all spliced, but 33% unaligned
# Scaffold264     splkf2eg7m1     mRNA    363123  388491  67      -       .       ID=Funhe2Eq7m000197t1;cov=67%,369/547;nex
#  on=4;splice=12;sense=-1;Target=Funhe2Eq7m000197t1 93 545;gaps=147,LGap:546-547,MGap:129-181,RGap:1-92,;gescore=67;splinfo
#  =-1116,nd4,sc83.058;aalen=127,69%,partial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1
# #d1.Scaffold264 splkf2eg7m1     mRNA    358922  388491  70      +       .       ID=Funhe2Eq7m000197t1;cov=70%,385/547;nex
#  on=5;splice=10;Target=Funhe2Eq7m000197t1 5 545;gaps=125,RGap:546-547,MGap:182-202,LGap:1-4,MGap:127-168,MGap:34-89,;gesco
#  re=53;splinfo=+1116,nd4,sc100.074;aalen=127,69%,partial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1
## >> 1 exon 98% aligned << GMAP picks this one; is  Transposon gene xpressed.
# #d2.Scaffold827 splkf2eg7m1     mRNA    75515   76058   98      -       .       ID=Funhe2Eq7m000197t1;cov=98%,534/547;nex
#  on=1;splice=0;Target=Funhe2Eq7m000197t1 1 546;gaps=1,RGap:547-547,;gescore=49;splinfo=+8227,nd4,sc19.502;aalen=127,69%,pa
#  rtial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1

	# ** BAD here? 1-exon lower cover scores above many-exon,higher cover model.. due to missing canon splices
	# .. maybe not bad, depends how much to weight canon-splices. 5-exon part align could be trash
	## >> bad case, 400 bp fragment likely TE-related, many valid aligns.. pick any
	# Scaffold2012:1096-16540:- FunhetEGm000040t3_G9;cov=58%,228/392;nexon=5;splice=12;Target=1 354;gaps=139;gescore=39;splinfo=-20302,nd9,sc120.736
  # Scaffold9953:1773835-1774051:- FunhetEGm000040t3_G6;cov=49%,193/392;nexon=1;splice=0;Target=1 213;gaps=179;gescore=44;splinfo=+24504,nd9,sc116.334
  # gmap got: Scaffold9941:1139485-1139880, 100% align, 1 exon
	## another problem case; gmap Scaffold10169:368500-368925:., 1exon
	# Scaffold151     11294   11509   47      +       ID=FunhetEGm000261t2;cov=47%,200/423;nexon=1;splice=4;Target=Fu
	# nhetEGm000261t2 147 395;gaps=174,RGap:1-146,LGap:396-423,;gescore=42;splinfo=-1038,nd9,sc121.769
	# #d1.Scaffold10107       555234  736229  58      -       ID=FunhetEGm000261t2;cov=58%,244/423;nexon=6;splice=6;T
	# arget=FunhetEGm000261t2 69 340;gaps=151,RGap:1-68,LGap:341-423,;gescore=23;splinfo=-26355,nd9,sc133.703
	# #d2.Scaffold10057       119821  1114429 72      +       ID=FunhetEGm000261t2;cov=72%,306/423;nexon=7;splice=2;T
	# arget=FunhetEGm000261t2 1 328;gaps=95,RGap:329-423,;gescore=19;splinfo=+12384,nd9,sc89.493
	## another .. gmap at Scaffold609:125478-125752:-, this is tiny frag. 
	# Scaffold609     124944  125176  52      +       ID=FunhetEGm000532t5;cov=52%,193/373;nexon=2;splice=4;Target=Fu
	# nhetEGm000532t5 1 203;gaps=170,LGap:204-373,;gescore=31;splinfo=-3552,nd3,sc103.482
	# #d1.Scaffold609 105775  125841  66      -       ID=FunhetEGm000532t5;cov=66%,248/373;nexon=7;splice=8;Target=Fu
	# nhetEGm000532t5 1 321;gaps=107,MGap:207-221,MGap:252-291,RGap:322-373,;gescore=28;splinfo=+17201,nd3,sc103.174;
	# #d2.Scaffold609 105775  125841  56      -       ID=FunhetEGm000532t5;cov=56%,209/373;nexon=5;splice=6;Target=Fu
	# nhetEGm000532t5 1 321;gaps=158,MGap:269-300,LGap:322-373,MGap:183-256,;gescore=25;splinfo=-17201,nd3,sc95.246;a

=cut

	## ^^ min spct for 0 splice is $nx / 5*nx .. 2/10, 3/15, 4/20 = .20;   
	## eg. hiscore 19 exons, 68 splscore, spct= .92 
	## gapspan included in pctalign;  
	# * $tscore? dont know what splign.log-tscore represents.. leave out?
	$genescore= int(0.5 + 100*$genescore);
	
	my $mattr="ID=$id;cov=$pctalign%,$maln/$mlen;nexon=$mnxon;splice=$msplscore$addat"; # $pctalign%, ?? dont need pctalign 2 places.. keep in col5
	$mattr.= ";Target=$trid $mrnab $mrnae"; # always add, so can measure partial aligns of mrna
	# $mattr.= ";trid=$trid" if($id ne $trid); # or Target=$trid 1 $mlen ?? add always??
	# $mattr.= ";splscore=$msplscore" if($mnxon>1); #? or always
	$mattr.= ";gaps=$gapspan,$gaps" if($gaps);	## prefix with gapspan and/or ngaps
	$mattr.= ";gescore=$genescore"; # insert in col5 instead of pctalign?

	# $mattr.= ";tscore=$tscore" if(1); # if($KEEP_SPLIGNSCORE);
	$mattr.= ";splinfo=$idir,nd$tsdup,sc$tscore" if($debug); # if($KEEP_SPLIGNSCORE);
	$mattr.= ";$trattr" if($trattr);
	
	## ?? replace print here with store in @gff, then filter out low-score 2ndary aligns per trid
	## use genescore AND mrnaspan to decide to keep, need mrnaspan parts for split-genes
	my $genomespan="$mref:$mgb-$mge:$mor"; # mor ?
	my $mrnaspan="$trid:$mrnab-$mrnae"; #?
	
	my $mrnagff= join("\t",$mref,$SRC,"mRNA",$mgb,$mge,$pctalign,$mor,".",$mattr)."\n";
	# print $OUTH join("\t",$mref,$SRC,"mRNA",$mgb,$mge,$pctalign,$mor,".",$mattr)."\n";
	$ntrout++;
	
	my @exongff=();
	foreach my $x (@xon) {
	  $x->[8] =~ s/Parent=$trid/Parent=$id/ if($id ne $trid);
	  $x->[8] =~ s/;align=.*$//; # dont need output? or leave to drop later?
	  $x->[6] = $mor if($flipor);
	  my $xout= join("\t",@$x)."\n"; push @exongff, $xout;
	  # print $OUTH join("\t",@$x)."\n";
	}
	
	# OPTION: output introns where adjacent exons have proper splices ..
	
	# FIXMEd: add CDS exons, given mRNA CDS offsets .. need more inputs
	my @cdsgff=();
	if($trattr =~ /offs=(\d+).(\d+)/) { # add cds-exons : option?  check $trattr=~/aalen=/ also ?
	  my($cb,$ce)=($1,$2); 
	  my $cdir= ($mor eq '-')?-1:1;
	  foreach my $x (@xon) {
	    my @c= @$x;
	    my($tb,$te)= $c[8] =~ m/Target=\S+ (\d+) (\d+)/;
	    if($te > $cb and $tb < $ce) {
	      $c[2]= "CDS";
	      ## fixme again; -strand needs c3,c4 swap also
        if($tb < $cb) { my $d=$cb-$tb; 
          if($cdir<0) { $c[4] -= $d; } else { $c[3] += $d; } # FIXME: - bad or offset; flip d sign, $or * $d
          }
        if($te > $ce) { my $d=$te-$ce; 
          if($cdir<0) { $c[3] += $d; } else { $c[4] -= $d; }
          }
	      $c[5]=1;
	      # $c[7]= $phase; ## FIXME: need phase calcs
	      $c[8] =~ s/;Target=.*//;
	  		my $xout= join("\t",@c)."\n"; push @cdsgff, $xout;
	      # print $OUTH join("\t",@c)."\n";
	    }
	  }
	}
	
	my @gene=($mrnagff,@exongff,@cdsgff);
	return ($trid,$genescore,$genomespan,$mrnaspan,\@gene); 
}


__END__

=item splign vs gmap

  info for kfish2p67vs.mrna gene set
  	.. overall pctalign is higher for splign
  	.. gmap for this mrna set has high proportion of split genes (chimera),
  		 but they are split at same locus, ie have gaps or mismatches to genome,
  		 and full span should be alignable .. maybe other gmap option will do that, dont know ..
  		 
  gmap nopath:14950  splign nomap:30106 of ntr=137660 ids
  .. check 'good' gmap chimera on 2 scaffolds versus splign
	
	split-gene FunhetEGm067171t5 failure in compart:
		- need to test other opts to get such, or replace w/ blastn  > splign..
		- this opt should pull many low-cover cases: -min_compartment_idty 0.50 => 0.20
					 
					 
	eg. spl-nomap hiqual
		FunhetEGm028320t1 4274,91%,complete	# gm: Scaffold10072:6500-8633,cov=15%  
			^^ end of scaf problem, split also?
		  ^^ there is more rna-seq at scaf:1-6500 should be gene asm parts; geno-align problem?
			^^ good full-length fish gene:  sacsin-like, 80% aln to 4284 aa, cichlids, tetraodon,xenopus,takifugu,..
		FunhetEGm067171t5 3577,62%,complete # gm split: Scaffold10067:190888-209862,cov=36%,3'part2 + Scaffold9902:670788-785051,cov=64,5'part
				# ^ this is true split-gene; note t1-4,t6 are mapped and longer **
		FunhetEGm025643t1 3222,61%,complete	# gm: Scaffold755:62292-137840,cov=52%
			^^ spans 75kb, with 2 smallish intronic genome gaps; other fish prots map at 52% here.
			>> there is 30kb gap at 5'end, start of scaf; this gene covers all of 150kb Scaffold755 (but top end, rev gene), 
			>> introns span gene.
			-- what part of tr/prot is not aligned? aug models are 1/2 size 1400 aa but same mapped exons.
			^^ good full-length fish gene: SZT2, 80% aln to 3698 aa, cichlids, tilapia,  ..
			>> Maylandia zebra, african cichlid is turning up top-hit for these, other kfish genes.
			# ^ these 3 not in compart output table, no align to genome?
		
		FunhetEGm018322t1 2941,93%,complete	# gm split: Scaffold630:10474-27099,cov=30,5'end + Scaffold848:4429-30540,cov=70,3'end
		FunhetEGm073739t1 2819,96%,complete # gm: Scaffold10038:1175876-1194030,cov=30
			^^ blastp align fragmented but full span for Maylandia, Oryzias fish top hits 2313aa, 40% cov/40% ident,  CDD:MDN1 AAA ATPase containing von Willebrand factor type A
			-- genome asm several big/small gaps both ends of mapped exons; missing 70% likely in 20kb gap.
			splign several alts map: FunhetEGm073739t13 cov=82%, 2/3 size; FunhetEGm073739t4 cov=75%, 2k-inc. 
						FunhetEGm073739t2 cov=88%, 2780 aa but exons spred beyond expected gene span.
		FunhetEGm058834t1 2741,90%,complete # gm: Scaffold10038:1365317-1389369,cov=46
			^^ blastp fragmented but full span for  Maylandia, tilapia fish, 3790aa, serine-rich adhesin for platelets-like
			splign 3 alts map here: FunhetEGm058834t8, 96% 1650aa-inc; FunhetEGm058834t10 ditto.. short aa
			.. problem alts: FunhetEGm058834t6, FunhetEGm058834t14 map past 5'end of t1..8: are these paralogs? evidence unclear
				.. prots are repetitive, but FunhetEGm058834t1 x FunhetEGm058834t6 only slight align.
			-- no genome gaps? ah, a large gap 30kb past 3' mapped end; 
			-- no other fish prots mapped here.. 
			
		FunhetEGm018012t1 2693,83%,complete	# gm: cov=47% Scaffold521:144238-196402 
			- several alts gm maps, spl misses; possible tandem dup gene, or one w/ exon mixup,
			- no compart align for FunhetEGm018012t*, but gmap got ~50% of many.
			- 52kb span, 1 small intronic scaf gap, but Scaffold521 ends at 3' exon .. end of scaf problem
			= 'Telomerase-associated protein 1' Funhe5EG032579t1 (1193aa), Funhe5EG032580t1(460 aa)
		FunhetEGm022232t1 2625,96%,partial3	# gm: cov=32% Scaffold9898:106970-129598 
			- strong express 5' end with messed up intron mapping: many overlaping, bi-dir.
			- Neuroblast differentiation-associated protein AHNAK,
			- splign maps FunhetEGm022232t2, 1180aa-inc, 100% cov; 
			- confused locus: several other trloci/diff prots, map here also; 3kb gap past 3' intron mess
			- expression suggests 2-3 long exons of 2k-5k, but introns map into these express exon spans.
		FunhetEGm009379t2 2605,73%,partial3	# gm: cov=32%
		
			* retest compart,splign on these cases of poor map (due to genogaps?), not splitscaf:
				FunhetEGm025643t1/Scaffold755, 
				FunhetEGm073739t1/Scaffold10038, FunhetEGm058834t1/Scaffold10038, 
				FunhetEGm018012t1/Scaffold521
					
			* blastp check proteins to see which of these are "real" full kfish genes?
			'FunhetEGm028320t1|FunhetEGm067171t5|FunhetEGm018322t1|FunhetEGm025643t1|FunhetEGm073739t1|FunhetEGm058834t1'
			
			* blastp check cases of gappy align to genome: is problem commoner w/ mrna or genome?
			
			
		FunhetEGm016085t4 2587,85%,complete	# t1 is mapped, very short 101aa, others t2-5 are 1500-2500aa
			# ^ what gives w/ evg tr2aacds calling short,partial 100aa FunhetEGm016085t1=kfish2qf7soapk21loc945695t1
			# t4 = kfish2qf7velvk71Loc2507t16
			# this is a NOMAIN, algo should have picked longest aa as t1: trevg67pt/publicset/kfish2p67vs.mainalt.tab
			# kfish2eg6velvk55Loc4503t9 NOMAIN  kfish2qf7soapk21loc945695t1/althi,kfish2eg6velvk55Loc4503t1/althi,
			#     kfish2eg6velvk55Loc4503t4/althi,kfish2qf7velvk71Loc2507t16/althi,kfish2eg6velvk61Loc4732t4/althi
 
=item retest missed compart

# test compart new opts for missed aligns:
cpdefault="-penalty 0.55 -min_idty 0.70"
cpopts="-penalty 0.25 -min_idty 0.25";  
test2: "-penalty 0.25" << no align
test3: "-min_idty 0.25"  << full align as 1st; use this min_idty (or lower?)
test4: "-min_idty 0.15" << same as .25
test5: "-min_idty 0.45" << same as .25

try this splign opt to match compart opt:
	splopt="-type mrna -min_compartment_idty 0.45"


CASE5: split FunhetEGm018322t1/s08split = Scaffold630,Scaffold848  
  compart cpopt=-min_idty 0.45
    >> only Scaffold848:4427-30540:+, tr:2832-9451
    
  compart cpopt=-penalty 0.25 -min_idty 0.25
    ++ adds Scaffold630:10474-27099:-, tr:1-2833
  
  compart cpopt=-penalty 0.25 -min_idty 0.45
    >> NO, same as 1.
    
  compart cpopt= -min_idty 0.25  
    ++ both, same as #2,  ie. skip -penalty, use low -min_idty
    
	** need new splign opts for this split case?  2nd part is lower score..
	a: splopt="-min_compartment_idty 0.45"
	b: splopt="-min_compartment_idty 0.25"  << this one to get 2nd part of split **
	   what of -compartment_penalty 0.55 def; No
	c: splopt="-min_compartment_idty 0.45 -compartment_penalty 0.25"   << no good
		"Multiple compartments will only be identified if they have at least this level of coverage"
			-- is that cover w/in the compart, or total cover of the mrna ?
	  -min_compartment_idty <Real, 0..1> Minimal compartment identity to align. Default = `0.7'

	
	$nbin/splign $splopt -comps $qfile1.compart8 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog8a -aln $qfile1.aln > $qfile1.splign8a
	#.....
		FunhetEGm018322t1.splog8a
		+1      FunhetEGm018322t1       Scaffold848     Ok      -504.258
		-1      FunhetEGm018322t1       Scaffold848     Ok      -805.335	

		FunhetEGm018322t1.splog8b
		+1      FunhetEGm018322t1       Scaffold848     Ok      -504.258
		-1      FunhetEGm018322t1       Scaffold848     Ok      -805.335
		+2      FunhetEGm018322t1       Scaffold630     Ok      -1193.53 << 2nd part here



CASE4: split  FunhetEGm067171t5/s27split = Scaffold10067,Scaffold9902 
	$nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  

  newspan: Scaffold9902:670788-801381/lend , Scaffold10067:108914-209862/rend
     159 FunhetEGm067171t5.mrna.compart5  : none orig run
  
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign

  FunhetEGm067171t5.mrna.splog
  +1	FunhetEGm067171t5	Scaffold10067	Ok	-1797.34  << best, 63/93 stdsplice,fwd-fwd
  -1	FunhetEGm067171t5	Scaffold10067	Ok	-1892.56  < worse, 23/89 stdsplice,rev-rev
  +2	FunhetEGm067171t5	Scaffold9902	Ok	-1920.22  << best, 69/92 stdsplice,fwd-fwd
  -2	FunhetEGm067171t5	Scaffold9902	Ok	-2175.28  <? 34/87 stdsplice, rev-rev


CASE3: FunhetEGm073739t1/Scaffold10038
	-min_idty 0.45; $nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  
	newspan: Scaffold10038:1175876-1216671:+  # gm is partial align here. this span just part of gene.. gaps
	
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign
		^^ splign result similar to t2+ alt mapping, ok, gaps due to geno gaps.
		
FunhetEGm073739t1.mrna.compart5
lcl|FunhetEGm073739t1   lcl|Scaffold10038       100     17      0       0       8631    8647    1216655 1216671 0       34
lcl|FunhetEGm073739t1   lcl|Scaffold10038       87.0968 31      4       0       8594    8624    1216621 1216651 0       62
lcl|FunhetEGm073739t1   lcl|Scaffold10038       96.2025 79      3       0       7109    7187    1214040 1214118 0       158
..
lcl|FunhetEGm073739t1   lcl|Scaffold10038       100     87      0       0       766     852     1176638 1176724 0       174
lcl|FunhetEGm073739t1   lcl|Scaffold10038       99.3141 717     4       0       34      750     1175909 1176625 0       1434
lcl|FunhetEGm073739t1   lcl|Scaffold10038       93.3333 30      2       0       1       30      1175876 1175905 0       60


CASE2: FunhetEGm018012t1/Scaffold521
	-min_idty 0.45; $nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  
	newspan:  Scaffold521:144238-188797:+ << this is gmap span
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign
	
FunhetEGm018012t1.mrna.compart5
lcl|FunhetEGm018012t1	lcl|Scaffold521	99	100	1	0	7410	7509	188698	188797	0	200
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	146	0	0	7264	7409	188279	188424	0	292
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	132	0	0	7010	7141	178766	178897	0	264
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	124	0	0	6885	7008	178507	178630	0	248
lcl|FunhetEGm018012t1	lcl|Scaffold521	98.8764	83	0	0	6797	6879	178313	178395	0	166
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	128	0	0	6668	6795	178097	178224	0	256
..
lcl|FunhetEGm018012t1	lcl|Scaffold521	99.4898	194	0	0	627	820	149559	149752	0	388
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	42	0	0	582	623	148245	148286	0	84
lcl|FunhetEGm018012t1	lcl|Scaffold521	98.6755	446	4	0	131	576	147791	148236	0	892
lcl|FunhetEGm018012t1	lcl|Scaffold521	93.0769	128	8	0	1	128	144238	144365	0	256


CASE1: FunhetEGm025643t1/Scaffold755
	newspan: Scaffold755:37619-137803:+
	
Lower opts:
$nbin/compart $cpopt -qdb FunhetEGm025643t1.mrna -sdb Scaffold755.gdna
	new compart opts gets Scaffold755:37619-77998,102301-137803 / tr:4511-8860,8863-15043
	^^ this is correct, introns cover this span, extends lower than gmap, catfish/zfish match full span.
	
Remap data created for Scaffold755.gdna; max offset = 146636
...
 Generating index volume: Scaffold755.gdna.424932660.p.v1 ... Ok
 Matching (strand = plus) ... 
lcl|FunhetEGm025643t1   lcl|Scaffold755 90      49      4       0       14995   15043   137755  137803  0       98
lcl|FunhetEGm025643t1   lcl|Scaffold755 82.2581 61      10      0       14933   14993   137688  137748  0       122
lcl|FunhetEGm025643t1   lcl|Scaffold755 94.1176 51      3       0       14870   14920   137624  137674  0       102
lcl|FunhetEGm025643t1   lcl|Scaffold755 95.8425 457     19      0       14398   14854   137151  137607  0       914
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.6974 1319    3       0       13077   14395   135830  137148  0       2638
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.8821 1648    0       0       11424   13071   134179  135826  0       3296
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     155     0       0       11165   11319   133979  134133  0       310
..
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.6109 253     0       0       8863    9115    102301  102553  0       506
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     185     0       0       8676    8860    77814   77998   0       370
..
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     128     0       0       5130    5257    51280   51407   0       256
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     284     0       0       4840    5123    37619   37902   0       568
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     158     0       0       4511    4668    3174    3331    0       316

Defaults:
$nbin/compart -qdb FunhetEGm025643t1.mrna -sdb Scaffold755.gdna
	new compart opts gets Scaffold755:37619-77998,102301-137803 / tr:4511-8860,8863-15043
 	.. 100% align for 37619-77998, 82%..100 for rest.
  gm: Scaffold755:62292-137840,cov=52%
 
$nbin/compart  -qdb $qfile1 -sdb $genome  | less
 Remap data created for Scaffold755.gdna; max offset = 146636
 Scanning 1 genomic sequences ... Ok
 Constructing FV ... Ok
 Remap data created for sequences; max offset = 15611
 Scanning sequences for N-mers and their positions.
 Generating index volume: qq.1505824399.p.v1 ... Ok
 Scanning Scaffold755.gdna for N-mers and their positions.
 Generating index volume: Scaffold755.gdna.1505824399.p.v1 ... Ok
 Matching (strand = plus) ... Ok
 Reading/transforming FV ... Ok
 Scanning sequences for N-mers and their positions.
 Generating index volume: qq.1505824399.m.v1 ... Ok
 Scanning Scaffold755.gdna for N-mers and their positions.
 Generating index volume: Scaffold755.gdna.1505824399.m.v1 ... Ok
 Matching (strand = minus) ... Ok

=cut

=item splign table

	spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splign
	
	+1	FunhetEGm046329t1	Scaffold0	0.994	344	1	344	2942702	2943045	  <exon>  	M89RMRM252   # msum= 252+1+89 = 342/344
	+2	FunhetEGm046569t1	Scaffold0	1	88	1	88	5831260	5831347	  <exon>GT	M88
	+2	FunhetEGm046569t1	Scaffold0	1	63	89	151	5831489	5831551	AG<exon>GT	M63
	+2	FunhetEGm046569t1	Scaffold0	1	212	152	363	5832998	5833209	AG<exon>  	M212
	+3	FunhetEGm047943t1	Scaffold0	-	5	1	5	-	-	<L-Gap>	-
	+3	FunhetEGm047943t1	Scaffold0	1	207	6	212	541156	541362	TT<exon>GT	M207
	+3	FunhetEGm047943t1	Scaffold0	1	425	213	637	544630	545054	AG<exon>GT	M425
	+3	FunhetEGm047943t1	Scaffold0	1	24	638	661	618363	618386	AG<exon>  	M24

=item splign logfile : useful?
	# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
	
	+1	FunhetEGm046329t1	Scaffold0	Ok	-8.248		: has only 1 exon partly aligned, thus low -score
	+2	FunhetEGm046569t1	Scaffold0	Ok	-17.252
	+3	FunhetEGm047943t1	Scaffold0	Ok	-28.458
	+4	FunhetEGm047950t1	Scaffold0	Ok	-28.65
	+5	FunhetEGm047996t1	Scaffold0	Ok	-16.632
	+6	FunhetEGm048121t1	Scaffold0	Ok	-18.912
	+7	FunhetEGm048618t1	Scaffold0	Ok	-151.9
	+8	FunhetEGm057611t1	Scaffold0	Ok	-333.571
	+9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: + = forward
	-9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: - = reverse

=item eg FunhetEGm056896t1
	# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
	-- several mappings, probably take highest -score, skip any 2ndary below <99% of top?
	+9	FunhetEGm056896t1	Scaffold0	Ok	-44.051
	-9	FunhetEGm056896t1	Scaffold0	Ok	-44.051
	+875	FunhetEGm056896t1	Scaffold348	Ok	-60.039    < best score, 
	-875	FunhetEGm056896t1	Scaffold348	Ok	-56.914
	+3190	FunhetEGm056896t1	Scaffold9988	Ok	-49.873
	-3190	FunhetEGm056896t1	Scaffold9988	Ok	-49.873
	+5626	FunhetEGm056896t1	Scaffold936	Ok	-47.645
	-5626	FunhetEGm056896t1	Scaffold936	Ok	-47.645
	+8236	FunhetEGm056896t1	Scaffold10172	Ok	-37.791
	-8236	FunhetEGm056896t1	Scaffold10172	Ok	-37.791
	
	# are these horrid long mismatch strings useful?
	#   M10RM9RM9RM16RMRM12RM23RM11RM6RM6RM51RM2I8MRM3IMRM2RMRMRMRM2D2M2RM33RM8RM3RM37RM62
	+875	FunhetEGm056896t1	Scaffold348	0.905	346	1	337	76738	77081	  <exon>GG	M10RM9RM9RM16RMRM12RM23RM11RM6RM6RM51RM2I8MRM3IMRM2RMRMRMRM2D2M2RM33RM8RM3RM37RM62
	+875	FunhetEGm056896t1	Scaffold348	-      36	338	373	-	-	<M-Gap>	-
	+875	FunhetEGm056896t1	Scaffold348	0.905	231	374	603	77082	77306	GA<exon>GT	M6RM15RM28RM24RM42DM41RM6IM2RMR2MD5M3R3M3R3M14RM23
	+875	FunhetEGm056896t1	Scaffold348	0.933	360	604	963	79954	80313	GG<exon>  	M16RM18RM28RM43RM16RM79RMRM12RM16RM2RM3R2MRM4RM2RM2RM20RM8RM5RM2RM5RM2RM12RM22RM17
	
	# .. rev map is slightly different: 77082 .. 77306,77313 
	-875	FunhetEGm056896t1	Scaffold348	0.932	353	963	611	80313	79961	  <exon>GT	M17RM22RM12RM2RM5RM2RM5RM8RM20RM2RM2RM4RMR2M3RM2RM16RM12RMRM79RM16RM43RM28RM18RM9
	-875	FunhetEGm056896t1	Scaffold348	0.908	238	610	374	77313	77082	TG<exon>TC	M30RM14R3M3R3M3R2M2RD5M3IM5RM41DM42RM24RM28RM15RM6
	-875	FunhetEGm056896t1	Scaffold348	-	     36	373	338	-	-	<M-Gap>	-
	-875	FunhetEGm056896t1	Scaffold348	0.905	346	337	1	77081	76738	CC<exon>  	M62RM37RM3RM8RM33RM2D2M2RMRMRMRM2RMIM3RM2I8MRM51RM6RM6RM11RM23RM12RMRM16RM9RM9RM10

=item eg. +8	FunhetEGm057611t1	Scaffold0	Ok	-333.571

	+8	FunhetEGm057611t1	Scaffold0	-	198	1	198	-	-	<L-Gap>	-
	+8	FunhetEGm057611t1	Scaffold0	1	176	199	374	2463919	2464094	AA<exon>GT	M176
	+8	FunhetEGm057611t1	Scaffold0	1	131	375	505	2465439	2465569	AG<exon>GT	M131
	+8	FunhetEGm057611t1	Scaffold0	1	505	506	1010	2468703	2469207	AG<exon>GT	M505
	+8	FunhetEGm057611t1	Scaffold0	1	771	1011	1781	2472316	2473086	AG<exon>GT	M771
	+8	FunhetEGm057611t1	Scaffold0	1	195	1782	1976	2474583	2474777	AG<exon>GT	M195
	+8	FunhetEGm057611t1	Scaffold0	1	164	1977	2140	2476325	2476488	AG<exon>GT	M164
	+8	FunhetEGm057611t1	Scaffold0	1	162	2141	2302	2476610	2476771	AG<exon>GT	M162
	+8	FunhetEGm057611t1	Scaffold0	1	77	2303	2379	2476868	2476944	AG<exon>GT	M77
	+8	FunhetEGm057611t1	Scaffold0	1	223	2380	2602	2478033	2478255	AG<exon>GT	M223
	+8	FunhetEGm057611t1	Scaffold0	1	234	2603	2836	2478338	2478571	AG<exon>GT	M234
	+8	FunhetEGm057611t1	Scaffold0	1	92	2837	2928	2478671	2478762	AG<exon>GT	M92
	+8	FunhetEGm057611t1	Scaffold0	1	128	2929	3056	2478935	2479062	AG<exon>GT	M128
	+8	FunhetEGm057611t1	Scaffold0	1	155	3057	3211	2479950	2480104	AG<exon>GT	M155
	+8	FunhetEGm057611t1	Scaffold0	1	135	3212	3346	2480194	2480328	AG<exon>GT	M135
	+8	FunhetEGm057611t1	Scaffold0	1	195	3347	3541	2480501	2480695	AG<exon>GT	M195
	+8	FunhetEGm057611t1	Scaffold0	1	183	3542	3724	2480929	2481111	AG<exon>GT	M183
	+8	FunhetEGm057611t1	Scaffold0	1	205	3725	3929	2481201	2481405	AG<exon>GT	M205
	+8	FunhetEGm057611t1	Scaffold0	1	147	3930	4076	2481498	2481644	AG<exon>GT	M147
	+8	FunhetEGm057611t1	Scaffold0	1	143	4077	4219	2481731	2481873	AG<exon>GT	M143
	+8	FunhetEGm057611t1	Scaffold0	1	204	4220	4423	2482561	2482764	AG<exon>GT	M204
	+8	FunhetEGm057611t1	Scaffold0	1	118	4424	4541	2482958	2483075	AG<exon>GT	M118
	+8	FunhetEGm057611t1	Scaffold0	1	103	4542	4644	2483719	2483821	AG<exon>GT	M103
	+8	FunhetEGm057611t1	Scaffold0	0.994	464	4645	5108	2483953	2484415	AG<exon>AA	M425RM7RM4DM25

=cut


=item genosplign.sh cluster script

	see evigene/scripts/rnaseq/genosplign.sh
	-------------
	#! /bin/bash
	### env genome=xxx mrna=mrna.fa datad=`pwd` qsub -q normal genosplign.sh
	# NCBI splign align mRNA to genome
	splopt="-type mrna -min_compartment_idty 0.5"
	
	# split mrna.fasta to ncpu parts, run ncpu cluster jobs, output to spl$i part folders
	qset=`/bin/ls $mrna.split.*.fa`
	i=0; for qfile in $qset; do {
		# setup inputs ..
	  # this is fork set:
  	echo "# splign $qfile1 x $genome in $ispldir";
  	( $nbin/splign -mklds ./; \
  	  ln -s ../$genome.* ./ ; \
      $nbin/makeblastdb -parse_seqids -dbtype nucl -in $qfile1; \  
      $nbin/compart  -qdb $qfile1 -sdb $genome  > $onam.cpart; \  
      $nbin/splign $splopt -comps $onam.cpart -blastdb $genome -ldsdir ./ -log $onam.splog > $onam.splign; ) &
		i=$(( $i + 1 ))
	} done
	wait
	#.................. 

=cut

=item compart options (blastn variant)


bin/compart -help
USAGE
  /bio/bio-grid/mb/ncbix/bin/compart [-h] [-help] [-xmlhelp]
    [-qdb qdb] [-sdb sdb] [-ho] [-penalty penalty] [-min_idty min_idty]
    [-min_singleton_idty min_singleton_idty]
    [-min_singleton_idty_bps min_singleton_idty_bps] [-max_intron max_intron]
    [-dropoff dropoff] [-min_query_len min_query_len]
    [-min_hit_len min_hit_len] [-maxvol maxvol] [-noxf] [-seqlens seqlens]
    [-N N] [-version-full] [-dryrun]

DESCRIPTION
   Compart v.1.35. Unless -qdb and -sdb are specified, the tool expects
   tabular blast hits at stdin collated by query and subject, e.g. with 'sort
   -k 1,1 -k 2,2'

OPTIONAL ARGUMENTS
 -h	 Print USAGE and DESCRIPTION;  ignore all other parameters
 -help	 Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -xmlhelp	Print USAGE, DESCRIPTION and ARGUMENTS in XML format; ignore all other parameters
 -qdb <String>	 cDNA BLAST database
 -sdb <String>	 Genomic BLAST database
 -ho		Print raw hits only - no compartments
 -penalty <Real, 0..1>	Per-compartment penalty   Default = `0.55'
 -min_idty <Real, 0..1>
   Minimal overall identity. Note: in current implementation  there is no
   sense to set different 'min_idty' and 'min_singleton_idty' (minimum is used
   anyway).
   Default = `0.70'
 -min_singleton_idty <Real, 0..1>
   Minimal identity for singleton compartments. The actual parameter passed to
   the compartmentization procedure is least of this parameter multipled by
   the seq length, and min_singleton_idty_bps. Note: in current implementation
    there is no sense to set different 'min_idty' and 'min_singleton_idty'
   (minimum is used anyway).
   Default = `0.70'
 -min_singleton_idty_bps <Integer>
   Minimal identity for singleton compartments in base pairs. Default =
   parameter disabled.
   Default = `9999999'
 -max_intron <Integer>
   Maximum intron length (in base pairs)
   Default = `1200000'
 -dropoff <Integer>
   Max score drop-off during hit extension.
   Default = `5'
 -min_query_len <Integer, 21..99999>
   Minimum length for individual cDNA sequences.
   Default = `50'
 -min_hit_len <Integer, 1..99999>
   Minimum length for reported hits in hits-only mode. No effect in
   compartments mode.
   Default = `16'
 -maxvol <Integer, 128..1024>
   Maximum index volume size in MB (approximate)
   Default = `512'
 -noxf
   [With external hits] Suppress overlap x-filtering: print all compartment
   hits intact.
 -seqlens <File_In>
   [With external hits] Two-column file with sequence IDs and their lengths.
   If none supplied, the program will attempt fetching the lengths from
   GenBank. Cannot be used with -qdb.
 -N <Integer>
   [With external hits] Max number of compartments per query (0 = All).
   Default = `0'
 -version-full
   Print extended version data;  ignore other arguments
 -dryrun
   Dry run the application: do nothing, only test all preconditions


=item splign options

bin/splign -help
USAGE
  /bio/bio-grid/mb/ncbix/bin/splign [-hits hits]
    [-comps comps] [-mklds mklds] [-blastdb blastdb] [-ldsdir ldsdir]
    [-query query] [-subj subj] [-disc] [-W mbwordsize] [-type type]
    [-compartment_penalty compartment_penalty]
    [-min_compartment_idty min_compartment_identity]
    [-min_singleton_idty min_singleton_identity]
    [-min_singleton_idty_bps min_singleton_identity_bps]
    [-min_exon_idty identity] [-min_polya_ext_idty identity]
    [-min_polya_len min_polya_len] [-max_intron max_intron]
    [-max_space max_space] [-direction direction] [-log log] [-asn asn]
    [-aln aln]

DESCRIPTION
   Splign: 1.39.8
   Package: public 12.0.0, build Jul 18 2013 16:30:51
    LLVMGCC_421-Debug64--x86_64-apple-darwin12.4.0-c_67_173_144_85
    
OPTIONAL ARGUMENTS
 -hits <File_In>
   [Batch mode] Externally computed local alignments (such as blast hits), in
   blast tabular format. The file must be collated by subject and query (e.g.
   sort -k 2,2 -k 1,1).
 -comps <File_In>
   [Batch mode] Compartments computed with Compart utility.
 -mklds <String>
   [Batch mode] Make LDS DB under the specified directory with cDNA and
   genomic FASTA files or symlinks.
 -blastdb <String>
   [Batch mode] Blast DB.
 -ldsdir <String>
   [Batch mode] Directory holding LDS subdirectory.
 -query <File_In>
   [Pairwise mode] FASTA file with the spliced sequence.
 -subj <File_In>
   [Pairwise mode] FASTA file with the genomic sequence.
 -disc
   [Pairwise mode] Use discontiguous megablast to facilitate alignment of more
   divergent sequences such as those from different organisms (cross-species
   alignment).
 -W <Integer>
   [Pairwise mode] Megablast word size
   Default = `28'
 -type <String, `est', `mrna'>
   Query cDNA type: 'mrna' or 'est'
   Default = `mrna'
 -compartment_penalty <Real, 0..1>
   Penalty to open a new compartment (compartment identification parameter).
   Multiple compartments will only be identified if they have at least this
   level of coverage.
   Default = `0.55'
 -min_compartment_idty <Real, 0..1>
   Minimal compartment identity to align.
   Default = `0.7'
 -min_singleton_idty <Real>
   Minimal singleton compartment identity to use per subject and strand,
   expressed as a fraction of the query's length.
 -min_singleton_idty_bps <Integer>
   Minimal singleton compartment identity to use per subject and strand, in
   base pairs. The actual value passed to the compartmentization procedure is
   the least of (min_singleton_idty * query_length) and
   min_singleton_identity_bps.
   Default = `9999999'
 -min_exon_idty <Real, 0..1>
   Minimal exon identity. Segments with lower identity will be marked as gaps.
   Default = `0.75'
 -min_polya_ext_idty <Real, 0..1>
   Minimal identity to extend alignment into polya. Polya candidate region on
   mRNA is detected first. Alignment is produced without the polya candidate
   region After that alignment will be extended into the polya candidate
   region to deal with case when initial polya detection was wrong
   Default = `1'
 -min_polya_len <Integer, 1..1000000>
   Minimal length of polya.
   Default = `1'
 -max_intron <Integer, 7..2000000>
   The upper bound on intron length, in base pairs.
   Default = `1200000'
 -max_space <Real, 500..4096>
   The max space to allocate for a splice, in MB. Specify lower values to
   spend less time stitching over large genomic intervals.
   Default = `4096'
 -direction <String, `antisense', `auto', `both', `default', `sense'>
   Query sequence orientation. Auto orientation begins with the longest ORF
   direction (d1) and proceeds with the opposite direction (d2) if found a
   non-consensus splice in d1 or poly-a tail in d2. Default translates to
   'auto' in mRNA and 'both' in EST mode
   Default = `default'
 -log <File_Out>
   Splign log file
   Default = `splign.log'
 -asn <File_Out>
   ASN.1 output file name
 -aln <File_Out>
   Pairwise alignment output file name
   
=cut