#!/usr/bin/perl
# altbest.pl : pull valid alt-transcript assemblies

=item about

 altbest.pl : pull valid alt-transcript assemblies
    input: bestgenes.gff + all_altasm.gff
    input exons must be marked with intr=..N1,N2, intron splice ids
    input gff must be gene-record ordered (mRNA,exons)..
    
 intron chain annotation:
   $workd/scripts/overlapfilter -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
   -mark intr -over $workd/intron/intron_good.gff.gz -in genes.gff \

 steps: 
   1. equivalence allalt.gff to bestgenes (CDS equiv >= 66 %?)
   2. pull intron chains of equivalent alts, filter only uniq, drop short chains

=item fixmees

	FIXME3: Nopath (unmapped) alts are kicked out, should not do that .. but no introns for such
			.. however, qlen=bogus from gmap-nopath, needed below.
			.. also may want to keep alts w/ mapto -introns but not annotated valid intr=
			
	** FIXME2: 
	a. done? revise to handle valid alternate exon ends (e.g. 1 alt continues w/ intron, 2nd ends cds)
	
  b. option to use genes.gff introns for valid intron chains, instead/with evd-intron annots.
     ^^b. do here?, or let caller make intron.tab and annotate genes, per above, would be simpler.
          i.e. gmap.out, gmap.gff output has valid intron splices for trgenes, use that?
  
  b2. input introns.gff/table instead of requiring genes.gff exon intr= annotation
      .. need parts from overlapfilter -intron2splice=error
          
 ** FIXME altbest.pl : reduce for partial prot vs shorter full prot of ~same rna models
 .. problem from pasa_asm > longest-orf-partials not as good as full prot alt models.
   
=cut

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use Getopt::Long;

our $EVIGENES="$FindBin::Bin";  # $evigene="/bio/bio-grid/mb/evigene/scripts";
my $MINID_CDS=66;
my $MINID_UTR=40; # this is really full transcript span percent, utr+cds
my $MIN_PCDS=40;  # 0=all, assume mRNA cxlen=cds,exon or calculate? use to filter falseutr models
my $MIN_CDS2MAIN= 0.55; # relative to main pcds, ie drop false long utr aberrations
my $MIN_ALTCHAINSIZE = 0.30; # OPTION: min proportion of introns in short alt, of longest
my $TAG_ALT='t'; # append/replace to prime id
my $CHANGEALTID=1; 
my ($maingenes, $altgenes, $eqtab, $altab, $debug,$keepnointron);

my $optok= GetOptions(
  "altgenes=s", \$altgenes, 
  "maingenes=s", \$maingenes,  
  "eqtab=s", \$eqtab, 
  "mincds=i", \$MINID_CDS,  
  "minexon|minutr=i", \$MINID_UTR,  
  "minintronchainsize=s", \$MIN_ALTCHAINSIZE,  
  "mincoding=i", \$MIN_PCDS,  
  "minpcds2main=s", \$MIN_CDS2MAIN,  
  "keepnointron!", \$keepnointron, # is this flag or value option? to say keep NOPATH/others w/o intr= annots
  "tagalt=s", \$TAG_ALT, 
  "CHANGEALTID!", \$CHANGEALTID,  # noCHANGEALTID turns off, default on?
  "debug!", \$debug, 
  ## separate debug opt to write invalidchain to ATAB ?
  );

die "usage: altbest.pl -main cacao11_bestgenes.gff -alt asm.all.gff
   options: -mincds=$MINID_CDS -minexon=$MINID_UTR -mincoding=$MIN_PCDS -nochangealtid
Note: input gff must have exon intr=..,N1,N2,N3 intron chain annots.
" unless($optok and (-f $maingenes and -f $altgenes));

$MIN_CDS2MAIN=$MIN_CDS2MAIN/100 if($MIN_CDS2MAIN>1);
$MIN_ALTCHAINSIZE=$MIN_ALTCHAINSIZE/100 if($MIN_ALTCHAINSIZE>=1);

# change to -type similarCDS / utr and filter input by min cds/utr
# my $eqcmd="$EVIGENES/overgenedup.pl -exon=CDS -type CDS -mincds=$MINID_CDS -slopexon=8 -act markid -mark overg "
# ." -in $maingenes -over $altgenes";

my $eqcmd="$EVIGENES/overgenedup.pl -type similarCDS -mincds=$MINID_CDS  -minutr=$MINID_UTR -act markid -mark overg "
." -in $maingenes -over $altgenes";

unless($eqtab) {
	$eqtab=$maingenes; $eqtab=~s/.gz//; $eqtab =~ s/.gff//; $eqtab .= ".equalalt.tab";
}

$altab= $eqtab; $altab =~ s/.equalalt.*.tab//; $altab.= ".altbest.tab";
my $valtgff= $altab; $valtgff =~ s/.tab//; $valtgff .=".gff";

# ** FIXME: in busy regions get some bizarre assemblies that cover several genes,
# .. with long false UTRs; filter somewhere here, using multi-maingene overlaps?

# -- measure cdslen, exonlen and drop out false utr where cds/exon < 0.3,0.4 ?

warn "#options: -mincds=$MINID_CDS -minexon=$MINID_UTR -minintronchain=$MIN_ALTCHAINSIZE"
	." -mincoding=$MIN_PCDS -codingOfMain=$MIN_CDS2MAIN -tag=$TAG_ALT\n" if $debug;

# make equiv table
unless( -f $eqtab) {
	warn "# equalalt $eqtab = $eqcmd\n" if $debug;
	open(EOUT,">$eqtab") or die "writing $eqtab";
	open(ECMD,"$eqcmd |") or die "error $eqcmd";
	while(<ECMD>) {
		if(/\tmRNA/) { 
		my($g)=m/ID=([^;\s]+)/; my($ov)= (m/overg=([^;\s]+)/) ? $1 : ""; next unless($ov);
		my $og= (m/oid=([^;\s]+)/) ? $1:"noid"; 
		my($r,$b,$e)=(split"\t")[0,3,4]; $r=~s/^Scaffold[_]?(\S+)/${1}sc/i; 
		print EOUT join("\t",$g,$og,$ov,"$r:$b-$e"),"\n" if($ov);
		}
	} 
	close(ECMD);  close(EOUT);
}


my (%gmain,%galt,%ichain,%allalt2main,%mainpcds,%gcds,%exonends);
my ($nerrmrna,$ninmrna,$nxin,$nxintr)=(0) x 10;

# read equiv table
open(EOUT,$eqtab) or die "open $eqtab";
while(<EOUT>){
  my($gid,$altids)=(split"\t")[0,2];
  # my @altids= map{ s,/.*$,,; $_ } split",",$altids;
  #v2: parse  ID/[CI]99.88,
  my @altids= grep /\w/, map{ 
    my($c,$x)=m,/[CI]?(\d+)[.]?(\d*),; s,/.*$,,;
    $c||=0; $x=$c unless($x);
    $_ if($c>=$MINID_CDS and $x>=$MINID_UTR); 
    } split",",$altids;
  #? should this track/log altids not retained in galt ? ie, cds/utr   <min 
  if(@altids) {
    $galt{$gid}= \@altids; # ok? eqtab gid col is uniq
    map{ $allalt2main{$_}= $gid; } @altids; # FIXME: alts map onto 2+ gids ; last not best **
  }
}

# collect exon intron chains, only for equiv genes
foreach my $gf ($maingenes, $altgenes) {
  my $ismain= ($gf eq $maingenes)?1:0;
  my $gop = ($gf =~ /.gz/) ? "gunzip -c $gf" : "cat $gf"; 
  # drop exon filter, test mRNA attribs.
  
  my ($gid,$keep,%gx,$iexon,$nexon); 
  open(G,"$gop |") or die "bad $gop";
  while(<G>){
    next unless(/^\w/);
    if(/\tmRNA/) {
      if($keep and $gid) {  ## FIXME3 here? for NOPATH, no intron alts; other mapped+nointron alts?
				my @xin=();
				if($gx{$gid}) { @xin= map{ $gx{$gid}{$_} } sort{ $a<=>$b } keys %{$gx{$gid}}; }
        # if($keepnointron and not @xin) { } # fake it? but how so we get some useful alts, but not lots o crap      
        $ichain{$gid}= join "; ",@xin if @xin; 
   		}
        
      ($gid)= m/ID=([^;\s]+)/; $ninmrna++;
      $keep=1; %gx= ();   #NO cant reset here: %gcds=();
      $iexon=0; ($nexon)= m/nexon=(\d+)/; $nexon||=1;
      
      ## FIXME: mRNA fields
      ## gmap: ID=Funhe2Exy3m110884t1;trg=Funhe2Exy3m110884t1 1 292;
      ## aalen=51;cov=93.0;indels=3/5;match=278;nexon=1;pid=93.6;qlen=314;cdsindel=-2;
      ## aalen=63,61%,complete;offs=17-208;oid=Fungr1EG3m041003t1
      ## fix2: 900/153061 with chimera=breakpoint.. are missing aalen= from gmap
      ## FIXME3: nopath alts should be kept, but have no introns, invalid clen=99 false
      # NOPATH  kf2xx9gmap  mRNA  1   69  .. ID=Funhe2Exx11m004617t1;altc=althi1;trg=Funhe2Exx11m004617t1 1 69;
      #	 cov=0;match=0;pid=1;qlen=99;path=0/0;cdsindel=1866;aalen=803,94%,complete;offs=4-2415;oid=Funhe2E6bm004599t6
      # NOPATH  kf2xx9gmap  mRNA  1   69  .. ID=Funhe2Exx11m004617t2;altc=main;trg=Funhe2Exx11m004617t2 1 69;
      #  cov=0;match=0;pid=1;qlen=99;path=0/0;cdsindel=691;aalen=847,92%,complete;offs=84-2627;oid=Funhe2E6bm004599t1
      
      my($cw,$xw)= m/cxlen=(\d+).(\d+)/; ## BAD expectation .. deal w/ missing field
    	unless($xw) { ($xw)=m/(?:qlen|clen)=(\d+)/;  }
     	unless($cw) { 
      	## my($ofs)=m/offs=([\d-]+)/; if($ofs) { my($ob,$oe)=split/-/,$ofs; $cw=1+$oe-$ob; }
      	my($aw)=m/aalen=(\d+),\d/; 
      	($aw)=m/aalen=(\d+)/ unless($aw); 
      	$cw=3*$aw if($aw); 
      }
      $nerrmrna++ unless($cw and $xw);
      die "ERR: Missing mRNA annots: cxlen= or clen/qlen= and aalen= for $nerrmrna/$ninmrna at ID=$gid in $gf"
      	if($nerrmrna/$ninmrna > 0.2);
      	
      #? FIXME here: check for alts mapped to diff locus than main ? should remove such from alt test.
       	
      if( $cw and $xw ) { 
        my $pcds= ($xw>0) ? int(0.5 + 100*$cw/$xw) : 100; 
        if($ismain) { $mainpcds{$gid}= $pcds; }
        else {
          my $mainid= $allalt2main{$gid}; # FIXME: 2+ mainids possible?
          my $mcds= ($mainid) ? $mainpcds{$mainid} : 100; $mcds||=0; # error?
          $keep=0 if($pcds < $MIN_PCDS or $pcds < $MIN_CDS2MAIN * $mcds); # maybe should filter relative to maingene pcds< YES
          # FIXME? record not-keeps IDs/reason
          }
        }
      
    } elsif(/\texon/ and $keep) { # assumes now is gene record sorted
    	## for altexonend, may need CDS end points also: $gcds{$pid}="begin,end" ? or exons b,e?
      my($pid)= m/Parent=([^;\s]+)/ or next; 
      my($gref,$b,$e)=(split"\t")[0,3,4];  $nxin++;  $iexon++;
      # my($in)= m/;intr=[^,]+,(N[^;\s]+)/; 
      my($inscore,$in)= m/;intr=([^,]+),(N[^;\s]+)/; ## inscore=123 ; =+33/-22 ; =whatelse?
      $gmain{$pid}++ if($ismain);
    
      ## Change here? for altexonend: add exon b,e to each in found .. to the right one ! that we dont know from inID
      ## ;intr=+2041/-3,N204822,N204851,N204852,N204853
      # odd intr bug: N25789,N25789 >> N25785; N25785,N25786; N25786,N25787; N25787,N25788; N25788,N25789,N25789
      if($in) { 
        my($intscore)= $inscore=~m/([+-]?\d+)/; $intscore ||= 1;
        my @in=split",",$in; my %din; @in=grep{ !$din{$_}++ } @in; 
        if(@in) { 
          my $inkey = join ",",@in; 
          $gx{$pid}{$b} = $inkey;
          $gcds{$pid}{$inkey}="$b,$e"; 
          $exonends{$inkey}{$b}+= $intscore; # not ++ count,
          $exonends{$inkey}{$e}+= $intscore;  # dont reset exonends for each mRNA
          $nxintr++; 
          }          
     	} elsif($keepnointron) { # fake something ?
     		#? Fix here for too many trivial alts of mapped tr : maybe exclude end exons (but for NOPATH)
     		# .. seems to work, altbest count about 1/2 between -nokeepnoin and last -keepnoin
     		## naltbest, nokeep=77131/nopath=0, lastkeep=143763/nopath=15130, thiskeep=102248/nopath=14703
     		## FIX-gsplign has lots of false long,middle introns for poormap cases; should not allow those as keepers here..
     		## my $localkeepnointron= ($keepnointron and not ($src =~ /$GSPLIGN_SOURCE/))?$keepnointron:0;
     		## .. or find other way for this; add some okintron annots to input.gff ?
     		
        my $inkey = "NOINT".$gref; ## ."xe$b$e"; # for NOPATH, b,e are all same
        my $intscore=1; my $keepnoin=1; 
        if($gref =~ /^NOPATH/) { $inkey = "NOINT".$pid; } else { $keepnoin=0 if($iexon==1 || $iexon>=$nexon); }
        if($keepnoin) {
        $inkey.= ($iexon==1)? "xb$e" : ($iexon>=$nexon)? "xe$b" : "x$b$e";
        $gx{$pid}{$b} = $inkey;
        $gcds{$pid}{$inkey}="$b,$e"; 
        $exonends{$inkey}{$b}+= $intscore; # not ++ count,
        $exonends{$inkey}{$e}+= $intscore;  # dont reset exonends for each mRNA
        }
      }
# #intronchains for mainID=Funhe2Exx11m004617t1, mainIn=0 == NOPATH w/ valid alts
# invalidchain    Funhe2Exx11m004617t1    z2      dropshort,0,    NOINTFunhe2Exx11m004617t1xe169
# invalidchain    Funhe2Exx11m004617t2    z3      dropshort,0,    NOINTFunhe2Exx11m004617t2xe169
#...      	
    }
  } close(G);
  
	if($keep and $gid) {  ## FIXME3 here? for NOPATH, no intron alts; other mapped+nointron alts?
		my @xin=();
	  if($gx{$gid}) { @xin= map{ $gx{$gid}{$_} } sort{ $a<=>$b } keys %{$gx{$gid}}; }
		# if($keepnointron and not @xin) { }# fake it? but how so we get some useful alts, but not lots o crap		
		$ichain{$gid}= join "; ",@xin if @xin; 
	}

}

warn "# intron chains: $ninmrna mRNA input, $nerrmrna errors skipped, $nxintr/$nxin exons w/ introns, from main, alt.gff\n" if $debug;


my %valclass=( 
	 1 => "okay", 2 => "okaltend", 3 => "oknewlocus",
	-1 => "dropdup", -2 => "dropshort", -3 => "dropsubchain", 0 => "drop" );
my (%altnum, %galtok, %galtsame, %alt2gid, %nclass, $ntabmain, $ntabalt);

warn "# $altab: valid-alts for main=$maingenes from alts=$altgenes\n" if $debug;
open(ATAB,">$altab") or die "$altab"; # STDOUT?
print ATAB "#altbest input main=$maingenes, alts=$altgenes\n";
foreach my $gid (sort keys %galt) {
  my (%did,@did,%inset);
  my @alt= @{$galt{$gid}} or next;
  my $dn1= $ichain{$gid} or next; # FIX here for no introns.
  my @dns= sort{ length($ichain{$b}) <=> length($ichain{$a}) or $a cmp $b } grep{ $ichain{$_} } @alt;
  my $anum=1; 
  my $invalnum=0;
  my $i1= $dn1=~tr/;/;/; 
  my $imincount = int($MIN_ALTCHAINSIZE*$i1); # if long = 10 introns, < 3 is too short; can be 0 for keepnointron 
  $imincount||=1 unless($keepnointron);
	$ntabmain++;
	%inset=(); 
	map{ $inset{$_}++ } split /\W+/,$dn1; ## add gid to inset to track old/new loci?
	#? map{ $inset{$_}= $gid } split /\W+/,$dn1; ## add gid to inset to track old/new loci?

  print ATAB "#intronchains for mainID=$gid, mainIn=$i1\n"; # needs to be parsable to link altids 
  foreach my $aid ($gid, @dns) { 
    my $ichain=$ichain{$aid}; my $inum= $ichain=~tr/;/;/; 
    my $ismain= ($aid eq $gid)?1:0;
    
    ## oknewlocus: handle like new gid here?  need to reset dn1, imincount, inset ..
    ## gets too messy .. should separate out newloci from alt.gff before this script.
# if(0) {    
#     my $gmain= $gid;
#     my %ins=(); map{ $ins{$_}++ } split /\W+/,$ichain; 
#     my %imain=(); map{ my $imain=$inset{$_}; $imain{$imain}++ if($imain); } keys %ins;
#     my($imainid)= sort{ $imain{$b} <=> $imain{$a} } keys %imain; # gmain replaces gid?
#    	if($gmain ne $imainid) { # $mainset{$imainid} and 
#    		$ismain=1; $gmain= $imainid; $dn1= $ichain; 
#  			$imincount = int($MIN_ALTCHAINSIZE*$inum); # if long = 10 introns, < 3 is too short
# 		}
# }		
		
    ## FIXME: here to keep alternate end exons: intrchain may be same as others, or may be short,
    ##   but some are valid alts, esp. when end exon is longish.
    ## need examples to check params. record ic-invalids ? check for alt-drops that have large cds/aalen
    
    my $icsameasmain= (!$ismain and $ichain eq $dn1)?1:0;
			#^ does icsameasmain imply icvalid=-1/-4 ? set  $did{$dn1}=1 ? 
			#^ ah, for aid in (gid,dns) does this.
			
    ## these lowqual intronchain tests may need adjust: $inum < int(0.3*$i1) or ichain subset in @did
    my $icvalid = 1; 
    $icvalid=-1 if($icvalid==1 && $did{$ichain});  # ichain is already done
    $icvalid=-2 if($icvalid==1 && ($inum < $imincount));  # ichain is too short ?? is this good test? drop: $inum==0 for keepnointron
    $icvalid=-3 if($icvalid==1 && @did and scalar( grep{ $_ =~ m/$ichain/ } @did )); # ichain is subset of valids

		## NOPATH,nointron alts get dropped here: inum==0 (only 1 exon), icvalid == -1, all alts have diff trID ichain tho
		## eg: invalidchain    Funhe2Exx11m004617t1    z2      dropshort,0,    NOINTFunhe2Exx11m004617t1xe169
		
		## FIXME: non-alt diff loci here, check for common introns.  dropshort/-2 wrong unless shared introns
		## ? do this test for all? no, but icvalid==1 can use it, 
		## mark new locus if alt ichain large and completely diff from main ichain
		if(($icvalid == -2 || $icvalid == 1) and $inum>0 and not $ismain) {
			my %ins=(); map{ $ins{$_}++ } split /\W+/,$ichain; 
			my $icommon=0; map{ $icommon++ if($inset{$_}) } keys %ins;
			## add gid to inset to track old/new loci?
			#? my $icommon=0; map{ $icommon++ if($inset{$_} eq $gid) } keys %ins;
			## for oknewlocus, record new inset to get common new loci?
			
			$icvalid=3 unless($icommon>0); # icvalid=3/oknewlocus ??
		}
		
    my $debugval="";
use constant TEST_ALTEND => 1;
if(TEST_ALTEND) {
		# FIXME? get some dupl subchains -3? No, looks like 2 same-altends w/ diff ichains .. should be -1 did ichain
		if($icvalid == -3 and $gcds{$aid}) { 
			## gcds == geneends : cds ends or exon ends?
			# my($cb,$ce)=split",",$gcds{$aid}; ## do what now?

      my $validaltend=0;
			my($icb,$ice)=(split"; ",$ichain)[0,-1];  # front, back intron chain ids
			# test both icb,ice end points of ichain intron chain for exon longer than common size for internal exon.
			for my $ici ($icb,$ice) { 
			  next unless($exonends{$ici});
			  
			  my @ends= sort{$a<=>$b} keys %{$exonends{$ici}}; 
			  next unless(@ends >=2);
			  # exonends are counted, most freq = common inner exons << No, sum of intr=score gets common inner exons
			  # this max test fails unless larger ichain has several alts to count .. count not # alts, but intr=score
			  my($end1max,$end2max)= sort{ $exonends{$ici}{$b} <=> $exonends{$ici}{$a} or $a <=> $b } @ends;			  
			  ($end1max,$end2max)= ($end2max,$end1max) if($end1max > $end2max);
			  my($end1left,$end2right)= @ends[0,-1]; # sorted @ends
			  
        my($end1this,$end2this)= split ",", $gcds{$aid}{$ici}; # == $exonendsid{$aid}{$ici};
        $end1this||=$end1max; $end2this||=$end2max;
        if($end1this < $end1max) { 
          $validaltend++ if($end1this == $end1left); # exon left-end longer; keep? only if end1this == end1left ?
          } 
        if($end2this > $end2max) { 
           $validaltend++ if($end2this == $end2right);  #? exon right-end longer; only if end2this == end2right?
          }
    
        $debugval.="$ici:$end1this/$end1max/$end1left-$end2this/$end2max/$end2right,"
          if($debug);
      ## need some debug output here; why reject this one?
      # invalidchain    Funhe2Exx6m124896t6     -3,6    N72309,N72313;...
      ## but keep this:
      # intrchain       Funhe2Exx6m124896t12    2,4     N72313,N72314,N72315;...
          
			}
			$icvalid=2 if($validaltend); #?
	
      # test2 icvalid counts:		 # report these?
      # 50071 -1, 7576 -2, 14550 -3, 59776 1, 9642 2  << retained 10k as altend .. but right ones?

			## retain Only if cb,ce exon is longer than bigger chain exon; shorter is trash.
			## if($cb < exon-end[i] in @mainic or $ce > exon-end[j] in @mainic) {} have new alt end
			## then icvalid=1 ; modify ichain any as flag?
			##  maybe solution is add exon end locs to each intron chain string? should be all same for inner exons
			
			## append to intron chain as new end pts? or use simple flag for alt exon end?
			## need to know if exon end differs from larger intron chain exons
		}
}
		
    # galtsame opt: collect alt == main set, ichain eq dn1 //did{ichain}, ie same intron chain as maintr, for asmrna check of main    
    push @{$galtsame{$gid}},$aid if($icsameasmain); ## ($aid ne $gid and $ichain eq $dn1);
    
    # unless($did{$ichain} or $inum==0 or $inum < int(0.3*$i1) or scalar(grep{ $_ =~ m/$ichain/ } @did)) 
    ## add new alt anum to ATAB
    ##? output icvalid as %lclass: 'okay$icv', 'drop$icv' instead of number?
    my $valclass= $valclass{$icvalid} || $icvalid;
    $nclass{$icvalid}++; $ntabalt++ unless($ismain); # report
    if($icvalid>0)
    {  
      $did{$ichain}++; push @did, $ichain; 
      ## Maybe add valid ichain introns to inset{} ; maybe not: throws off above test
      #? if($icvalid == 3) { } # new gid?
      #  else { map{ $inset{$_}= $gid } split /\W+/,$ichain; }
      unless($ismain) { 
      	push @{$galtok{$gid}},$aid; $alt2gid{$aid} = $gid; $altnum{$aid}= ++$anum; 
      	
     		if($icvalid == 3) { $alt2gid{$aid} = "$gid;newlocus=1" ; } # oknewlocus handle how? 
      	## may be several alts of gid go to new but same locus, 
      	## what of newlocus2, newlocus3? leave that? or try to guess from intron equivs
      }
      print ATAB join("\t","intrchain",$aid,$TAG_ALT.$anum,"$valclass,$inum,$debugval",$ichain)."\n"; # this is long table, changed cols.
   	} elsif($debug) {  # debug or always?
   	  ++$invalnum; # track  invalid alt num?
      print ATAB join("\t","invalidchain",$aid,'z'.($invalnum+$anum),"$valclass,$inum,$debugval",$ichain)."\n"; # this is long table, changed cols.
   	}
   }
}

my $nclass= join", ", map{ $valclass{$_}.'='.$nclass{$_} } sort{ $b <=> $a } keys %nclass;
my $classrep= "# sum.intrchains classified $nclass for $ntabmain mains, $ntabalt alts as to $altab\n";
warn $classrep if $debug;
print ATAB $classrep;

=item ic -3 case

	#intron chains of alts for Funhe2Exx6m124896t1, inmain=10
	intrchain       Funhe2Exx6m124896t1     1,10,   N72304,N72305;...
	intrchain       Funhe2Exx6m124896t9     1,11,   N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t7     1,10,   N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t4     1,9,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t2     1,8,    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t14    -3,7,N72295,N72306,N72307,N72308:235875/235875/235551-236175/236175/236175,
		N72318:244848/244848/244848-245291/245291/245291, N72295,N72306,N72307,N72308;...
	intrchain       Funhe2Exx6m124896t3     1,9,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t8     1,7,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t5     1,8,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t15    1,6,    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t16    -1,6,   N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t17    -1,6,   N72295,N72306,N72307,N72308;...
	
	** gmap/gsplign mapping differs for t6, ichains differ & affect this test; above is all? gspl, later pub9set has gmap t6
	>> t6 end 237589 isn't longest, t7 is
	invalidchain    Funhe2Exx6m124896t6     -3,6,N72309,N72313:237589/237694/237549-237813/237813/237813,N72319:246562/24
	6562/245302-246834/246834/246834,       N72309,N72313;...
	
	invalidchain    Funhe2Exx6m124896t10    -3,5,N72309,N72313:237589/237694/237549-237813/237813/237813,N72319:246562/24
	6562/245302-246834/246834/246834,       N72309,N72313;...
	
	t7/Funhe2Eq7m077478t4/Funhe2Exx9m124883t16 longest exon at N72309,N72313: .. but doesnt end there. 
	grep 'Parent=Funhe2Exx6m124896t' kf2mixx_dedup.analt.gff | grep exon | grep 237549
	Scaffold10150   splkf2eg37mxx   exon    237549  237813  0.978   -       .       Parent=Funhe2Exx6m124896t7;trg=Funhe2Exx6m124896t7 1374 1641;splice=AGTC-
		;intr=+77/-46,N72309,N72313
	
	t6/Funhe2E6bm008467t4/Funhe2Exx9m124883t6
	grep 'Parent=Funhe2Exx6m124896t' kf2mixx_dedup.analt.gff | grep exon | grep 237589
	Scaffold10150   splkf2eg37mxx   exon    237589  237813  0.982   -       .       Parent=Funhe2Exx6m124896t6;trg=Funhe2Exx6m124896t6 1344 1567;splice=AGAG-
		;intr=+77/-46,N72309,N72313
	Scaffold10150   kfish2evg367    exon    237589  237813  98      -       .       Parent=Funhe2Exx6m124896t10;trg=Funhe2Exx6m124896t10 1224 1447;ix=6
		;intr=+77/-46,N72309,N72313

=item icvalid/invalid work

    ## ^ic -3 is commonest invalid, also where alt exons likely.
    ## kf2mixx count: 45270 -1(did) 7331 -2(short) 29196 -3(subset)
    ## eg: invalidchain    Funhe2Exx6m124896t6/Funhe2E6bm008467t4/Funhe2Exx9m124883t6     -3      N72309,N72313;...
		## may well be case of valid alt 3' exon end, longer than spliced exon at that spot of larger alts.
		## t1,t2,t4 have same exon w/ splice in middle..

	> t6 3'end exon, 105b longer exon, cds 69b longer
	Scaffold10150   splkf2eg37mxx  exon    237589  237813  0.982   -       .       Parent=Funhe2Exx6m124896t6;trg=Funhe2Exx6m124896t6 1344 1567;splice=AGAG-;
		intr=+77/-46,N72309,N72313 << this is same inID set as t2, but -46 means one N72309? maps inside exon ..
	Scaffold10150   splkf2eg37mxx  CDS     237625  237813  1       -       .       Parent=Funhe2Exx6m124896t6
	> t2 same exon
	Scaffold10150   kf2evg367mixy  exon    237694  237813  100     -       .       Parent=Funhe2Exx6m124896t2;trg=Funhe2Exx6m124896t2 1351 1470;ix=8;
		intr=123,N72309,N72313
	Scaffold10150   kf2evg367mixy  CDS     237694  237813  100     -       .       Parent=Funhe2Exx6m124896t2
	Scaffold10150   kf2evg367mixy  exon    237420  237539  100     -       .       Parent=Funhe2Exx6m124896t2;trg=Funhe2Exx6m124896t2 1471 1590;ix=9;intr=386,N72308,N72309,N72310,N72311,N72312
	
	grep Funhe2Exx6m124896t  kf2mixx_dedup.altbest.tab | sed 's/;.*/;.../;' | less
	#intron chains of alts for Funhe2Exx6m124896t1
	intrchain       Funhe2Exx6m124896t1     N72304,N72305;...
			^^ exon 234272  234420 ; intr=+113/-3,N72304,N72305 3'end exon;..
	intrchain       Funhe2Exx6m124896t9     N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t7     N72303,N72304,N72305;...
			^^ exon 234128  234420 ; intr=+113/-38,N72303,N72304,N72305 : 150b longer 3'end, hits new intron
	intrchain       Funhe2Exx6m124896t4     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t2     N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t14    -3      N72295,N72306,N72307,N72308;...
	intrchain       Funhe2Exx6m124896t3     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t8     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t5     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t15    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t16    -1      N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t17    -1      N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t6     -3      N72309,N72313;...							<<* likely valid altend
	invalidchain    Funhe2Exx6m124896t10    -3      N72309,N72313;...
	invalidchain    Funhe2Exx6m124896t12    -3      N72313,N72314,N72315;...
	altsof  Funhe2Exx6m124896t1     Funhe2Exx6m124896t9,Funhe2Exx6m124896t7,Funhe2Exx6m124896t4,Funhe2Exx6m124896t2,Funhe2
	Exx6m124896t3,Funhe2Exx6m124896t8,Funhe2Exx6m124896t5,Funhe2Exx6m124896t15
	
	# start-5': N72318,N72319; N72319 (but for shorter 5' start N72317,N72318; N72318 inside longer : is it valid? see t15)
	intrchain       Funhe2Exx6m124896t1     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t9     N72303,N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318; N72319; N72319
	intrchain       Funhe2Exx6m124896t7     N72303,N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t4     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t2     N72295,N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t14    -3      N72295,N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	intrchain       Funhe2Exx6m124896t3     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t8     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t5     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t15    N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t16    -1      N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t17    -1      N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t6     -3      N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t10    -3      N72309,N72313; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t12    -3      N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319

=cut
		

print ATAB "#valid-alts for main=$maingenes, alts=$altgenes\n";
foreach my $gid (sort keys %galtok) { 
  my @v= @{$galtok{$gid}};
  print ATAB "altsof\t$gid\t",join(",",@v),"\n";
}  

# upd main equiv
print ATAB "#altsame as main=$maingenes, alts=$altgenes\n";
foreach my $gid (sort keys %galtsame) { 
  my @v= @{$galtsame{$gid}};
  print ATAB "altsame\t$gid\t",join(",",@v),"\n";
}  

close(ATAB);


# also create altbest.gff from these.. changing alt ID to main.tN...
# ADD option to NOT change alt id.
# ADD option not to write this valtgff .. make sure ATAB is sufficient to regenerate alt.gff
# FIXME: oknewlocus class needs to be flagged, ID changed to non-alt, per CHANGEALTID option ..
#   .. If this is accurate reclass, which seems to be mostly.
# kf2pub11_both.altbest.tab: valid-alts for main=kf2pub11_both.anmain.gff from alts=kf2pub11_both.analt.gff
# sum.intrchains classified oknewlocus=6787, okaltend=15150, okay=63817, dropdup=63739, dropshort=22672, dropsubchain=26996 for 25364 mains, 173797 alts as to kf2pub11_both.altbest.tab
## newloc eg
#intronchains for mainID=Funhe2Exx11m000009t1, mainIn=94
# intrchain       Funhe2Exx11m000009t1    t1      okay,94,        N158164;..
# intrchain       Funhe2Exx11m000009t2    t2      okay,93,        N158164;..
# intrchain       Funhe2Exx11m000009t3    t3      okay,93,        N158164;..
# ..
# intrchain       Funhe2Exx11m000009t26   t16     oknewlocus,45,  N99571,N99577;..
# .. N99571,N99577; N99578,N99579,N99580; N99580,N99583,N99584,N99586,N99587; N99586,N99588; N99587,N99588,N99589; N99589,N99590; N99590; N99592; N99592,N99593; N99593,N99594; N99594,N99595; N99595,N99596; N99596,N99598; N99598,N99599; N99599,N99600; N99600,N99601; N99601,N99602; N99602,N99603; N99603,N99606; N99606,N99607; N99607; N99608,N99609; N99608,N99609,N99610,N99611; N99611,N99612,N99618; N99618,N99619; N99616,N99620; N99617,N99620,N99621,N99622; N99622,N99623; N99623,N99624; N99624,N99625,N99626; N99626,N99629; N99629,N99630; N99630; N99631; N99631,N99632; N99632,N99633; N99633,N99634; N99634,N99635; N99635,N99637; N99637,N99638; N99638,N99639; N99639,N99641; N99641,N99642; N99642,N99643; N99643,N99644; N99644
# invalidchain    Funhe2Exx11m000009t29   z24     dropshort,34,   N158241;..
# intrchain       Funhe2Exx11m000009t25   t17     oknewlocus,39,  N99580,N99583,N99
# intrchain       Funhe2Exx11m000009t28   t18     oknewlocus,39,  N99580,N99583,N99
# invalidchain    Funhe2Exx11m000009t27   z27     dropshort,30,   N158169,N158170,N
# intrchain       Funhe2Exx11m000009t19   t19     oknewlocus,23,  N99631;..
# intrchain       Funhe2Exx11m000009t18   t20     oknewlocus,21,  N99631;..
# invalidchain    Funhe2Exx11m000009t32   z32     dropshort,13,   N158259;..
# intrchain       Funhe2Exx11m000009t30   t21     oknewlocus,10,  N99646;.. << dropshort-newlocus
# intrchain       Funhe2Exx11m000009t35   t22     oknewlocus,1,   N99658;.. << dropshort-newlocus


my $gf= $altgenes;
my $gop = ($gf =~ /.gz/) ? "gunzip -c $gf" : "cat $gf"; 
open(G,"$gop |") or die "bad $gop";

warn "# $valtgff: valid-alts for main=$maingenes from alts=$altgenes\n" if $debug;
open(AGFF,">$valtgff") or die "$valtgff"; 
print AGFF "##gff-version 3\n";
print AGFF "# valid-alts for main=$maingenes from alts=$altgenes\n";
my ($altid,$isalt,$isnewmain,$galt,$gmain);
while(<G>) {
  next unless(/^\w/);
  if(/\tmRNA/) { 
    my($tid)=m/ID=([^;\s]+)/;
    $gmain= $alt2gid{$tid} || "";
    $isalt= ($gmain)? 1: 0;
    
    $isnewmain= ($gmain =~ /newlocus/)?1:0;
    if($isnewmain) {
      $galt=$tid; ($altid=$gmain) =~ s/;.*//; $altid.="new"; # what for new loc? may be several alts > new
      my $ai= 1;  
      if($altid =~ /${TAG_ALT}1$/) { $altid =~ s/${TAG_ALT}1$//; }  
      $altid .= $TAG_ALT.$ai;
      if($CHANGEALTID) { s/ID=$tid/ID=$altid;newlocus=1;oldaltid=$galt/; } 
      else { s/$/;/ unless(m/;/); s/;/;newlocus=$altid;/; }
      
    } elsif($isalt) { 
      $galt=$tid; $altid=$gmain; # change t1 to taltnum, but if gmain == t2..tn then problem.
      my $ai= $altnum{$tid} || 99;  
      if($altid =~ /${TAG_ALT}1$/) { $altid =~ s/${TAG_ALT}1$//; } # else problem or not?
      $altid.= $TAG_ALT.$ai;
      
      ## FIXME: oid=OtherOldId already exists; change to oid2=? origid=? lastid=?
      if($CHANGEALTID) { s/ID=$tid/ID=$altid;oldaltid=$galt/; } 
      else { s/$/;/ unless(m/;/); s/;/;newaltid=$altid;/; }
      
      unless(m/;alttr=1/) { s/$/;alttr=1/; }
      }  
  } elsif(/\t/ and $CHANGEALTID and ($isalt||$isnewmain)) {
    s/Parent=$galt/Parent=$altid/;
  } 
  print AGFF $_ if ($isalt||$isnewmain);
}
close(G); close(AGFF);

__END__

=item intronchains parsing

/bio/bio-grid/kfish2/rnas/kf2evgr/trevg367mixx/alt11
nam=kf2pub11_both

$evigene/scripts/overlapfilter \
 -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum -mark intr \
 -over $kfish2/intron/intron_okay.gff.gz -in $nam.gff.gz > $nam.an2.gff

 wc -l $nam.an2.*chains
   33075 kf2pub11_both.an2.galtchains
   11772 kf2pub11_both.an2.gsamechains
  183032 kf2pub11_both.an2.ichains

## intron chains table:

grep intr= $nam.an2.gff | perl -ne \
'($d)=m/Parent=(\w+)/; ($in)=m/intr=([^;\s]+)/;
($iv,$in1,$in2)=split",", $in; $iv=~s,/.*,,; if($iv>0) { map{
$gin{$d}{$_}++; $ic{$_}{$d}++ } grep/\w/, ($in1,$in2); } END{ for $g
(sort keys %gin) { @in=sort keys %{$gin{$g}}; %lg=(); $ni=@in; $niv=0;
for $i (@in) { $niv+=$gin{$g}{$i}; map{ $lg{$_}++; }  keys %{$ic{$i}};
} ($gd=$g)=~s/t\d+$//; $minc=$ni*0.20; @lg= map{ my $c=$lg{$_};
$c.="NO" if($c<$minc);  s/$gd/id_/; "$_:$c"; } grep{ $_ ne $g }
sort{$lg{$b}<=>$lg{$a} or $a cmp $b} keys %lg;
$ng=join",",grep/NO/,@lg; $lg=join",",grep{ not /NO/} @lg; $ng||="na";
$lg||="na"; print join("\t",$g,$ni,$niv,$lg,$ng)."\n";  }}  ' \
  > $nam.an2.ichains

cat $nam.an2.ichains | cut -f1,2,4 | perl -ne \
'($d,$ni,$ad)=split; ($g,$t)= $d=~m/(\w+)t(\d+)$/; @at=map{ s/id_t//;
s/:\w+//;  $_; } grep /id_/,split",",$ad; @ot=sort{$a<=>$b}($t,@at);
$ot=$ot[0]; map{ $sg{$g}{$ot}{$_}++ } @ot; if($lg and $lg ne $g) {
@t=sort keys %{$sg{$lg}}; %did=(); for $t (@t) { @at=sort{$a<=>$b}
keys %{$sg{$lg}{$t}}; $new=1; map{$new=0 if($did{$_}++); }@at;
if($new) { $at=join",",@at; print join("\t",$lg,"t".$t,$at)."\n"; } }
} $lg=$g; ' \
 > $nam.an2.galtchains

cat $nam.an2.ichains | cut -f1,2,4 | perl -ne \
'($d,$ni,$ad)=split; next if($ad eq "na"); ($g,$t)=
$d=~m/(\w+)t(\d+)$/; @at=map{  s/:\w+//;  $_; } grep { not /id_/ }
split",",$ad; next unless(@at); @ot=sort{$a<=>$b}(@at); $ot=$t; map{
$sg{$g}{$ot}{$_}++ } @ot; if($lg and $lg ne $g) { @t=sort keys
%{$sg{$lg}}; %did=(); for $t (@t) { @at=sort keys %{$sg{$lg}{$t}};
$new=1; map{$new=0 if($did{$_}++); }@at; if($new) { $at=join",",@at;
print join("\t",$lg,"t".$t,$at)."\n"; } } } $lg=$g; ' \
> $nam.an2.gsamechains

==> kf2pub11_both.an2.galtchains <==
  n=4855 loci (20%) have diff locus alts (>t1), ng=29135 loci have some altchain, ng1=26624 loci have t1 chain (vs 27492 below)
Funhe2Exx11m000001      t1      1
Funhe2Exx11m000002      t1      1
Funhe2Exx11m000003      t1      1
Funhe2Exx11m000003      t3      3
Funhe2Exx11m000004      t1      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
..
Funhe2Exx11m000012      t1      1,2,3,4,5,9,10,13,14,15,16,18  << Scaffold10067
Funhe2Exx11m000012      t17     17   << Scaffold9902
Funhe2Exx11m000012      t6      6,7,8,11,12   << Scaffold9981
  
==> kf2pub11_both.an2.gsamechains <==
  n=9988 loci (33%) have other at same locus, ng1=7312 t1 have same-locus
Funhe2Exx11m000009      t18     Funhe2Exx11m000013t1,Funhe2Exx11m000013t2,Funhe2Exx11m000013t3,Funhe2Exx11m000013t4,Funhe2Exx11m000013t7
   ^^ Scaffold14 09t18,013t (other 09t1.. at Scaffold323)
Funhe2Exx11m000012      t11     Funhe2Exx11m000093t1,Funhe2Exx11m000093t2,Funhe2Exx11m000093t3
Funhe2Exx11m000012      t17     Funhe2Exx11m000027t1,Funhe2Exx11m000027t10,Funhe2Exx11m000027t11,Funhe2Exx11m000027t12,Funhe2Exx11m000027t13,Funhe2Exx11m000027t14,Funhe2Exx11m000027t15,Funhe2Exx11m000027t16,Funhe2Exx11m000027t17,Funhe2Exx11m000027t18,Funhe2Exx11m000027t19,Funhe2Exx11m000027t2,Funhe2Exx11m000027t22,Funhe2Exx11m000027t25,Funhe2Exx11m000027t3,Funhe2Exx11m000027t4,Funhe2Exx11m000027t5,Funhe2Exx11m000027t6,Funhe2Exx11m000027t7,Funhe2Exx11m000027t8,Funhe2Exx11m000027t9
Funhe2Exx11m000013      t1      Funhe2Exx11m000009t18,Funhe2Exx11m000009t19,Funhe2Exx11m000009t25,Funhe2Exx11m000009t26,Funhe2Exx11m000009t28
	gsamechains t1 vs kf2pub11_both.main.eqgene:
	  of samechain t1=7312, noeqgene n=3631; eqcds>=20% n=2190; eqgene > 0 n=4410;
		of main.eqgene w/o samechain, n=74894, nno=68374 no eqgene, 20%eqcds n=2451
	main.eqgene total n=79899, nno=69957,  20%eqcds n=4641, 
		
 >> main.eqgene w/ no samechain, 20%eqcds n=2451 should be marked same locus likely missing valid introns (but check some?)
 >> samechain w/ no main.eqgene, noeqgene n=3631: what effects? : many are t1 x t2+, not tested in main.eqgene
   * check qual of shared intronchains for same-loci: weak vs strong? some are utr overlap junk
   * probably should use only equalgene cds overlap to call same-locus, otherwise utr overlaps pile up.
   * redo equalgene for all main+alts.
   
egrep '(Funhe2Exx11m000052|Funhe2Exx11m040951)' kf2pub11_both.an2.ichains | cut -f1,2,4 
Funhe2Exx11m000052t1    14      id_t2:10,id_t3:10,id_t5:10,id_t6:10,id_t8:9,id_t9:9,id_t4:8,id_t7:8,
		> utrover? Funhe2Exx11m040951t1:3,Funhe2Exx11m040951t2:3,Funhe2Exx11m040951t3:3,Funhe2Exx11m040951t4:3,Funhe2Exx11m040951t5:3,Funhe2Exx11m040951t6:3
		> no main/alt.eqgene cds or utr overlap for g052 x g40951
Funhe2Exx11m040951t1    4       id_t2:4,id_t3:4,id_t4:4,id_t5:4,Funhe2Exx11m000052t1:3,id_t6:3
	> share UTR exons, no CDS, w/ g52t1: Scaffold10165:612509-612543,642313-642387

==> kf2pub11_both.an2.ichains <==
	nt=183032 tr have ichains, and ng=29136 evasm loci, of nt=322715,ng=92004 in kf2pub11 mRNA; ng1=27492 loci with t1 in ichain
Funhe2Exx11m000001t1    65      127     na      Funhe2Exx11m007299t12:1NO,Funhe2Exx11m035699t12:1NO,Funhe2Exx11m035699t7:1NO
Funhe2Exx11m000002t1    66      131     na      Funhe2Exx11m000007t1:3NO,Funhe2Exx11m000007t2:3NO,Funhe2Exx11m000007t3:3NO,Funhe2Exx11m000007t4:1NO
Funhe2Exx11m000003t1    19      35      na      na
Funhe2Exx11m000003t3    4       6       na      na
Funhe2Exx11m000004t1    131     255     id_t2:129,id_t3:126,id_t4:126,id_t5:124,id_t6:124,id_t8:124,id_t7:123,id_t9:122,id_t10:121,id_t11:117,id_t12:113,id_t13:106,id_t14:99,id_t15:96,id_t18:94,id_t16:91,id_t17:85,id_t25:83,id_t19:82,id_t20:82,id_t21:79,id_t22:77,id_t23:77,id_t27:76,id_t24:72,id_t26:71,id_t28:67       id_t29:12NO,id_t30:12NO,id_t31:8NO,Funhe2Exx11m028732t1:8NO,Funhe2Exx11m052498t1:5NO,Funhe2Exx11m067108t1:4NO,Funhe2Exx11m049230t1:1NO
Funhe2Exx11m000004t10   121     239     id_t1:121,id_t2:121,id_t3:121,id_t4:121,id_t5:121,id_t6:121,id_t8:121,id_t7:120,id_t9:119,id_t11:117,id_t12:111,id_t13:103,id_t14:96,id_t18:89,id_t15:88,id_t16:88,id_t19:79,id_t17:77,id_t23:77,id_t25:76,id_t22:75,id_t20:74,id_t21:74,id_t24:70,id_t26:68,id_t27:67,id_t28:64        id_t29:12NO,id_t30:12NO,id_t31:8NO,Funhe2Exx11m028732t1:8NO,Funhe2Exx11m052498t1:5NO,Funhe2Exx11m049230t1:1NO,Funhe2Exx11m067108t1:1NO


=cut
