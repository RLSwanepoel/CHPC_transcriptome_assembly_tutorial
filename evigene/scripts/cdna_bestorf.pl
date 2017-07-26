#!/usr/bin/perl
# cdna_bestorf.pl cut from genefindcds.pl

=item about
  
  finds best orf, adds proteins and CDS to transcript assembly cdna
    cdna_bestorf.pl  -cdna transcript.fasta  > transcript.aa

=item todo
   
   FIXMEd: read all tr 1st is bad idea, huge.tr set chokes
   fixed: new option: noutrorf=turn off utrorf output, see splitutrorfs

	 ** add hdr report of longest revaa (2nd longest?) unless fwd/rev spec.
	 -- some (rare?) cases of shorter revaa is true prot; homology &/or gmap intron splicing validate this
	 
   clean up "action" options; separate opts
   add frameshift detect option (cdna_proteins.pm)
   
=item author
  
  don gilbert, gilbertd near indiana edu, 2011
  part of EvidentialGene, evigene/scripts/
  cut from genefindcds.pl, which is drawn largely from Brian Haas's PASA scripts

=cut

use FindBin;
use lib ("$FindBin::Bin");

use strict;
use warnings;
use Getopt::Long;
use cdna_proteins;

use constant VERSION  => '20130818';  # 0228'20120721'; 

my ($cdnaseq, $nin, $skipsubsetseq, $splitutrorfs, $noutrorf, $revaa_report)= (0) x 20;

my $action="fasta";  # or cdnafasta  cdsfasta ..
my $debugin= undef;
my %cdnaseq; my %cdnahead;

# add -output option to file updated genes
my($oformat,$ostrand,$opart,)=("") x 10;
my($aaout,$cdsout,$cdnaout)= (undef) x 3;

my $optok= GetOptions(
  "cdnaseq=s", \$cdnaseq,  
  "aaseq|outaa|output:s", \$aaout, #? make default -out=infile.cds unless -out stdout
  "cdsseq|outcds:s", \$cdsout, #? make default -out=infile.cds unless -out stdout
  "outcdna:s", \$cdnaout, #? make default -out=infile.cds unless -out stdout
  "action=s", \$action,    # need to specify actions: 
      # now packing on:  "fwd|rev" strand + cdna|cds|noaa,fasta|notfa + full|longpart|bestpart
   "oformat=s", \$oformat,    # replace action: 
   "ostrand=s", \$ostrand,    # replace action: 
   "opart=s", \$opart,    # replace action: 
     
  "minaa=i", \$MINAA,  # not used ?
  "fullorf:s", \$ORF_FULLvPART,  
  "nostopcodon", \$NoStopCodon, 
##  "skipsubsetcdna!", \$skipsubsetseq, 
  "splitutrorfs!", \$splitutrorfs, 
  "noutrorf", \$noutrorf, 
  "revaareport", \$revaa_report, # dont need special opt? debug/verbose/action or default?
  "goodlen!", \$USEGOODLEN,  # add MINGOOD
  "Selenocysteine|selc!", \$USESelenocysteine,  
  "goodmin=s", \$MINGOOD,  
  "debug:i", \$debugin, 
  );

die "usage: cdna_bestorf.pl -cdna cdna.fa > protein.aa
 opts: -minaa=$MINAA  -nostopcodon -goodmin=$MINGOOD -splitutrorf   
      -action=fasta|table, [cdna|cds][only]fasta, [fwd|rev|both]strand, [full|long|best]part
" unless($optok and  $cdnaseq);

# push @input, "stdin" unless(@input);
if(defined $debugin) { $DEBUG=($debugin>0)?$debugin:1; } else { $DEBUG=0; }

if($action =~ /reportrevaa|reportrev|report/) { # messy place to add this; revfasta conflict
	$revaa_report=1; $action =~ s/reportrevaa|reportrev|report//;
}
$action=$ostrand."strand$action" if($ostrand);
$action=$opart."part$action" if($opart);
$action=$oformat.$action if($oformat);

# "action=s", # now packing on:  "fwd|rev" strand + cdna|cds|noaa,fasta|notfa + full|longpart|bestpart
# .. split action into format=fasta|table + strand=fwd|rev|both + part=full|long|best
# .. -act cdnanoaafasta  or -act cdnaonlyfasta
#  my $whichstrand= ($action =~ /rev/)? "rev" : ($action =~ /fwd/) ? "fwd" : "both";
#  my $fullpart= ($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
#  my $outtype = ($action =~ /fasta/)?"fasta":"table";
my $noprotaction= ($action=~/cdna|cds/ and $action=~m/noaa|only/)?1:0;  # and $outtype eq "fasta"

# FIXed: add option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
## protein=MLSHQLLEDSTMMQMKHGLRQGRENICQGSRLLLIGNVLVDNXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
$MINGOOD= $MINGOOD/100.0 if($MINGOOD > 1); 

##was $ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART > 1); 
 # ^^ -full=0 means only full orfs?; -full=1 means never?? change to -full=1 means full only; -full=0 longest part always
$ORF_FULLvPART=1 unless($ORF_FULLvPART =~ /\d/);
$ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART >= 1); 
$ORF_FULLvPART=1 if($ORF_FULLvPART == 0); # means always longest

useSelenocysteine() if($USESelenocysteine);

sub openout {
  my($oname,$osuf,$ohandle)= @_;
  my $OUTH= undef;
  if(defined $oname and not $oname) { # use cdnaseq name
    $oname= ($cdnaseq =~ /stdin|^\-$/)? "cdna_bestorf_out" : $cdnaseq;
    $oname=~s/\.gz//; $oname=~ s/\.\w+$//; $oname.= $osuf;
  }
  if($oname) { 
    $oname.= "out$osuf" if($oname eq $cdnaseq);
    open($OUTH, ">$oname") or die "ERR writing $oname"; $ohandle= $OUTH; # *OUTH; 
  }
  return($ohandle,$oname);
}

# outputs:
## drop use of "-action=xxxcds -action=xxxcdna for outputs
$cdsout="" if($action =~ /cds/i and not defined $cdsout); # define output
$cdnaout="" if($action =~ /cdna/i and not defined $cdnaout); # define output

my($OUTAA,$onameaa)= openout($aaout,".aa",*STDOUT);
my($OUTCDS,$onamecds)= openout($cdsout,".cds",undef);
my($OUTCDNA,$onamecdna)= openout($cdnaout,".cdna",undef);


MAIN: { 
  my $ovh; my $ok= 0; 
  my($nin,$tnin,$tngood,$tnsplit,$tnskip)=(0) x 9;
  if($cdnaseq =~ /.gz$/) { $ok= open(OVR,"gunzip -c $cdnaseq |");  $ovh= *OVR; }
  elsif($cdnaseq =~ /stdin|^\-$/) { $ok=1; $ovh= *STDIN; }
  elsif($cdnaseq) { $ok= open(OVR,$cdnaseq); $ovh= *OVR; }
  die "bad -cdnaseq=$cdnaseq" unless($ok);
  my $id="none";
  while(<$ovh>) {
    if(/^>(\S+)/) { $id=$1; 
      do{ my($ngo,$nspl,$nsk)= cdna_bestorf(); $tngood+=$ngo;$tnsplit+=$nspl;$tnskip+=$nsk;
         %cdnaseq= %cdnahead=(); $nin=0; } if($nin>49);
      $cdnaseq{$id}="";  $nin++;  $tnin++;
      if(m/ +(\S.+)$/){ 
      	my $h=$1; $h =~ s/\s*\b(len|cf|nt)=\S+//g; 
      	$h =~ s/\baa c offs=/offs=/;	## fixup old mess: 'aa c offs=xxx' from bad (aa|c)len= cut
      	$cdnahead{$id}=$h; }
      ## FIXME: cut (len|cf|nt)= is bad? need option/spec of atts to drop, esp. replace new atts
    } elsif(/\w/) { chomp; $cdnaseq{$id}.= uc($_); } } 
  close($ovh);

  my($ngo,$nspl,$nsk)= cdna_bestorf(); 
  $tngood+=$ngo;$tnsplit+=$nspl;$tnskip+=$nsk;
  warn"#cdna_bestorf from cdnaseq:$cdnaseq n=$tnin, ngood=$tngood , nsplit=$tnsplit, nskip=$tnskip\n" 
     if $DEBUG;
  $onameaa||=""; $onamecds||=""; $onamecdna||=""; #nogood# $^W=0; # $WARNING=0;
  warn"#cdna_bestorf output to aa:$onameaa  cds:$onamecds cdna:$onamecdna \n" if $DEBUG;

}  

#..................

# ## this is slow.. not best way to do this. cd-hit-454 -c 1.00 will be faster
# my %containedin;
# sub cdna_containedin  
# {
#   my @cids= sort { length($cdnaseq{$b}) <=> length($cdnaseq{$a}) or $a cmp $b }  keys %cdnaseq; 
#   my $allseq=";"; # may be large; likewise %cdnaseq
#   foreach my $cid (@cids) {
#     my $seq = $cdnaseq{$cid};
#     my $at= index($allseq, $seq); # need to do revcomp?
#     if($at < 0) {
#       my $rseq= revcomp($seq);
#       $at= index($allseq, $rseq);  
#       $seq=$rseq if($at >= 0);
#     }
#     
#     # while($at >= 0) .. get all superseq ids
#     if($at >= 0) {  
#       my $c= rindex($allseq,":",$at);  my $b= rindex($allseq,";",$c);
#       my $inid= substr($allseq,$b+1,$c-$b-1);
#       $containedin{$cid}= $inid; ## need all superseq ids ?
#       }
#     $allseq.="$cid:$seq;";
#     }
#   return @cids;
# }


sub cdna_bestorf
{
  # FIXmaybe: add other option before output: skip subset cdna contained in *good* longer cdna
  # .. change sorting of @cdnain to long > short-containedin
  my @cdnaid; 
     @cdnaid= sort keys %cdnaseq; 
  
  my %isgoodseq;
  my $ncdna= @cdnaid;
  # warn"#cdna_bestorf from cdnaseq:$cdnaseq n=$ncdna\n" if $DEBUG;
  my @outadd=(); push @outadd, "cds" if(defined $OUTCDS); push @outadd, "cdna" if(defined $OUTCDNA);

  my ($ngood,$nsplit,$nskip)= (0) x 10; 
  my $MINTR= $MINAA * 3;
  my @rev=($action =~ /rev/)? (1) : ($action =~ /fwd/)? (0) : (0,1);
  my $whichstrand= ($action =~ /rev/)? "rev" : ($action =~ /fwd/) ? "fwd" : "both";
  my $fullpart= ($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
  $fullpart= $whichstrand."strand".$fullpart;
  
# ADD maybe: detect joined/fused genes using 
#  1. cds span < ~ 1/2 trspan, leaving long utr, << annotate these cases, require aa complete/partial5
#  2. check the utr for long cds
 
  foreach my $id (@cdnaid) {
    my $cdnain= $cdnaseq{$id}; 
    my $cdnasize= length($cdnain); 
    my $cdnabest= $cdnain;
    my($orfprot, $prostart5, $proend3, $bestorf, $utrorf, $orflen, $strand, $isrev, $aalen, $pcds, $compl)= (0) x 20;
    my $doskip= 0;
    my($fahead,$orfprotfa,$revaalen,$revinfo)=("") x 9;
    
    if($cdnasize < $MINTR) { 
      $doskip=1; $nskip++; # warn??
    } else {

    #...... new : both dirs ...

    # FIXed: respect action =~ rev or fwd here: fullpart fixed
    my $allorfs;
    ( $bestorf, $allorfs)= getBestProt2( $fullpart, $cdnain); # , undef, $oldStart_b,$oldStart_e
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
    ( $orfprot, $prostart5, $proend3, $orflen, $strand)= orfParts($bestorf, [qw(protein start stop length orient)]);
    $isrev= ($strand eq '-') ? 1 : 0;
    #NO, always input strand now# $cdnabest= revcomp($cdnain) if($isrev); #?? maybe not, only for output here.
    
    ($utrorf)= utrorf_test( $bestorf, $allorfs, $cdnasize); # noutrorf option here or below?
    
		if($revaa_report) {
			($revaalen,$revinfo)= revorf_report( $bestorf, $allorfs, $cdnasize); # return info only for longest reversed orf
			# revinfo == "aarev=55%[rev2fwd],aa99,90%,complete,s-,o9-309"
		}

#    ## FIXME: check on 200 of 3000 that have good pCDS : 70..80, but utr long enough to pick amino..
# >veln5ptf014k33Loc11t2cdna aalen=1067,76%,complete; cutlen=3301,at=3301; clen=4164; strand=+; offs=54-3257; 
# >veln5ptf014k29Loc21t5cdna aalen=1480,76%,complete; cutlen=4486,at=4486; clen=5798; strand=+; offs=36-4478; 
# >veln5ptf013k41Loc3t1cdna aalen=1024,63%,complete; cutlen=3457,at=1408; clen=4865; strand=-; offs=4552-1478; 
# >veln5ptf009k29Loc3t1cdna aalen=1710,81%,complete; cutlen=5445,at=5445; clen=6322; strand=-; offs=5418-286; 
# >veln5ptf004k41Loc2t3cdna aalen=2209,61%,complete; cutlen=9712,at=1006; clen=10718; strand=+; offs=1062-7691; 

###    #.........
#### old vers 
#     for my $rev (@rev) {
#       if($rev==1){  $cdnain= revcomp($cdnain); } 
#       my($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= getBestProt($fullpart, $cdnain); # , undef, $oldStart_b,$oldStart_e
# 
#       my $change=0;
#       ($change,$orfprot,$bestorf) = bestorf_test($bestorf,$bestorfi);
#       if($change) { 
#         $isrev= $rev; $cdnabest= $cdnain; $utrorf= $utrorfi;
#         $strand=($isrev==1)?"-":"+";  # crev
#         ( undef, $prostart5, $proend3) = orfParts($bestorf);
#         }
#       }
#     if($utrorf) {
#       my($orfok) = bestorf_test(undef,$utrorf);
#       $utrorf=0 unless($orfok);
#     }
#  # old
      
    } # MINTR
 
     
    #above# my($fahead,$orfprotfa,$revaalen,$revinfo)=("") x 9;
    if($orfprot) {
      ($aalen,$pcds,$compl,$orflen,$fahead,$orfprotfa)
          = proteindoc($bestorf,$cdnasize,$strand,($action =~ /fasta/));

		  ## note: fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;";
      #  if( $skipsubsetseq ) {}   #.. drop this, use other methods to remove subset seq
      $isgoodseq{$id}++ unless($doskip);
      $ngood++ unless($doskip); # unless($utrorf) ??

			$fahead.=" $revinfo;" if($revaa_report and $revaalen >= $MINAA and $revinfo);
				# revinfo == "aarev=55%[rev2fwd],aa99,90%,complete,s-,o9-309"
        ## revaa_report : do we want option to also print OUTAA >id.rev\nrevprotfa ??
        ## OR option to replace fwdaa with revaa output (aa and cds), IF %revaa/fwdaa >= minrevaa
    }
    
    if($doskip) {
      # warn? any output?  $nskip++ above
      
    } elsif($action =~ /fasta/) {  # FIXME, do fasta or table here, format at print...
      my $cdnah= $cdnahead{$id}||"";
			$cdnah =~ s/\s*\b(aalen|clen|strand|offs)=[^;\s]+[;]?//g;
			## above fixup old mess: 'aa c offs=xxx'
			# 2013feb: add output files: aaout, cdsout, cdnaout?
      
      if($action =~ /all/ or $orfprot) {
      
        unless($noprotaction) {
        print $OUTAA ">$id $fahead $cdnah\n";
        print $OUTAA $orfprotfa,"\n";
        }
        ## revaa_report : do we want option to also print OUTAA >id.rev\nrevprotfa ??
        ## OR option to replace fwdaa with revaa output (aa and cds), IF %revaa/fwdaa >= minrevaa
        
        my $utrorfseq=""; my $utrcut=0;
        # FIXME for utrorf: find way to split cdna b/n 1st, 2nd orf: midway?
        if($utrorf) {
            # main:  $orfprot, $prostart5, $proend3, $orflen, $strand
          my( $uprostart5, $uproend3, $uorflen, $ustrand)= orfParts($utrorf, [qw( start stop length orient)]);
          my($aalen,$pcds,$compl,$orflen,$fahead,$orfprotfa)
            = proteindoc($utrorf,$cdnasize,$ustrand, 1);
          unless($noprotaction or $noutrorf) {
          print $OUTAA ">",$id,"utrorf $fahead $cdnah\n";
          print $OUTAA $orfprotfa,"\n";
          }
          
          if($splitutrorfs) {
            # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
            # *** matters here if cdnabest= revcomp(cdnain) ***
            my($b1,$e1)= ($prostart5 > $proend3) ? ($proend3,$prostart5) : ($prostart5,$proend3);
            my($b2,$e2)= ($uprostart5 > $uproend3) ? ($uproend3,$uprostart5) : ($uprostart5,$uproend3);
            my $cut= ($b2 >= $e1) ? $e1 + int(($b2 - $e1)/2) 
                   : ($e2 <= $b1) ? $e2 + int(($b1 - $e2)/2) : 0;
            if($cut>0) {
              $utrcut  = $cut;
              my $seq1 = substr($cdnabest,0,$cut);
              my $seq2 = substr($cdnabest,$cut);
              $cdnabest=  ($prostart5 > $uprostart5) ? $seq2 : $seq1;
              $utrorfseq= ($prostart5 > $uprostart5) ? $seq1 : $seq2;
              # $cdnasize=length($cdnabest); #?? change or not; need change prostart/end also

## was BAD cut: rev, cant get 651aa from 857cds << reverse cseq, orfseq
# >veln5ptf001k29Loc11t1 aalen=651,64%,complete; clen=3019; strand=-; offs=2949-994; 
# >veln5ptf001k29Loc11t1.utrorf aalen=239,23%,partial5-utrbad; clen=3019; strand=+; offs=2-721; 
# >veln5ptf001k29Loc11t1.cdna aalen=651,64%,complete; cutlen=857,at=857; clen=3019; strand=-; offs=2949-994; 
# >veln5ptf001k29Loc11t1.utrorf.cdna aalen=239,23%,partial5-utrbad;cutlen=2162,at=857; clen=3019; strand=-; offs=2-721; 

              }
            }
          }
  

# FIXME: remove tag from ID, add as type=tag
        ## my @outadd=(); push @outadd, "cds" if(defined $OUTCDS); push @outadd, "cdna" if(defined $OUTCDNA);
        foreach my $tag (@outadd) {
          my $seq=  ($tag eq "cds") ? $bestorf->{sequence} : $cdnabest; 
          my $outh= ($tag eq "cds") ? $OUTCDS : $OUTCDNA;
          
          my $cutadd=($utrcut>0)? " cutlen=".length($cdnabest).",at=$utrcut;" : "";
          $seq  =~ s/(.{60})/$1\n/g;
          print $outh ">",$id, " type=$tag; aalen=$aalen,$pcds%,$compl;$cutadd clen=$cdnasize; strand=$strand; offs=$prostart5-$proend3; \n";
          print $outh $seq,"\n";
 
          if($utrorf) {
          my $uoseq=  ($tag eq "cds") ? $utrorf->{sequence} : $utrorfseq; 
          if($uoseq) {
            my( $uprostart5, $uproend3, $uorflen, $ustrand)= orfParts($utrorf, [qw( start stop length orient)]);
            my($aalen,$pcds,$compl) = proteindoc($utrorf,$cdnasize,$ustrand, 1);
            my $cutadd=($utrcut>0)? " cutlen=".length($uoseq).",at=$utrcut;" : "";
            $uoseq  =~ s/(.{60})/$1\n/g;
            unless($noutrorf) {
            print $outh  ">",$id,"utrorf", " type=$tag; aalen=$aalen,$pcds%,$compl;$cutadd clen=$cdnasize; strand=$ustrand; offs=$uprostart5-$uproend3; \n";
            print $outh $uoseq,"\n";
            }
            $nsplit += 2; #not now# $ngood--;
          }
          }
        }
        
#.... old ....      
#         if($action =~ /cds|cdna/i) {
#           my $tag= ($action=~/cds/i) ? "cds" : "cdna";
#           my $seq= ($action=~/cds/i) ? $bestorf->{sequence} : $cdnabest; 
#           
#           my $cutadd=($utrcut>0)? " cutlen=".length($cdnabest).",at=$utrcut;" : "";
#           $seq  =~ s/(.{60})/$1\n/g;
#              # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
#           print $OUTCDS ">",$id, " type=$tag; aalen=$aalen,$pcds%,$compl;$cutadd clen=$cdnasize; strand=$strand; offs=$prostart5-$proend3; \n";
#           print $OUTCDS $seq,"\n";
#           if($utrorfseq) {
#             my( $uprostart5, $uproend3, $uorflen, $ustrand)= orfParts($utrorf, [qw( start stop length orient)]);
#             my($aalen,$pcds,$compl) = proteindoc($utrorf,$cdnasize,$ustrand, 1);
#             # my $usize= length($utrorfseq);
#             my $cutadd=($utrcut>0)? " cutlen=".length($utrorfseq).",at=$utrcut;" : "";
#             $utrorfseq= $utrorf->{sequence} if($action=~/cds/i);
#             $utrorfseq  =~ s/(.{60})/$1\n/g;
#             unless($noutrorf) {
#             print $OUTCDS  ">",$id,"utrorf", " type=$tag; aalen=$aalen,$pcds%,$compl;$cutadd clen=$cdnasize; strand=$ustrand; offs=$uprostart5-$uproend3; \n";
#             print $OUTCDS $utrorfseq,"\n";
#             }
#             $nsplit += 2; #not now# $ngood--;
#             }          
#           }
          
      }
      
    } else { # move this outformat above with >fasta

			if($revaa_report and $revaalen >= $MINAA and $revinfo) {
    		$compl.=":$revinfo" ; # HACK, not a good place, add Notes column before orfprot ?
    		# revinfo == "aarev=99,90%,complete,-,9-309"
			}
			
      print $OUTAA join("\t", $id, $aalen, $pcds, $compl, $cdnasize, $strand, $prostart5, $proend3, $orfprot),"\n"; 
      if($utrorf and not $noutrorf) {
        my( $orfprot, $prostart5, $proend3, $orflen, $strand)= orfParts($utrorf, [qw(protein start stop length orient)]);
        my($aalen,$pcds,$compl) = proteindoc($utrorf,$cdnasize,$strand, 0);
        print $OUTAA join("\t", $id."utrorf", $aalen, $pcds, $compl, $cdnasize, $strand, $prostart5, $proend3, $orfprot),"\n"; 
      }
    }

  }
  # warn"#cdna_bestorf found $ngood good, $nsplit split, $nskip skip of $ncdna\n" if $DEBUG; 
  # fixme, write split count ngood
  return($ngood,$nsplit,$nskip);
}

__END__

### MOVED to cdna_proteins.pm ........................................................
# convert to array handling? _max(@list) 
# sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
# sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
#
# sub getBestProt {}
# 
# ## revise longest_orf_finder
# 
# sub getAllOrfs {}
# 
# sub get_orfs {}
# 
# sub identify_putative_starts {}
# 
# sub identify_putative_stops {}
# 
# sub revcomp {}
# 
# sub revcomp_coord {}
# 
# # parts from PASA/PasaLib/Nuc_translater.pm 
# use vars qw ($currentCode %codon_table);
# 
# # sub codon1 {}
# 
# sub translate_sequence {}
# 
# BEGIN {}
# 
# 1;
