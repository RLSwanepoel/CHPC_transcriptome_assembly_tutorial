# cdna_proteins.pm

# package cdna_proteins;
package main;

use strict;
use warnings;

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (translate_sequence get_protein reverse_complement);

use vars qw ($DEBUG $MINAA $MINEXON $MINGOOD $USEGOODLEN 
      $AA_cdna_GT_genome $ORF_FULLvPART $KEEPSAMECDS $NoStopCodon
      $pCDSbad $pCDSpoor $MINUTRORF  $USESelenocysteine
      );

our $DEBUG=1;
our $MINAA= 30;  # used 
our $MINEXON= 60; # for cut
our $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
# our $USE_CDSEXONS = 0;
our $USEGOODLEN=1;
our $AA_cdna_GT_genome= 1.50;  # for prot(cdna) > prot(genome) test; option?
# FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
our $ORF_FULLvPART = 0.85;
our $KEEPSAMECDS= 0; # prefer keep same CDS exons but can extend/shorten protein bounds
our $NoStopCodon=0;
    ## need global params for utrbad/poor
our $pCDSbad = 35;
our $pCDSpoor= 60;
our $MINUTRORF= 300; #?

# our @stop_codons = qw(TAA TAG TGA); # allow changes, esp TGA => SelC
# parts from PASA/PasaLib/Nuc_translater.pm 
use vars qw ($currentCode %codon_table @stop_codons %amino2codon_table);


# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


sub bestorf_test
{
  my($bestorf,$nextorf) = @_;
  my ($bestprot)= orfParts($bestorf);
  my ($nextprot)= orfParts($nextorf);
  my $MinCDS= 3*$MINAA;

       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
  
  my $oki= ( $nextprot =~ /\w/ 
    && ($nextorf->{goodlen} >= $MinCDS)  # dang, goodlen,leng are cdslen not aalen
    && ($nextorf->{goodlen}/$nextorf->{length} >= $MINGOOD)) ? 1 : 0;
    
  if($oki) {
    if( $bestprot ) {
      ## problem here for rev=0,1 missing huge partial rev for short fwd .. but get huge if just rev
      my $arat = $nextorf->{goodlen} /  $bestorf->{goodlen}; # ok,> MINAA
      my $adiff= $nextorf->{goodlen} -  $bestorf->{goodlen};  

      if( $arat > 1.0 and  $arat < $AA_cdna_GT_genome and $bestprot =~ m/XX/) { ##  $bestorf->{goodlen}/$bestorf->{length} < 0.99
        (my $bp= $bestprot) =~ s/X/./g; 
        # if(length($bp) > length($nextprot) { } # chomp some?
        return(1, $nextprot, $nextorf) if($nextprot =~ m/$bp/);
      }
      
      if( $arat > $AA_cdna_GT_genome
        or ( $adiff >  0 and ($nextorf->{complete}==3 or $bestorf->{complete} <3) ) 
        or ( $adiff >= 0 and $nextorf->{complete}==3 and $bestorf->{complete} <3 ) )
       { 
        return(1, $nextprot, $nextorf);
       }
    } else {
      return(1, $nextprot, $nextorf);
    }
  }
  return( 0, $bestprot, $bestorf);
}

sub utrorf_test
{
  my($bestorf, $allorfs, $cdnasize) = @_;

  my ($utrorf,$utrosize)= getUtrOrf($bestorf, $allorfs, $cdnasize);  
  if ($utrorf) {
    my ($orfok) = bestorf_test(undef,$utrorf);
    return($utrorf,$utrosize) if($orfok);  
  }
  return(0,0);  
}

sub revorf_report
{
  my($bestorf, $allorfs, $cdnasize, $issorted) = @_;
  my ($revorf,$revosize)= getRevOrf($bestorf, $allorfs, $cdnasize, $issorted);  
  if ($revorf) {
    my($aalen,$pcds,$compl,$orflen,$fahead)  
     	= proteindoc($revorf,$cdnasize); # ,$strand,0
		 ## note: fahead= "aalen=99,90%,complete; clen=350; strand=-; offs=9-309;";
		 ## aarev=99,90%,complete,strand:-,offs:9-309;
		 ## aarev=99,90%,complete,-,9-309;  << use this?  or aarev=99,...,9-309:-; 
		 ## add pct-of-longest to report?  maybe drop strand,offset to simplify
		 ## aarev=55%,aa99,90%,complete,s-,o9-309;  << use this?  or aarev=99,...,9-309:-; 
		my $longsize= $bestorf->{goodlen} || 1; # revosize == goodlen; orflen == length
		my $plong= int(0.5 + 100*$revosize / $longsize);
		$fahead =~ s/aalen=/aa/; $fahead =~ s/clen=d+;//; 
		$fahead =~ s/strand=/s/; $fahead =~ s/offs=/o/; 
		$fahead =~ s/=/:/g; $fahead =~ s/ //g; $fahead =~ s/;\s*$//;  $fahead =~ s/;/,/g; 
		return ($aalen,"aarev=$plong%,$fahead"); #?
	}
  return(0,0);  
}

sub proteindoc
{
  my($orf, $cdnalen, $cdnarev, $forFASTA) = @_;
  my($aalen, $compl, $pcds) = (0) x 10;
  my( $orfprot, $prostart, $proend, $orflen, $orient) = orfParts($orf);
  $cdnarev ||= $orient; # fill in blank; use always? shouldnt need to pass cdnarev as param.
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  if($orfprot) {
    $pcds  = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
    my $urev= ($prostart>$proend)?1:0;
    my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
    my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;
    
    if( $NoStopCodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    $aalen-- if(substr($orfprot,-1) eq '*'); # annoyance
    $compl= $orf->{complete};
    $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
    ##? not bad if partial? if u1len or u2len == 0
    ## need global params for utrbad/poor
    if($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
    elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
    
    if($forFASTA) { $orfprot =~ s/(.{60})/$1\n/g; } #? only for fasta output
  } else {
    $orfprot="X"; # ? make dummy orfprot?
  }

  ## fixme: add goodlen or gaps count to aadoc    
  # my $aagap= $orf->{length} - $orf->{goodlen};
  # my $aagap= $orfprot =~ tr/X/X/;
  
  my $fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;";
  return($aalen,$pcds,$compl,$orflen,$fahead,$orfprot);
}



sub getBestProt2
{
  my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my($orfprot,$prostart5,$proend3)=("",0,0,"");
  
  # fix this to return longest full prot, and longest partial (if longer)
  # .. test which is best.
  # FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
  #   my $longorf= $longest_orf_finder->get_longest_orf($cdna); # == hash
  # FIXME3: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
  
  # my ($longorf,$longfull,$orfs) = getAllOrfs($cdna); 
  
  my $strands= ($ptype =~ /(both|fwd|for|rev)/i)? $1 : "fwd";
  my ($longorf,$longfull,$orfs) = getAllOrfs($cdna, $strands);  
  #  where strands=null/fwd/rev/both, and strands=rev means do revcomp on cdna
  # Ooops, orf->{start,stop} are reversed for strand=rev; start>stop
  
  if(ref($longorf)) {
    # $orfprot= $longest_orf_finder->get_peptide_sequence();
    # ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();

    my $lookmore= ($ptype =~ /long/)?0:1;
    if($KEEPSAMECDS and $oldStart_b > 0) { # may not be right yet.
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
      my $lsize= $longorf->{$keylen}; 
      my $fsize= $longfull->{$keylen};
      
      # FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
      if($fsize >= $ORF_FULLvPART * $lsize) {
        if(ref($exongff)) {
        my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff2( $exongff, $longorf);
        my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff2( $exongff, $longfull);
        $longorf= $longfull if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
        } else {
	      $longorf= $longfull;
	      }
      }
    }
    
    ($orfprot,$prostart5,$proend3)= orfParts($longorf);
  }

  return($longorf, $orfs); # return allorfs now
  # old# return($orfprot, $prostart5, $proend3, $longorf, $orfs); ##, $utrorf 
}

sub getBestProt # old 
{
  # my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my ($longorf, $orfs)= getBestProt2( @_ );
  my ($orfprot,$prostart5,$proend3)= orfParts($longorf);
    # my $getutrorf=1; # always test for utr protein? moved out of here
    # my ($utrorf,$utrosize)= getUtrOrf($longorf, $orfs, length($cdna)); # ($getutrorf) ? xxx() : (0,0);
    # ^^ remove from here, add to utrorf_test
  return($orfprot, $prostart5, $proend3, $longorf); ##, $utrorf 
}

sub getUtrOrf
{
  my ($longorf, $orfs, $cdnalen)= @_;
  my ($utrorf,$utrosize)=(0,0);
  my $lsize= $longorf->{length};
  # my $lgood= $longorf->{goodlen};
  return($utrorf,$utrosize) unless(($cdnalen - $lsize >= $MINUTRORF));  # test even if lsize/cdna > 60% ?  
  # my $dotest= (($cdnalen - $lsize > $MINUTRORF) || ((100*$lsize/$cdnalen) < $pCDSpoor));
  #  ^^ this is bad test; some large tr with pCDS > pCDSpoor are hiding other orf, usually in 3+utr exons on 1 end
  if(1) { 
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
    my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
    foreach my $orf (@$orfs) {  # orfs can be rev of longorf, start/stop are fwd tho
      my($ob,$oe,$ogood,$osize)= ($orf->{start},$orf->{stop},$orf->{goodlen},$orf->{length},);
      ($ob,$oe)=($oe,$ob) if($ob>$oe);
      if(($ob > $le or $oe < $lb) and ($ogood>=$MINUTRORF) and $ogood>$utrosize) { # size or $osize > 0.5*$lsize ??
        $utrorf= $orf; $utrosize= $ogood;        
        }
      }
    }
  return($utrorf,$utrosize);    
}

sub getRevOrf # or revorf_test  # find/report on longest/best orf in reverse of longorf
{
  my ($longorf, $orfs, $cdnalen, $issorted)= @_;
  my $MinCDS= 3*$MINAA; $issorted||=0;
  my ($revorf,$revosize)=(0,0);
  my $lorient= $longorf->{orient};
  # my $lsize= $longorf->{length}; # goodlen ?
	##  orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
	## my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
	foreach my $orf (@$orfs) {  # expect long sorted, or check all? 
		my($oorient,$ogood,$osize)= ($orf->{orient},$orf->{goodlen},$orf->{length});
		# $ob,$oe, = $orf->{start},$orf->{stop}, # ($ob,$oe,$oorient)=($oe,$ob,'-') if($ob>$oe);
		if($oorient ne $lorient and $ogood>=$MinCDS and $ogood > $revosize) {  
			$revorf= $orf; $revosize= $ogood; 
			last if $issorted;
			}
   }  
  return($revorf,$revosize);    
}

# replace old sub getCDSgff($exons,$orfprot,$prostart5,$proend3) w/ getCDSgff2($exons,$orf)
sub getCDSgff2
{
  my($exons,$orf)= @_;
  
  # my($exons,$orfprot,$prostart5,$proend3,$addutr) = @_;
  my ($orfprot,$prostart5,$proend3)= orfParts($orf);
  # Ooops, orf->{start,stop} are reversed for strand=rev; start>stop
    
  my ($cds5,$cds3)=(0,0);
  my @cds= ();
  my @utr= ();
  ## for phase; need reverse @exons
  my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3, $strand)= (0) x 10;
  $cdna1= 1;
  $nt_length= 0; # $prostart5 % 3; #??
  
  # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
  # .. must make next exon(if exists) part of CDS stop
  
  foreach my $exon (@$exons) {
    my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
    $strand= $xo;
    ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
    
    my $xd= abs($xend3 - $xend5); # ?? +1 for width
    
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
      
    } elsif(1) { 
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
    
    $cdna1= $cdna2+1;
  }       
  
  
  my $trlen= $cdna1; # = $nt_length
  use constant forFASTA => 0;
  my($aalen,$pcds,$compl,$orflen) ## ,$fahead,$orfprot2
     = proteindoc($orf,$trlen,$strand,forFASTA);
  
#  #....
#   my $aalen=length($orfprot); 
#   $aalen-- if(substr($orfprot,-1) eq '*');
#   my $clen= $aalen * 3; # can be off by -1,-2 here. 
#   my $ap=int(100 * $clen/$trlen);
#  #....

    
  my $mattr="cxlen=$orflen/$trlen;aalen=$aalen,$pcds%,$compl";
  $mattr.= ";protein=$orfprot" if($orfprot);
#   if($orfprot) {
#     my $p5= (substr($orfprot,0,1) eq 'M')?0:1;
#     my $p3= (substr($orfprot,-1,1) eq '*')?0:1;
#     my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
#     $mattr.= ",$prostat;protein=$orfprot";
#   }
  
  $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
  ## mattr keys: cxlen,aalen,protein,utrx
  
  # FIXME: resort @cds by loc, not reversed. : let caller do
  # @cds = sort _sortgene @cds;

  ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
  ## ($cdslong, $attrL, $pcodeL, $maxutrL)
  return (\@cds, $mattr, $pcds, _max($nu5,$nu3), \@utr); 
}

  
sub getAllOrfs { # added strands to return both, or rev
  my ($input_sequence, $strands,$flags) = @_;
  $strands ||= "fwd";
  $flags ||="";
  # add fwd,rev orfs here ? gff-caller wants only 1 strand at a time, need option 
  
  return undef unless ($input_sequence or length ($input_sequence) >= 3) ;
  $input_sequence = uc ($input_sequence); # was lc() change all to uc()
  
#  ## fixme: this screws seq indices: start,stop; instead of chomp, restrict @starts,@stops
#   if($flags =~ /dropnnn|chomp/i) {
#     $input_sequence =~ s/^N+//;
#     $input_sequence =~ s/N+$//; 
#   }

  my (@starts, @stops, @orfs);
  unless($strands =~ /^(r|\-)/) {  # forward_strand_only();
    @stops  = identify_putative_stops($input_sequence);
    @starts = identify_putative_starts($input_sequence,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $input_sequence, '+');
    }
  if($strands =~ /^(b|r|\-)/) { # rev|r|both
    my $revseq= revcomp($input_sequence); 
    @stops  = identify_putative_stops($revseq);
    @starts = identify_putative_starts($revseq,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $revseq, '-');
  }

  if (@orfs) {  # get both longest complete, longest partial
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

sub orfParts
{
  my($orf,$parts) = @_;
  if(ref $orf) {
    # $parts= [ qw(protein start stop length) ] unless(ref $parts); 
    if(ref $parts) { return @$orf{ @$parts }; }
    return($orf->{protein}, $orf->{start}, $orf->{stop}, $orf->{length}, $orf->{orient}); # added, $orf->{orient}
  } else {
    return("",0,0);
  }
}


sub get_orfs {
  my( $starts_ref, $stops_ref, $seq, $direction) = @_;
  
  # $direction here only affects start,stop relative to forward cdna seq .. use it.
  # .. input seq is already rev(seq) for -dir
  # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  # .. fix some places, leave that way others?
    
  # unless ($starts_ref && $stops_ref && $seq && $direction) { warn "get_orfs: params not appropriate"; }  
  # want only max orf, complete + partial
  # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
   
	my %last_stop_part = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
	my %last_stop_full = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
  my @orfs;
  my $seq_length = length ($seq);
  my $norf=0;
  
  ##  $seq .= "###" if($DEBUG); # why're we getting bad prot at partial end? isnt prot but cds range needs adjust for -2,-1
  
  foreach my $start_pos (@{$starts_ref}) {
		my $start_pos_frame = $start_pos % 3;
		# includes partial5: starts at 0,1,2 unless real start
		# ** PROBLEM? -- for internal mis-assemblies, should
		#   test partial-starts also following each stop? dont require ATG start inside?
		
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
				  print STDERR "#ERR: innerstop\n"; #? last unless($DEBUG);
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
        my $nxxx= $protein =~ tr/X/X/;
        ##my $isbad= ($nxxx > 0.5*length($protein)) ? 1: 0;
        my $goodlen= $orflen - 3*$nxxx;
        $norf++;
        # opt to drop *stopcodon, here or where, 
        
        if($DEBUG>1){ ##my $aalen=length($protein);
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
  #global# my @stop_codons = qw(TAA TAG TGA);
## add frameshift detect option here?  only if remaining seq >> nnn, only if %cds/utr falls below pCDSpoor ?
## need indel max to test: offby -2,-1,1,2 only to shift $i
  
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


sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}

sub backtranslate_protein {
  my ($sequence) = @_;
  $sequence = uc ($sequence); # redundant ??
  my $seq_length = length ($sequence);
  my $cds_sequence="";
  for (my $i = 0; $i < $seq_length; $i++) {
      my $codon;
      my $amino = substr($sequence, $i, 1);  # deal with non-alpha
      if (exists($amino2codon_table{$amino})) {
        $codon = $amino2codon_table{$amino};
      } elsif($amino =~ /[A-Z]/) {
        $codon = 'NNN'; # fixme
      } else {
        $codon= ''; # eat it?
      }
      $cds_sequence .= $codon;
  }
  return($cds_sequence);
}
      

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

## add frameshift detect option here?  if codon = stop_codons ; need it above other places
## need indel max to test: offby -2,-1,1,2 only to shift $i

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

sub useSelenocysteine
{
  @stop_codons = grep !/TGA/, @stop_codons; #  qw(TAA TAG); # not TGA
  $codon_table{'TGA'} = 'u';
  $currentCode = "universalSelC";
  #from# $USESelenocysteine=1;
}

BEGIN {
  $currentCode = "universal";
  @stop_codons = qw(TAA TAG TGA);
  # stops: TAG  TGA  TAA; Note: TGA also == Selenocysteine valid non-stop
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
  TGT => 'C', TGC => 'C', 
    TGA => '*',  # alternate UGA => 'U' Sec / SeC Se-Cys
  TGG => 'W',
  CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',  CGN => 'R',
  AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R',
  GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',  GGN => 'G', 
  );
  
  # $USESelenocysteine=0;
  useSelenocysteine() if($USESelenocysteine);

  # should use freq/codon or other choice of best codon per aa.
  foreach my $codon (sort keys %codon_table) {
    my $aa= $codon_table{$codon};
    $amino2codon_table{$aa}= $codon unless($amino2codon_table{$aa}); # multiple codons/aa, what to do?
  }
  
}

1;
