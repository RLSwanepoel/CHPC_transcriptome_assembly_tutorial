#!/usr/bin/perl
# arp_condef.pl

=item notes

  make consensus def for gene groups
  cut from arp-blastpgp.pl
  
  inputs:
    aabugs3_omclgn.tab : table of ARPid, geneid
    *.aa.deflines.gz (just grep ^> from .aa )

  gunzip -c pro3/*.aa.deflines.gz | cat aabugs3_omclgn.tab  - \
   | perl arp_condef.pl > aabugs3_omclgn.consensus_def.txt
  
  
=cut

use strict;

my $poorannots = 'xxxx'; ##'amel|apis';
my $goodannots = 'culex|aedes|ixodes|pediculus';
my $debug=0;
my (%gdef, %ar);
my $GTAG = $ENV{gtag} || "ARP";
my $OGPRE= $ENV{idprefix}||"${GTAG}1_G";

while(<>) {
  if(/^>(\S+)/) {
    my $gn=$1; my $an=$_; 
    $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    $gdef{$gn}= $an;
    
  } elsif(/^$GTAG(\d+)\s+(\S+)/) { 
    my $ar=$1; my $gn=$2; 
    push( @{$ar{$ar}}, $gn); # $gg{$gn}=$ar; 
  }
}

foreach my $ad (sort{$a <=> $b} keys %ar) {
  my @gdef= map{ $gdef{$_} } @{$ar{$ad}};
  my $condesc= get_condesc( $OGPRE.$ad, \@gdef);
  print $condesc,"\n";
}


=item condesc / consensus description

  Still need to clean out species names here:
  ARP1_G15.consensus Nasonia vitripennis hypothetical protein  (LOC100113581); src=nasonia_NCBI_hmm8074
  ARP1_G18.consensus Acyrthosiphon pisum hypothetical protein  (LOC100168668); src=acyr1_ncbi_hmm118344

  And drop (LOCnnnn) and like IDs, but leave in some (xxxx) descriptions

=cut

sub get_condesc {
  my($cluid, $genedb)= @_;
  my $condesc="$cluid.consensus";
  
  my @df= @$genedb;
  my (%didspp,%de);
  foreach (@df) {
    next if(/by Gnomon/); #?? or not
    #  Gene predicted by Gnomon |  Partial gene predicted by Gnomon ..
    ##next if(/ixodes|pediculus/); # hack
    chomp;
    s/^>(\S+)\s+//; 
    my $id=$1;
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : $id;

    s/PREDICTED:\s*//i;
    s/\w+ \w+ hypothetical/hypothetical/i;
    s/\w+ \w+ similar to\s*/ /i;
    s/similar to\s*/ /i;
    s/\(LOC\d+\)//; 
    s/\,\s.*//;  
    s/\[[^\]]+\]//; # [Drosophila melanogaster] ...
    s/\(Fragment\)//i;
    
    #  $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    if(s/ID=FBpp\w+;?//) { s/,?FBtr\d+//; s/parent=FBgn/geneid=FBgn/; }
    s/gid=\w+;?// if($spp =~ /tribolium/);
    
    my $de= $_;
    my $dbx= ($de =~ s/dbxref=(\S+)//i) ? $1 : "";
    if($de =~ s/desc=//) { $de =~ s/(\w+)=(.+)$//; } # extra jack
    $de =~ s/^\s+//; $de =~ s/\s+$//;  $de =~ s/  +/ /g; $de =~ s/;\s*;/;/g;
    $de =~ s=\w+ \w+/== if($spp =~ /daphnia/);
    $de =~ s=\,.+$== if($spp =~ /daphnia/); #?

    $de{$id}= $de if($de =~ /\w+/); # and not $didspp{$spp} #?? change: keep all/spp but downweight 2+ below by 0.3?
    $didspp{$spp}++;
  }
  
  #now what? need one readable/valid line from these; not all agree
  my (%kw, %kwd, %wtd);
  foreach my $id (sort keys %de) {
    # warn "# $de{$id}\t[$id]\n" if $debug;
    my $w= $de{$id};  # $w =~ s/[,;=].*$//;
    my @w=  split /\W+/, $w;
    ## $kw{$w} += @w; $kwd{$id}{$w} += @w; # score whole phrase? 
    foreach (@w) { 
      next if(/^CG\d/ || /^AGAP\d/); 
      $_= lc $_; $kw{$_}++; $kwd{$id}{$_}++; 
      }
  }

  my $nspp= scalar(keys %didspp);
  %didspp=();
  my @wt= sort{ $kw{$b} <=> $kw{$a} } keys %kw;
  foreach my $id (sort keys %kwd) {
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : $id;
    my $wt=0; my $wi= @wt;
    my $nw= int( 0.7 * $wi);
    $wi *= 3;
    foreach my $w (@wt[0..$nw]) { 
      $wt += $wi if( $kwd{$id}{$w} ); 
      $wi -= 3;
      }
    $wt = $wt * 0.3 if($didspp{$spp}++);  
    $wt = $wt * 0.4 if($id =~ /$poorannots/); # down/up weight others?
    $wt = $wt * 1.2 if($id =~ /$goodannots/); # down/up weight others? : add more good annots
    $wtd{$id}= $wt;
  }
  
  # should require at least 2 cases of agreement, otherwise decline consensus desc
  my($bestid)= sort{ $wtd{$b} <=> $wtd{$a} } keys %wtd;
  my $de= $de{$bestid};
  $de =~ s/;\s*$//;
  if($de=~/\w/){ $condesc .= "\t$de; src=$bestid"; } else { $condesc .= "\tna"; }
  
  # warn "# $condesc\n" if $debug;
  return $condesc;
}


