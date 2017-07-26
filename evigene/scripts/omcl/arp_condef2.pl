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
  
 Fish cleans:
  -- see also evigene2genbank for more name cleans per NCBI Genbank requirements
perl -pe 's/; src=.*//; 
s/\s+homolog//; 
s/H.sapiens //i; s/human and mouse\s*//i; s/human\s*//i; 
s/Novel\s*//; s/Putative\s+//i; 
s/\s*\([^\)\n]+\) *$//;  # many (Species) trailers in these names
s/\s*[\dA-Z]+$//;  # leave this trailing alphanum, cut for viewx
  
=cut

use strict;

my $debug=0;
my $poorannots = $ENV{poor} || 'xxxx'; ##'amel|apis';
my $goodannots = $ENV{good} || 'yyyy'; # 'culex|aedes|ixodes|pediculus';
my (%gdef, %ar);
my $GTAG = $ENV{gtag} || "ARP";
my $OGPRE= $ENV{idprefix}||"${GTAG}1_G";
my $NAME_NONE = "Unknown|Uncharacterized|Hypothetical";  
my $NAME_UNK  = "Uncharacterized protein"; # uniprot  

while(<>) {
  if(/^>(\S+)/) { #?? require >
    my $gn=$1; my $an=$_; 
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    if($an=~/\t/) { my @an=split"\t",$an; $an=join("\t",@an[0,1]); }
    $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    $gdef{$gn}= $an;
    
  } elsif(/^$GTAG(\d+)\s+(\S+)/) { 
    my $ar=$1; my $gn=$2; # dang: species_gene here, species:gene names
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    push( @{$ar{$ar}}, $gn); # $gg{$gn}=$ar; 

  } elsif(/^(\w\S+)\s+/) { # more names; may be \t table, split off \textra ?
    my $gn=$1; my $an=$_;
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    if($an=~/\t/) { my @an=split"\t",$an; $an=join("\t",@an[0,1]); }
    $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    $gdef{$gn}= $an;
  }
}

foreach my $ad (sort{$a <=> $b} keys %ar) {
  my @gdef= map{ $gdef{$_} } @{$ar{$ad}};
  my $condesc= get_condesc( $OGPRE.$ad, \@gdef);
  print $condesc,"\n";
}


=item condesc / consensus description

  FIXME: use  evigene/scripts/prot/protein_names.pm

  Still need to clean out species names here:
  >> use evigene/scripts/namecleangb.pl
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
    my ($isunk,$isput,$islike)= (0,0,0);

    next if(/by Gnomon/); #?? or not
    #  Gene predicted by Gnomon |  Partial gene predicted by Gnomon ..
    chomp;
    s/^>//; s/^(\S+)\s+//; 
    my $id=$1; $id=~s/_/:/ unless($id=~/:/); #DAMMIT.need.for.all.now# 
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : $id;

    s/PREDICTED:\s*//i;
    s/\w+ \w+ hypothetical/hypothetical/i;
    s/\w+ \w+ similar to\s*/ /i;
    s/similar to\s*/ /i;
    s/\(LOC\d+\)//; 
    s/\,\s.*//;  
    s/\[[^\]]+\]//; # [Drosophila melanogaster] ...
    s/\(Fragment\)//i;
    
	# Fish cleans; #
    #below# s/\s+homolog//;
    s/H.sapiens //i; s/human and mouse\s*//i; 
    #below# s/Novel\s*//; s/Putative\s+//i;
    s/\s*\([^\)]+\) *$//;  # OPTION: many (Species) trailers in these names

	# evigene2genbank nameclean()
  # typos and britglish > amglish; need config list
  s/Uncharacterised/Uncharacterized/; 
  s/dimerisation/dimerization/;  s/luminium/luminum/g;  # Aluminium
  s/signalling/signaling/; # Two Ls in British English, one in American English. 
  s/onoxygenase/onooxygenase/; # [Mm]onox..
  s/sulphide/sulfide/ig; s/sulphur/sulfur/ig;
  s/Tumour/tumor/ig;  
  s/haemoprotein/hemoprotein/i;  s/\bhaem/heme/ig; # Quinohaemoprotein>Quinohemoprotein
  s/characteris/characters/;  
  s/proteine/protein/;  s/\bcomponenet/\bcomponent/;

  if(s/^(Predicted|Conserved|Expressed|Novel)\s+//i) { }  # maybe set $isunk?
  if(s/^(Putative|Possible|potential|probable)\s+//i) { $isput=1; }   
  s/(homology|similar|Similarity|related) to\s*//i;
  s/\s\(Fragment\)//; s/\s[Ii]soform\s*.*$//; #? leave in isoform xx ?
  unless(m/protein\-protein/) { s/^(ORF|protein)\s*//i; }

  s/\b(Arabidopsis|thaliana|yeast|human|Staphylococcal|complete sequence|complete|genome|pseudogene)[,\.\s]*//ig;
  s/\s*\((?:InterPro|TAIR):[\w\.-]+\)//ig; #  (TAIR:AT1G22000.1); (InterPro:IPR010678)
  if(s/paralog of //i) { $islike=1; }

  s/\s*[Hh]omolo(gy|gue|g)\s+\d+//g; s/\b[Hh]omolo(gy|gue|g)\s*//g; # ? set $isput ? set -like? # add 'ortholog'
  s/\s*[Oo]rtholo(gy|gue|g)\s+\d+//g; s/\b[Oo]rtholo(gy|gue|g)\s*//g; 
  s/^(of|with)\s+//ig; # bad leading words
  # Add protein to the end when ending in:  'binding|domain|like|related'
  s/\b(binding|domain|like|related)\s(\W*)$/$1 protein $2/;  
  if( s/[,\s]*putative//g ) { $isput=1; } #s/putative, putative/putative/;

  # punctuation
  s/[\|]/:/g; s/#/n/g; 
  # s/_/ /g; # or leave uscores ??
  s/\s*$//; s/^\s$//; # lead/end spaces
  s/\s*[\/\.,;_-]*$//; # trailing punc
  s/[.] /, /g; # no sentences?
  s/^\W+//; # no leading crap

  # unbalanced brackets; # add {} ? not used; <> ? not brackets
  if(/[\(\)\[\]]/) {
    my($nb,$ne,$d);
    $nb= tr/\[/\[/; $ne= tr/\]/\]/; $d=$nb - $ne;
    while($d>0) { s/\[//; $d--; } while($d<0) { s/\]//; $d++; }
    $nb= tr/\(/\(/; $ne= tr/\)/\)/; $d=$nb - $ne;
    while($d>0) { s/\(//; $d--; } while($d<0) { s/\)//; $d++; }
  }

  # SEQ_FEAT.ProteinNameEndsInBracket:  Phosphoenolpyruvate carboxykinase [ATP] << change [] for ncbi
  if(/\]$/) {  s/\[/\(/g; s/\]/\)/g;}
  $_=$NAME_UNK unless(/\w\w/); # 'Conserved protein' becomes blank

  $isunk= (m/^($NAME_NONE)/i)?1:$isunk; # ($pi < $MIN_IDENTITY or  m/^($NAME_NONE)/i)?1:$isunk; 
  $isput= (!$isunk and $isput)?1:0; ## ($pi < $MIN_CERTAIN or $isput))?1:0;

  if($islike) { unless(/family/ or /\blike/ or /protein kinase/) { s/\s+protein//; $_ .= "-like protein"; }  }
  s/protein protein/protein/ig; # other stutters?
  $_ .= ", putative"  if($isput and not $isunk);

  s/^([a-z])([a-z][a-z])/\u$1$2/; # upcase 1st let
  #......

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
  # bug: missing bestid..
  my($bestid)= sort{ $wtd{$b} <=> $wtd{$a} or $a cmp $b } keys %wtd;
  my $de= $de{$bestid};
  $de =~ s/;\s*$//;
  if($de=~/\w/){ $condesc .= "\t$de"; $condesc.="; src=$bestid" if($bestid); } else { $condesc .= "\tna"; }
  
  # warn "# $condesc\n" if $debug;
  return $condesc;
}


