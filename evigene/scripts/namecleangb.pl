#!/usr/bin/perl
# evigene/namecleangb.pl  from evigene2genbanktbl.pl

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/prot/"); # == evigene/scripts/

use strict;
use Getopt::Long;
use protein_names; # Evigene nameclean() etc

use constant VERSION  => '20130504'; # '20130326';  

## moved to protein_names; all evigene.conf options
# my $NAME_NONE = "Unknown|Uncharacterized|Hypothetical";  
# my $NAME_UNK  = "Uncharacterized protein"; # uniprot  
# my $NAME_UNKADDLOCUS= 1; # policy, add/not the gene id to end of UNK names
# ## my $MIN_ESTIDENT = 10; # min align% to keep  EST/Rna evidence ; not used for naming
# my $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
#   # -- note MIN_ID is used for both naming and keeping prot evidence, split this? use lower IDLIKE for evid?
#   # MIN_IDENTITY == MIN_NAMEIDENT now
# my $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
# my $MIN_CERTAIN  = 0; ## 60;  # putative naming, what?
# my $MIN_PROTIDENT= 0; # $MIN_IDENTITY; .. only for putative, no prot id info here...

our $NAME_KEEPID= 0; # change default
our $USE_TENAME; # in protein_names
my $FILECLEAN=0;
my $DEBUG=0;

my $optok= GetOptions(
#  "config=s", \$config,
#  "input=s", \$input, # for stdin vs files
#  "output=s", \$output,
  "KEEPID!", \$NAME_KEEPID, 
  "TENAMES!", \$USE_TENAME, 
  "FILECLEAN=s", \$FILECLEAN, 
  "debug!", \$DEBUG, 
  );

die "usage: app mygene.names > mygene.nameclean 
 opts: -fileclean=uniprot|cdd.names -KEEPID -TENAMES"
 unless($optok);
 
## sub cleannamefile() handles CDD also .. change to that?
## but specialized for uniprot/cdd naming table w/ id,size cols
if($FILECLEAN) {
  my $infile= $FILECLEAN; # or all of @ARGV;
  my $iscdd= ($infile =~ /cdd/i)?1:0;
  if($infile and -f $infile) {
    my($nrefname,$outnames)= cleannamefile($infile,"",$iscdd) ;
    warn "#cleannamefile(infile=$infile,,iscdd=$iscdd) => ($nrefname,$outnames)\n";
  } else {
    warn "#cleannamefile(infile=$infile,,iscdd=$iscdd) => missing infile\n";
  }
  
} else {  
while(<>) {
  unless(/^\w/) { print; } #? or bad names?
  chomp;
  #old# my($id,$name)= split"\t",$_,2; # maybe maybenot
  my($id,$name,@other)= split"\t",$_; # maybe maybenot
  unless($name and $id=~/^[\w.]+/) { $name=$id; $id=""; }  
  if(not $name or $name eq "na") { $name= $NAME_UNK; }
  
  ## pi is not present in group name sets: tair, uniref, omcl .. need to dig out of genes.gff
  # my $pi= ( $name =~ s/\s+\((\d+)%.*\)// ) ? $1 : $MIN_CERTAIN;  # trailing pctident
  ## pi: look in @other for '33%,nnn/nnn,nnn' ?
  my $pi= $MIN_CERTAIN;
  if( $name =~ s/\s+\((\d+)%.*\)// ) { $pi=$1; }
  elsif(@other) {
  	my($po)= grep /^\d+%/, @other; if($po) { ($pi=$po) =~ s/%.*//; }
  }
  
  my ($newna,$lowqualname,$diff)= nameclean( $name, $pi );
  print "$id\t" if($id);
  print $newna;
  ##print "\tlowqual=$lowqualname" if($lowqualname); ## FIXME: this adds column; need to preserve input @other
  print "\t",join("\t",@other) if(@other);
  print "\t", ($lowqualname ? "lowqual=$lowqualname" : ""); # add 1 col always at end
  print "\t", (($diff) ?  "update" : "same") if($DEBUG);  
  print "\told:$name" if($DEBUG and $diff > 1);
  print "\n";
}
}

#.........

=item   moved subs to evigene/scripts/prot/protein_names.pm
  see also nameclean.perl
  sub nameclean {}

=cut

=item add config

# evigene_config($config, \@configadd); # always even if $config null
# evigene_cacao3_gbsubmit.conf
#   nameless    = Unknown|Uncharacterized|Hypothetical  
#   nameunknown = Uncharacterized protein   # uniprot 2011
#   nameunkaddid = 0 # turn off Unc.. locus Thecc11111 addition for gbsubmit complaints
#   nameuncertain = putative # at end; 
#   pctuncertain  = 60  #  min ident% for nameuncertain
#   pctunknown    = 35  #  min similar% for protein evid align/equivalence, for naming; JCVI uses 35% MIN_PROIDENT
#   pctproevidence = 10 #  min similar% for protein evid ; not used for naming 
#   pctrnaevidence = 10 #  min align% to keep  EST/Rna evidence ; not used for naming MIN_ESTIDENT
#   nameidpatt  = (Os|At|AT)\\d{1,2}g\\d{3,} << fix, see below
#   namedrops   = Arabidopsis|thaliana|yeast|complete sequence|complete|genome|pseudogene
 
=cut

