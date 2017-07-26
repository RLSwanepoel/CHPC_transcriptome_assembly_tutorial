#!/usr/bin/perl
# gmapgffix.pl

$notarget=1; #opt?
$src=$ENV{src};
while(<>) {
  s/Name=[^;]+;//; s/^###/#/; s/\.(path|mrna)(\d+)/p$2/g; 
  s/\t\S+/\t$src/ if($src);
  if(/\tgene/){s/.*//;} 
  elsif(/\t(CDS|exon)/) { 
    s/(ID)=[^;\n]+;?//g; 
    s/Target=/targ=/ if($notarget);
    } 
  elsif(/\tmRNA/){ 
    s/Parent=/gene=/; ($cv,$pi)=m/Coverage=([^;\s]+);Identity=([^;\s]+)/; ($d)=m/ID=([^;\s]+)/; 
    $cp=int($cv * $pi/100); s/\t\.\t/\t$cp\t/; $d=~s/p\d+$//; 
     ## $s=$sc{$d}; s/$/;vscore=$s/ if($s); 
  }
  print;
}


=item filter junk models

## use CDS:1 score if .gff has CDS; dont use score: unless that is meaningful
#...  -mrna mRNA -exon 'CDS,exon'

 gzcat $rnagenes.gff.gz | $evigene/scripts/overbestgene2.perl -in stdin \
 -alttr -typeover exon -noOVEREXON2 \
 -scoretype='many.score:3,CDS:2,UTR:1' \
 -genescore  -trivial 10  -pctover 10  -summarize -noskip \
  > $rnagenes.best1.gff

=cut

