#!/usr/bin/perl
# bam2bw.pl

use strict;

my $dolog= (defined $ENV{log}) ? $ENV{log} : 1;
my $doave= (defined $ENV{ave}) ? $ENV{ave} : 1;  # of n-bams
my $samtools=$ENV{samtools} || "samtools";
my $sopt= $ENV{opt} || "-d 49999";

my $dgenomesize= shift @ARGV;
#? add genome.faidx to get genome base vs read bases?
my @inbam= grep /\.bam/,  @ARGV;
my $inbam= join " ", @inbam;
my $nbam=@inbam;

### this is true, but the seq col has <<>> for introns (N) and other info to count..
# warn "FIXME: new samtools mpileup DOES NOT count reads/base; it is counting read span? i.e. NO INTRONS\n
# need older samtools pileup one.bam, or other tool to turn reads.bam to readcover.bed counts";

## mpileup opt -d maxreads-per-base-bam may be needed; default autosets (depend on bam count?) 
##    [mpileup] 7 samples in 7 input files : Set max per-file depth to 1142 
if($dgenomesize =~ m/\.fai$/ and $sopt !~ m/\-f/) { (my $dg=$dgenomesize) =~ s/.fai$//;  $sopt .= " -f $dg"; }

die "usage: bam2bw.pl /path/to/genome.chr_size.tab in1.bam [in2.bam ..]"
  unless(-f $dgenomesize and -f $ARGV[0]);

open(S,$dgenomesize) or die; 
my %rsize=(); while(<S>){ my($r,$n)=split; $rsize{$r}=$n; } close(S);

my $hasref= ($sopt =~ m/-f \S/) ? 1 : 0;
my($lr,$lb,$le,$lc,$li,$lm);
open(IN, "$samtools mpileup $sopt $inbam |") or die "samtools mpileup $sopt $inbam";
while(<IN>) {
  my($r,$b,$xn,@v)=split"\t"; 
  my ($c,$in,$mc)=(0,0,0); 

  for(my $i=1; $i<@v; $i+=3) { # i=1 offset to seq
    $mc += $v[$i-1];  # mpileup: triplets per in.bam file: count,seq,qual
    my $s= $v[$i];
    my($nc,$ni);
    if($hasref) { $nc = $s=~tr/.,/.,/; }
    else { $nc= $s=~tr/ACGTacgt/ACGTacgt/; }  # OR count "." if using genome.faidx
    $ni= $s=~tr/<>/<>/;
    $c += $nc;
    $in += $ni;
    }

  if($doave and $nbam>1) { map{ $_ = ($dolog) ? $_/$nbam : int($_/$nbam); } ($c,$in,$mc);  }
  if($dolog) {
    map{ $_= int(10 * log(1+$_))/10; } ($c,$in,$mc);
    # $c= int(10*log(1+$c))/10; $in= int(10*log(1+$in))/10; $mc= int(10*log(1+$mc))/10;
  }

  if($r eq $lr and $c == $lc and $in == $li and $le == $b-1 ) { $le=$b; } 
  else { putb($lr,$lb-1,$le,$lc,$li,$lm);  $le=$lb=$b; }  
  ($lr,$lc,$li,$lm)= ($r,$c,$in,$mc); 
} close(IN);

putb($lr,$lb-1,$le,$lc,$li,$lm); # last

sub putb { 
  if($_[0] && $_[2] > $_[1] && ($_[3] >= 0.5 || $_[4] >= 0.5)) { 
  my $s= $rsize{$_[0]}||0; return if($_[2] > $s or $_[1] >= $s); 
  print join("\t",@_),"\n"; }  
}

__END__

# merge subset.bam to one bed/bw set
for i in 3 4 16; do   { 
 for j in 1 8 C; do {
  brna="merge-$i-$j" 
  inbam=`ls $i-${j}R*-dmag2.bam`
  echo bam2bw $brna FROM $inbam
  
  bam2bw.pl $dgenomesize $inbam > $brna.bed
  ## nasty fails when read beyond end if dgenomesize.. filter
  bedGraphToBigWig $brna.bed $dgenomesize  $brna.bw
  
 } done
} done  
