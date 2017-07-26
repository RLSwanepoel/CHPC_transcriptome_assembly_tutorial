#!/usr/bin/perl
# cds2aadiff.pl : calc genes.cds2aa <> genes.aa 
# cat $genes.cds2aa $genes.aa  | cds2aadiff.pl  > $genes.diffa
# fixme: cds2aadiff.pl -cds $genes.cds2aa -aa $genes.aa > $genes.diffa

=item
  faTrans pub3i.cds stdout |\
  perl -ne 'if(/^>(\S+)/){$d=$1; puta() if($g); $g=$_; $aa=""; } \
  else { chomp; $aa.=$_; } END{puta();} sub puta{ $aa=~s/Z$/\*/; $aa=~s/(.{60})/$1\n/g; \
  $aa.="\n" unless($aa=~m/\n$/); print $g,$aa; }'\
   > pub3i.cds2aa     
=cut


$SHOWOK=$ENV{ok}||0;
$CUTP= $ENV{cut} || 'aalen|Name|oid';
@ERRS=qw(OK DIFF SMAL BIGG MISS MISC);

while(<>) {
  if(/^>(\S+)/){ $d=$1;  puta(); $g=$d; $ina=1 if(/$CUTP/); $aa=""; } 
  elsif(/^\S/) { chomp; $aa.=$_; } 
}
puta(); 
cmiss();
print "#diffaa: ", (join", ",map{ "n$_=".($nerr{$_}||0) } @ERRS),"\n"; 
 
sub cmiss {
  foreach my $g (sort keys %ca) { 
    next if($gotca{$g});
    $ca=$ca{$g};  $lca=length($ca);
    $et="MISC"; $er="$et"; $nerr{$et}++; 
    print "$er\t$g\t$lca\n" if($lca>9);   
  }
}

sub puta { 
  return unless($g);
  unless($ina){ $ca{$g}=$aa; } # cds2aa here
  else{ 
    $ca=$ca{$g}; # MISS handled
    if($ca ne $aa){ $ca=~s/Z/\*/g; 
    ($wca,$waa)= map{ length($_) } ($ca,$aa); $dap=""; $da= $wca - $waa; 
    ($xca,$xaa)=($ca,$aa); map{ s/\*$//;} ($xca,$xaa); 
    unless($ca){  $et="MISS"; $er="$et"; $nerr{$et}++; }
    elsif( $da>0 and ($i=index($xca,$xaa)) >= 0) {  $et="BIGG";  $er="$et.$i,$da"; $nerr{$et}++; $dap=difa(0,$aa,$ca); } 
    elsif(($i=index($xaa,$xca))>=0) { $et="SMAL";  $er="$et.$i,$da"; $nerr{$et}++; $dap=difa($i,$aa,$ca); } 
    else{ $et="DIFF"; $er="$et.$da"; $nerr{$et}++; $dap=difa(0,$aa,$ca); } 
    print "$er\t$g\n$dap";  $gotca{$g}= $et;
  } else { $et="OK"; print "$et\t$g\n" if($SHOWOK); $gotca{$g}= $et; $nerr{$et}++; } 
  } 
} 

sub difa { 
  my($io,$aa,$ca)=@_; my @ca=split"",$ca; my @aa=split"",$aa;
  $co=($io>0)?$io:0; unshift(@ca, ("-") x $co) if($co);  @ba=@ca; for my $i (0..$#aa) {
  my ($e,$c)=($aa[$i],$ca[$i]); $ba[$i]="." if($e eq $c); } my $ba=join"", @ba; 
  return "Da: $aa\nDc: $ba\n"; 
} 


