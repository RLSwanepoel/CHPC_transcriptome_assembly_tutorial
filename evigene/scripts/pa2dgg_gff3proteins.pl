#!/usr/bin/env perl

## 2011.aug: mostly replaced by evigene/scripts/genefindcds.pl
## ...............................

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Longest_orf; # dgg add

my $usage = "\n\nusage: $0 gff3_file genome_db [prot|CDS|cDNA|gene|checkprot|addprot|addprotpart,default=prot] [flank=0]\n\n";
$usage .="from PASA misc_utilities/gff3_file_to_proteins.pl by Brian Haas with updates by D. Gilbert\n";
$usage .="# 2011.aug: See improvement at evigene/scripts/genefindcds.pl\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $seq_type = $ARGV[2] || "prot";
my $flank = $ARGV[3] || 0;

my ($upstream_flank, $downstream_flank) = (0,0);

if ($flank) {
	if ($flank =~ /:/) {
		($upstream_flank, $downstream_flank) = split (/:/, $flank);
	}
	else {
		($upstream_flank, $downstream_flank) = ($flank, $flank);
	}
}

if ($upstream_flank < 0 || $downstream_flank < 0) {
	die $usage;
}



unless ($seq_type =~ /^(prot|CDS|cDNA|gene|checkprot|addprot|addprotpart)$/) {
    die "Error, don't understand sequence type [$seq_type]\n\n$usage";
}

my $longest_orf_finder = undef;  # dgg add
if($seq_type =~ /^(checkprot|addprot|addprotpart)/) {
  $longest_orf_finder = new Longest_orf();
  $longest_orf_finder->forward_strand_only();
  # $longest_orf_finder->allow_partials(); # need this NO, bad; but want in few cases of truncated genes
  # $longest_orf_finder->allow_3prime_partials();
}

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
## dgg: this can chew up memory, don't need all at once; 
#  - read contig/scaf/chr ids from fasta_db: >(\S+)
#  - iterate those, w/ index_GFF3_gene_objs(gff, objhref, contigid)

my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my %params;
        if ($seq_type eq "gene") {
            $params{unspliced_transcript} = 1;
        }
        
        $gene_obj_ref->create_all_sequence_types(\$genome_seq, %params);
        
		
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
			my $orientation = $isoform->get_orientation();
			my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $isoform->get_gene_span();
			
            my $isoform_id = $isoform->{Model_feat_name};
            my $com_name = $isoform->{com_name} || "";
	    if ($com_name eq $isoform_id) { $com_name = ""; } # no sense in repeating it
            
            my $seq = "";

            if ($seq_type eq "prot") {
                $seq = $isoform->get_protein_sequence();
            }
            elsif ($seq_type eq "CDS") {
                $seq = $isoform->get_CDS_sequence();
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $model_lend, $model_rend, $orientation, \$genome_seq);
				}
			}
            elsif ($seq_type eq "cDNA") {
                $seq = $isoform->get_cDNA_sequence();
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $gene_lend, $gene_rend, $orientation, \$genome_seq);
				}
			}
            elsif ($seq_type eq "gene") {
                $seq = $isoform->get_gene_sequence();
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $gene_lend, $gene_rend, $orientation, \$genome_seq);
				}
			}

      # dgg: add "addprot" : make prot from cDNA, print gene model w/ CDS added, mRNA prot=AAA attribute
      
      elsif($seq_type eq "addprot" || $seq_type eq "addprotpart") {
          my $gffcdna=  $isoform->get_cDNA_sequence();
          my $trlen= length($gffcdna) || 1;
          # dont bother w/ shorties? trlen < 30?  <60? <120?
          my $gattr = $isoform->{pub_comment} || ""; # dgg hack put mrna attributes here
          my($part5,$part3)=(0,0);
          if($seq_type eq "addprotpart") { ($part5,$part3)=(1,1); }
          elsif($gattr =~ m/chimera/i) {
	          my($cn)= $gattr =~ m/chim(\d)=/; 
            my($pn)= $gattr =~ m/path=(\d+)/;
	          if($pn == 1) { $part3=1; } elsif($pn == 2) { $part5=1; }
            elsif($cn == 2) { $part3=1; } elsif($cn == 1) { $part5=1; }
	        }

=item parse this chimera mess for 5/3 partial

Scaffold105     gmapx   mRNA    320388  320621  17      +       .       ID=ACYPI000251-RA_C1;Target=ACYPI000251-RA 1 23
4;cov=15.8;match=259;nexon=1;pid=100.0;qlen=1483;path=1/3;chim2=Scaffold190:491825-499988:-;chimera=exon-exon boundary 
(sense) at 234;desc=RefSeq XM_001947776 PREDICTED: similar to prolyl 4-hydroxylase alpha subunit 1, putative (LOC100158
825), partial mRNA;cxlen=24/234;aalen=8,10%;protein=MLLNILMH*
1 234;ix=1
Scaffold190     gmapx   mRNA    491825  499988  84      -       .       ID=ACYPI000251-RA_C2;Target=ACYPI000251-RA 235 
1483;cov=84.2;match=1249;nexon=8;pid=100.0;qlen=1483;path=2/3;chim1=Scaffold105:320388-320621:+;chimera=exon-exon bound
ary (sense) at 234;desc=RefSeq XM_001947776 PREDICTED: similar to prolyl 4-hydroxylase alpha subunit 1, putative (LOC10
0158825), partial mRNA;cxlen=123/1249;aalen=41,9%;protein=MTLFMTKRLKLLLIWLRRTCPTQHITLTVRLLCWTTNVWVS*

=cut

          $longest_orf_finder->{ALLOW_5PRIME_PARTIALS} = $part5; #($gffprot =~ /^M/)?0:1;
          $longest_orf_finder->{ALLOW_3PRIME_PARTIALS} = $part3; #($gffprot =~ /\*$/)?0:1;
          # should allow when gene ends near gap or contig end
          
          my $longorf= $longest_orf_finder->get_longest_orf($gffcdna); # == hash
          my($orfprot,$prostart5,$proend3,$prostat)=("",0,0,"");
          
          if(ref($longorf)) {
          $orfprot= $longest_orf_finder->get_peptide_sequence();
          ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();
	        if($part5 or $part3) {
            my $p5= ($orfprot =~ /^M/)?0:1;
            my $p3= ($orfprot =~ /\*$/)?0:1;
            $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
            }
          }
            # even if no orfprot
          my $gffnew= getNewGFF($isoform,$trlen,$orfprot,$prostart5,$proend3);
          foreach my $gf (@$gffnew) {
            if($prostat and $gf =~ /\tmRNA/) {
              $gf =~ s/;protein=/,$prostat;protein=/;
              }
            print $gf,"\n";
            }
          # print join("\n",@$gffnew),"\n";
          next;
      }
      
			# dgg here: compare gff.prot to Longest_orf
			elsif($seq_type eq "checkprot") {
         my $gffprot = $isoform->get_protein_sequence(); # fixme: case of no prot/no CDS, add them
         my $gffcdna=  $isoform->get_cDNA_sequence();
         my $trlen= length($gffcdna) || 1;
         $gffprot ||= "";
         
         $longest_orf_finder->{ALLOW_5PRIME_PARTIALS} = ($gffprot =~ /^M/)?0:1;
         $longest_orf_finder->{ALLOW_3PRIME_PARTIALS} = ($gffprot =~ /\*$/)?0:1;

         my $longorf= $longest_orf_finder->get_longest_orf($gffcdna); # == hash                  
         my $orfprot= "";         
         if(ref $longorf) { $orfprot=$longest_orf_finder->get_peptide_sequence(); }
         
         # watch out for "XXXXXX..." crap, check thru @orfs for longest - XXX ?
         unless( $orfprot ) {
           print "#checkprot: none, $isoform_id $asmbl_id:$model_lend-$model_rend:$orientation\n";
         
      	 } elsif( $orfprot ne $gffprot) { 
			    
			      # check if long contains short ?
			      # if shorter, check if gffprot has inner '*'
			      # need to allow partial orfs
			      
           my ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();
           my ($lenorf, $lengff)= (length($orfprot),  length($gffprot));
           my $diff= $lenorf - $lengff;
           my $nxxx= $orfprot =~ tr/X/X/;
           my $notpart=0;
           if($lengff > 0) {
           if($diff > 0) { $notpart= index($orfprot,substr($gffprot,0,-1)) < 0; }
           else { $notpart= index($gffprot,substr($orfprot,0,-1)) < 0; }
           }

	         ##$lengff||=1;
	         #? allow notpart GOOD sometimes
           my $bad= 0;
           if($lengff > 0) {
           $bad= $notpart || (($diff < 0) && ( 0.5 >= $lenorf/$lengff ));
           }
           $bad= 0 if($bad && index(substr($gffprot,0,-1),'*') >= 0);
           $bad= $bad || ( $diff != 0 and  ($nxxx / abs($diff)) > 0.4 );
           my $tag=($bad)? "POOR": "GOOD";
           $diff="+$diff" if($diff>0);

use constant CHECK2GFF => 1;           
        if(CHECK2GFF) {
          if($bad) { # dont do these..
            print "#checkprot: bad=$diff, $isoform_id $asmbl_id:$model_lend-$model_rend:$orientation\n";
            next;          
          } else {
            print "#checkprot: update_good=$diff, $isoform_id $asmbl_id:$model_lend-$model_rend:$orientation\n";
            my $addattr="checkprot=$diff,update_good";
            my $gffnew= getNewGFF($isoform,$trlen,$orfprot,$prostart5,$proend3,$addattr);
            print join("\n",@$gffnew),"\n";
            next;
          }
        } else {
           # change to gff output? put original CDS + new ?
           $seq= $orfprot;
	         $seq .="\n>OLD.$isoform_id\n$gffprot"; #??
           $com_name .= " CHECKPROT=$diff,$tag,$prostart5-$proend3/$trlen";
           # and continue on
        }
    
      } else {
         print "#checkprot: same, $isoform_id $asmbl_id:$model_lend-$model_rend:$orientation\n";
         next; # no report?
      }
  	}

            $seq =~ s/(\S{60})/$1\n/g; # make fasta format
            chomp $seq;
            
            # my $com_name = $isoform->{com_name} || "";
	    # if ($com_name eq $isoform_id) { $com_name = ""; } # no sense in repeating it

			my $locus = $isoform->{pub_locus};
			my $model_locus = $isoform->{model_pub_locus};
			
			my $locus_string = "";
			if ($model_locus) {
				$locus_string .= $model_locus;
			}
			if ($locus) {
				$locus_string .= " $locus";
			}
			if ($locus_string) {
				$locus_string .= " "; # add spacer
			}
			
            ##dgg.fmt#print ">$isoform_id $gene_id $locus_string $com_name $asmbl_id:$model_lend-$model_rend($orientation)\n$seq\n";
            print ">$isoform_id $gene_id $locus_string $com_name $asmbl_id:$model_lend-$model_rend:$orientation\n$seq\n";
        }
    }
}


exit(0);
#......................................

sub getNewGFF
{
  my($isoform,$trlen,$orfprot,$prostart5,$proend3,$addattr) = @_;
    
  my ($cds5,$cds3)=(0,0);
  my(%CDS, %mRNA);
  my @exons = $isoform->get_exons();
  my $cdna1= 1;
  foreach my $exon (@exons) {
    my ($xend5, $xend3) = $exon->get_coords(); # == get_mRNA_exon_end5_end3()
    $mRNA{$xend5}= $xend3;
                
    my $xd= abs($xend3 - $xend5);
    my $cdna2= $cdna1 + $xd;
    if($cdna1 < $proend3 and $cdna2 > $prostart5) { # overlap
                
      my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
      my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      
      my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
      my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
  
      $CDS{$c5}= $c3;
      $cds5=$c5 if($cdna1 <= $prostart5);
      $cds3=$c3 if($cdna2 >= $proend3);
    }
    
    $cdna1= $cdna2+1;
  }          
  $isoform->populate_gene_obj( \%CDS, \%mRNA,  );
  
  my $aalen=length($orfprot); 
  $aalen-- if(substr($orfprot,-1) eq '*');
  my $clen= $aalen * 3; 
  my $ap=int(100 * $clen/$trlen);
  my $mattr="cxlen=$clen/$trlen;aalen=$aalen,$ap%";
  $mattr.= ";protein=$orfprot" if($orfprot);
  $mattr.= ";$addattr" if($addattr);
  my @gffnew= split(/\n/, $isoform->to_GFF3_format());
  foreach (@gffnew) { if(/\tmRNA/) { s/$/;$mattr/; last; } }
  return \@gffnew;
}

####
sub add_flank {
	my ($seq, $upstream_flank, $downstream_flank, $lend, $rend, $orientation, $genome_seq_ref) = @_;
	
	my $far_left = ($orientation eq '+') ? $lend - $upstream_flank : $lend - $downstream_flank;
	
	if ($far_left < 1) { $far_left = 1; }
	
	my $flank_right = ($orientation eq '+') ? $downstream_flank : $upstream_flank;

	my $left_seq = substr($$genome_seq_ref, $far_left - 1, $lend - $far_left);

	my $right_seq = substr($$genome_seq_ref, $rend, $flank_right);
	
	if ($orientation eq '+') {
		return (lc($left_seq) . uc($seq) . lc($right_seq));
	}
	else {
		return (lc(&reverse_complement($right_seq)) . uc($seq) . lc(&reverse_complement($left_seq)));
	}
}



=item checkprot dgg             

problems: 

-- some gffprot have internal stops? check
-- some shorter orfprot are better, some not

-- crappy prot? from gaps; avoid this

>m7AUGepi4p1s10g69t1 m7AUGepi4p1s10g69t1  /PROT_CHANGE:399-1421 Scaffold10:954530-954754:+
MYKXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXIYYYILNYGNNIKI*
>OLD.m7AUGepi4p1s10g69t1
MKVRDQFDRNKVNTITCDDIVMTFVLENKHGLEDTIIGENVMDDLFACTDDKPMDLDNEM
DPADDDIFSNLFLD*

=item result1 per pred

grep '^>' genes/bestgenes_of7.an7s.checkprot | grep -v '>OLD' | grep _GOOD | perl -pe's/>//; s/ .*$//; s/Gsc.
*//; s/[sp]\d\S+/$1/; s/inity\S+/inity/;  ' | sort | uniq -c | sort -k1,1nr | head -30

Change counts:
   2510 aphid_cuf8r27  / 2623 : almost all bad; fix source CDS
   1196 AUGepi4        / 11801 .. 10% bad for AUG ..
   1174 AUGepir16b     / 9697
    907 AUGepir2       / 6713
    164 AUGepir9
    153 AUGepir10
    126 aphid2Trinity  / 150 : almost all bad
    107 AUGepir3
     86 AUGepi5
      1 ACYPI000049-RA
      1 ACYPI29102-RA

Same counts:
  10228 AUGepi4
    983 AUGepi5
    975 AUGepir10
   8287 AUGepir16b
   5675 AUGepir2
    772 AUGepir3
   1220 AUGepir9
   8140 PASAgasmbl    << all ok
     24 aphid2Trinity
    107 aphid_cuf8r27

grep '^>' genes/bestgenes_of7.an7s.checkprot | grep -v '>OLD' | grep _GOOD | perl -ne'm/>(\S+)/; $d=$1; ($s,$
n)= m/NEWPROT_\w+=(.)(\d+)/;  $d=~s/Gsc.*//; $d=~s/[sp]\d\S+/$1/; $d=~s/inity\S+/inity/; $sum{$d}{$s}+=$n; $cnt{$d}{$s
}++; $ss{$s}++;  END{ @s=sort keys %ss; print join("\t","Class________", map{"N$_"}@s, map{"Av$_"}@s),"\n"; for $d (re
verse sort keys %sum) { for $i (0..$#s) { $s=$s[$i]; $cn[$i]= $c= $cnt{$d}{$s} || 1; $av[$i]= int($sum{$d}{$s}/$c); }
 print join("\t",sprintf("%12s",$d), @cn, @av),"\n"; } }' | head -30

Class________   N+      N-      NAv+    NAv-
aphid_cuf8r27   2504    6       2       30940 < -bogus; Nav+ mostly off by 1
aphid2Trinity   126     1       1       0     < mostly off by +1 
    AUGepir9    137     27      28      98    .. AUG + missing CDSexon parts; NAv- maybe mostly inner stops
    AUGepir3    86      21      20      169   ..
    AUGepir2    619     288     21      180
  AUGepir16b    841     333     19      166
   AUGepir10    129     24      27      91
     AUGepi5    64      22      18      121
     AUGepi4    782     414     19      183
ACYPI29102-RA   1       1       1       0
ACYPI000049-RA  1       1       1       0


# SUMMARY of scores/source, scoretype=homolog:7,paralog:1,ref:3,est:3,pro:3,rseq:2,intr:20,nintron:40,inqual:15,maqual:5,terepeat:-1,UTR:1,CDS:1, 
# Source        N       Ave     UTR     est     homolog inqual  intr    maqual  nintron pro     rseq
# 1_Kept         43620    805   -28.4   391.2   160.7   403.3   782.6   499.3     2.3   421.7   724.6
# 2_Dropped     115071    754   -45.3   590.4   206.7   471.4   1257.1  633.2     3.7   546.4   949.1
# AUGepi4        11801    669    10.2   115.5    88.6    39.2   186.4    52.4     0.3   278.6   144.0
# AUGepi5         1070    956    92.8   992.7   288.1   1003.6  2639.2  1004.0    5.9   693.9   2322.0
# AUGepir10       1142    975    29.2   1055.7  304.5   1102.0  2587.5  1090.4    6.8   726.0   2736.9
# AUGepir16b      9697    623   -64.3   191.2   105.1   129.7   264.3   123.8     0.9   323.6   260.7
# AUGepir2        6713    758   -46.5   193.1   127.3    83.4   330.3    98.7     0.6   403.9   216.4
# AUGepir3         887   1128    60.0   1267.8  384.9   1235.9  3051.6  1220.6    7.6   948.4   2957.1
# AUGepir9        1395    806    10.5   869.8   269.2   906.9   2011.1  938.7     5.8   670.8   2661.0
# acypi              2   3807     2.0   255.5    71.0     0.0   8947.5    0.0     1.5   359.5   146.5
# ars17trinity     150    415     3.1   391.4   192.2   779.3   564.9   571.8     4.3   384.7     0.0
# ars27cuf8       2623    953     7.1   1123.2  438.9   1222.2  2302.8  2144.5    6.8   899.7   1633.1
# pasa2_aphid3    8140    743   -79.2   606.6   188.8   896.4   1194.7  1090.3    4.5   428.6   1191.2

=item cases

  >AUGepir16bp1s10g107t1 AUGepir16bp1s10g107t1   NEWPROT_GOOD=-216,386-526,Scaffold10:1759229-1759369:- Scaffold10:17370
  42-1759754:-
  MICLLRIAPVSVLAGWKRQWRIDGDGVGRQQQQQDNLTLQTTSLSL*
    vv--- crappy prot likely transposon disrupted;
    better model= m6AUGepir3p1s10g115t1, hits some rna-seq but no good rna-asm
    AUG cheats? calling these '*' as X
    > NAPVSVLA is consensus after * ==? transposon fragments
  >OLD.AUGepir16bp1s10g107t1
  MICLLRIAPVSVLAGWKRQWRIDGDGVGRQQQQQDNLTLQTTSLSL*NAPVSVLAGWKRQ
  WRIDGDGVGRQQQQQDNLTLQTTSLSL*NAPVSVLAGWKRQWRIDGDGVGRQQQQQDNLT
  LQTNSLSL*NAPVSVLAGWKRQWRIDGDGVGRQQQQQDNLTLQTNSLSL*NAPVSVLAGW
  KRQWRIDGDGVGRQQQQQDNLTLQTNSLSLWERKDGIKKWLCTKKSCYASILTNGDLNTI
  ELSRHYICRIVGLCPNMLFLVR*

  >AUGepir3p1s10g115t1
	MFHLNAPCRFFAFSTVIVARSPGARAPYRYHDNVMISTTVPMKYFIIV X< IAPVSVLAGWK
	RQWRIDGDGVGRQQQQQDNLTLQTTSLSL X< NAPVSVLAGWKRQWRIDGDGVGRQQQQQDN
	LTLQTTSLSL X< NAPVSVLAGWKRQWRIDGDGVGRQQQQQDNLTLQTNSLSL X< NAPVSVLA
	GWKRQWRIDGDGVGRQQQQQDNLTLQTNSLSL X< NAPVSVLAGWKRQWRIDGDGVGRQQQQ
	QDNLTLQTNSLSL X< NLNKSNNSSDEDNYVSTSFFTWSQKKSENVSVENELSEFFKKGPTK
	KLNVLNSMPTLKKVFIQYNTPLPVSAWPWPLGHKRSLGPSVFLESRTFFYLLVVCLVYRE
	SKDPVKERKDGIKKWLCTKKSCYASILTNGEYVHETINEHHHPENSQQSIERQVLREASK
	YVKYDGDRSFANKESGLRNVWNVLYSPEKLSFWC

>AUGepir16bp2s10g25t1 AUGepir16bp2s10g25t1   NEWPROT_GOOD=-163,257-2956,Scaffold10:2574708-2577407:+ Scaffold10:257445
2-2584157:+
MSKINLQKLIAVLKENGCTPCQWKNVKSGDDIIIKSCKACFNQNTCLKNRLRIRRILVDR
RNDILKLYNQFEPEKKNFTNKVRAVLNELPFEVKRIDSEIVIEICKRLGLVWNADKWKVQ
ALKNHDKVNSALRHIFKQDSKKVKTPVEKVMLYENIGQVSVGDDFLSSSRLMETSITVKE
ITDHGSLLLKPNLEESIVDEFITPTKRNEVIDIVECVNDQLLTHNCNSDIKGKISQVENL
...
>OLD.AUGepir16bp2s10g25t1
MSKINLQKLIAVLKENGCTPCQWKNVKSGDDIIIKSCKACFNQNTCLKNRLRIRRILVDR
RNDILKLYNQFEPEKKNFTNKVRAVLNELPFEVKRIDSEIVIEICKRLGLVWNADKWKVQ
ALKNHDKVNSALRHIFKQDSKKVKTPVEKVMLYENIGQVSVGDDFLSSSRLMETSITVKE
ITDHGSLLLKPNLEESIVDEFITPTKRNEVIDIVECVNDQLLTHNCNSDIKGKISQVENL
..
QVELQIPKRICAIKIKADCDEDNLHLNTPETWRKKGKRPVFLRNNNKILGKSCGENQQN* << inner stop
SVHKNTEKRNYWKDRVLNVICGVPKTWKPKQNVGDIKKWKDMFVQIPYQDDGFNCGVFML
YFANQLMNNKRIINVFNPNSYRLNLQDLILSSSKCMKNICLICGKDEGSFKQNYNCEMDS
TMVQCGSCTRWLHISCLPAMEQKEFENPDWVCGLCSNQCPTI*

=cut


