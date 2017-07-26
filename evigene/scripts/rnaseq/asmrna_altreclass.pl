#!/usr/bin/perl
# altreclass.pl : evigene alt-reclassifier

## FIXME: *** input aa-size in trclass is wrong, no gap adjust; 
##     original trclass/pubids sorting is right from using aa.qual table **

=item about altreclass.pl

 altreclass.pl : evigene alt-reclassifier
 evigene/scripts/rnaseq/asmrna_altreclass.pl 
 	- follows use of prot/tr2aacds.pl, rnaseq/asmrna_dupfilter, evgmrna2tsa.pl
 	  which create inputs of project/.trclass and project/publicset/.pubids
 	  
  examine evigene mRNA asm publicset and mrna.trclass, to
  identify/remove trival alternates (can be large amount, tends to collect as more trasm added)
  -- drop fragment althi subset: shorter aa, partial/fragment flags
  -- option altrenum: renumber alts so t1 main is longest aa (should have been, not always),
      and t2..tn ordered by alt aasize.  public gene ID is preserved.
  -- option maxaltsame: drop excess complete althi that are same aasize as others 
  
=item usage

  evigene/scripts/rnaseq/asmrna_altreclass.pl -maxalt=19 -altrenum  -trclass kfish2evg367mixx.trclass -out 

 	output: publicset/kfish2evg367mixx.realt_pubids 
  #altreclass: nin/out=499047/499047, ndrop/keep=141000/358047, nrenum=297357, ngenediff=35964/149901
      
=item todo

	merge with evigene asmrna_dupfilter (main classifier)
	offer option to rewrite publicset/ : only need to remove drop sequences (to 2nd files)
		.. but may need logic from evgmrna2tsa to do in compatible way (several tables to update).
	
	- add option for more agressive junk removal?
	-- drop incomplete/frag alts; drop noclass/incomplete|short<60aa 
		(largest fraction of noclass are minimal aa size, with few bases diff.  This is bad
		consequence of using %align criteria, need to add minimal base criteria to it.
		e.g. 5% diff for 1000bp = 50 significant, 5% diff for 100bp = 5 ns.
		
	cat  $pt.realt_pubids | grep -v drop | env noclassmin=60 perl -ne\
	'BEGIN{ $NOCLMIN=$ENV{noclassmin}; } $ok=(/\tmain/ or /complete/)?1:0; @v=split; ($aw)=$v[5]=~m/(\d+)/;
	 $ok=0 if(/noclass/ and $aw<$NOCLMIN); print if($ok);' > $pt.hivalue_pubids
	
	for tp in "aa cds mrna"; do {
 		gunzip -c $pt.${tp}_pub.fa.gz | cat $pt.hivalue_pubids - | env idpre=Pita perl -ne\
 		'BEGIN{ $IDPRE=$ENV{idpre}; } if(/^$IDPRE/) { ($d)=split; $ok{$d}=1; } else { 
 		if(/^>(\S+)/) { $d=$1; $ok=$ok{$d}; } print if($ok); }' > $pt.hivalue.$tp
  }		
		
=item notes  

  require sorted by gene at least, also by aa-size? or do that here?
  	-- input pubids should be gene-sorted
  * revised to skip pre-pubidx table creation .. input .trclass + publicset/.pubics
  
 grep okay *.trclass | cat - publicset/kfish2evg367mixx.pubids | perl -ne \
'if(/\tokay/) { @v=split"\t"; ($od,$cl,$pia,$aq)=@v[0,2,4,5]; \
($piav)= $pia =~ m,^(\d+/\d+),; $oda{$od}="$aq\t$piav\t$cl"; } \
elsif(/^Fun/) { ($pd,$od,$td,$ti)=split; $oda=$oda{$od}||0; s/$/\t$oda/; print; }' |\
 sort -k3,3 -k5,5nr -k1,1 > publicset/kfish2evg367mixx.pubidx2

 cat publicset/kfish2evg367mixx.pubidx2 | env altrenum=0 altreclass.pl > reclass.tab
 # redo for .trclass .pubids input, skip pubid2x

=cut

use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use cdna_evigenesub;  
# use cdna_proteins;
use File::Basename qw(basename dirname fileparse);

use constant VERSION => '2013.08.16';  
my($SAME_PI,$SAME_PA,$DIFF_PA)= (98,97,80); # samecds= $pi>=98 and $pa>=97; diffcds= $pa<=80
my($DROPSHORT_AAPART,$DROPSHORT_AAFULL,$DROPSHORT_ANTISENSE)= (0.9,0.4,0.8); # dropsfull = 0.6 ? 0.5?
my $REVISEPUBFILES=0;
my $arenum=$ENV{altrenum}||0; 
my $maxsame=$ENV{maxsame}||999; # for alts complete, same size, drop excess
my $noclasscut=$ENV{noclasscut}||0; # minaa size for noclass tiny excess,
my $debug=0;
my ($trclass,$pubids,$output,$mapsensetab);

my $optok= GetOptions(
  "class|trclass=s", \$trclass, # require
  "pubids=s", \$pubids, # option
  "output:s", \$output, # option
  "mapsensetab=s", \$mapsensetab, # option, or call opt -mapsensetab ? -gmaptable? or gmap.gff instead?
  #"logfile:s", \$logfile, # option
  "maxsame|maxaltsame=s", \$maxsame, 
  "noclasscut=s", \$noclasscut, 
  "altrenum!", \$arenum, 
  "revisepublicset!", \$REVISEPUBFILES, 
  "debug!", \$debug, 
  );

die "EvidentialGene altreclass VERSION ",VERSION,"
  examine Evigene mRNA publicset, mrna.trclass to reclassify trival alternates

Usage: altreclass.pl -trclass mymrna.trclass > mymrna.reclass.pubids
opts: -pubids publicset/mymrna.pubids  -out mymrna.reclass.pubids 
      -mapsensetab=gmapinfo.table -[no]altrenum  -maxaltsame=$maxsame\n"
  unless($optok and $trclass);  

my $outh=*STDOUT;
unless($pubids) { ($pubids=$trclass) =~ s,(\w+).trclass,publicset/$1.pubids,; } # default evg file set
if(defined $output and not $output) {
	($output=$pubids) =~ s/.pubids//; $output.=".realt_pubids";
}
if($output)	{
  rename($output,"$output.old") if(-f $output);
  open($outh,'>',$output) or die "ERR: writing $output"; 
}

# openloggit($logfile,$trclass);
# loggit(1, "EvidentialGene altreclass.pl VERSION",VERSION);
# loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);
# loggit(0, "altreclass: in trclass=$trclass, pubids=$pubids, out=$output") if($debug);

my($lgd,%trscore,%trv,%trline,%keepdrop,%newid); # globals now for reclassAlts/putgene

my $trinfo= readTrClass($trclass);

my $aaqualhash= readAaQual($trclass,$trinfo);

my $mapinfo= readMapsenseTab($mapsensetab); # if($mapsensetab);

my ($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene)= reclassAlts($pubids,$trinfo,$aaqualhash,$mapinfo);
close($outh);

## need logic from evgmrna2tsa make_pubseq(), make_annotab() to do in compatible way (several tables to update) ..
# revisePublicset() if($REVISEPUBFILES);

warn "#altreclass: nin/out=$nin/$nout, ndrop/keep=$ndrop/$nkeep, nrenum=$nrenum, ngenediff=$ngdiff/$ngene\n";

#---------------------

=item readMapsenseTab

	aligntab3.sh  mrna.gmap.gff > mrna.alnsense.tab
	OR add here code of aligntab3.sh ..

kfish2evg367mixx_realt.alnsense.tab
GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   tag
Scaffold0       5067-10366      -       Funhe2Exy3m032549t1     1-679   678     679     100.0   99.9    0       0       3       6       214,94%,complete        19-663  220     0       gmap
Scaffold0       25351-26047     .       Funhe2Exy3m069279t1     1-697   697     697     100.0   100.0   0       0       1       0       115,49%,complete-utrpoor        295-642 115     0       gmap
Scaffold0       30698-31095     -       Funhe2Exy3m078047t1     1-374   360     374     100.0   95.5    0       1/3     2       0       106,85%,partial3        55-372  45      -1      gmap
     ^^ asense but no valid splice 
Scaffold0       131984-136801   -       Funhe2Exy3m021876t4     1-1261  1254    1261    100.0   99.2    0       0/3     9       18      279,66%,complete        422-1261        270     -1      gmap
     ^^ antisense.valid.ids: hicov, full valid splice/9exons
	 
=cut 

sub readMapsenseTab {
	my($maptab)= @_;
	return {} unless($maptab);
	my($ok,$hin)= openRead($maptab);  die "ERR: reading $maptab" unless($ok);
	my %maps=();
	## GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   tag
	my @hd; my @hset=qw(QueryID qlen cov nexon splice aalen sense);
	while(<$hin>) {
		next unless(/^\w/);
		chomp; my @v=split"\t";		
		my($td,$ql,$cov,$nx,$nspl,$aw,$sens);
		if(/\tcov/ and /\tsense/ and /\tsplice/) {
			@hd=@v; map{ s/AQuery/Query/; } @hd; 
			my %hset=map{ $_=>1 } @hset; for my $h (@hd) { $hset{$h}=2 if($hset{$h}); } # check fields
			my @miss= grep{ $hset{$_} == 1 } @hset;
			if(@miss) { die "ERR: missing Mapsense table fields: @miss\n"; }
			next;
		} elsif(@hd) {
			my %v=(); for my $i (0..$#v) { my $h=$hd[$i]; $v{$h}=$v[$i]; }
			($td,$ql,$cov,$nx,$nspl,$aw,$sens)= @v{@hset}; # okay??
		} else {
			($td,$ql,$cov,$nx,$nspl,$aw,$sens)=@v[3,6,7,11,12,13,16]; # or require this?
		}
		my $nsx=($nx<2)?0:$nspl/$nx; 
		my $as1=($sens<0 and $nx>3 and $nsx> 1.5 and $cov>90)?1:0; # certain?
		my $as2=($sens<0 and not $as1 and (($nsx >= 1.5 and $cov>95) or ($nsx >= 1.8 and $nx > 5 and $cov>85)))?1:0; # likely
		my $asense=($as1 or $as2)?"antisense":"ok";
		my $quals="cov=$cov,nexon=$nx,splicex=$nsx,sense=$asense";
		$maps{$td}= $quals;
		} close($hin);
	return \%maps;
}

## FIXME: *** input aa-size in trclass is wrong, no gap adjust; 
##     original trclass/pubids sorting is right from using aa.qual table **

sub readTrClass {
	my($trclass)= @_;
	my($ok,$hin)= openRead($trclass);  die "ERR: reading $trclass" unless($ok);
	my %oda=();
	while(<$hin>) {
		next unless(/^\w/ and /\tokay/);
		chomp; my @v=split"\t"; 
		## trclass cols:
		##  oid,okay/drop,class,idbestmatch,pIdAln,aaqual,aaref/flags
		my($od,$cl,$pia,$aq,$aaref)=@v[0,2,4,5,6];  
		
		## FIXME: add aaref = aablastp ref from flag column 6, if there, as new qual score.
		## asmrna_dupfilter2.pl -ablast blastp.tall4, adds bitscore,refid, and refbest/refok/refgood flag
		## empty col6 == 0,0,pflag:0  [pflag == poor if > 0]
		## aaref col6 == 250,arath:AT4G28760.2,refok,pflag:0
		## 164.4,arath:AT3G42860.1,pflag:0
		## 224,arath:AT1G69530.4,aadup:1AB-I11_VO_k30Loc10139t3,pflag:0
		if($aaref =~ /^0,0/) { $aaref="0"; }
		else { $aaref =~ s/,pflag:.*//; $aaref=~s/,aadup:.*//;  } # trim other like ,aadup ? $aaref="aaref:$aaref";
		
		# my($piav)= $pia =~ m,^(\d+/\d+),;  #? full form: pi/pa[/optbestmatchid]
		$oda{$od}="$aq\t$pia\t$cl\t$aaref"; 
		} close($hin);
	return \%oda;
}

# readAaQual: inputset/*.aa.qual, or use evigenesub:fasize_nogaps(publicset/*.aa_pub.fa.gz,okc='A-WYZa-wyz',gapc='Xx*')
sub readAaQual { 
	my($trclass,$trinfo)= @_;
	my %aaqual=();

  my($name,$path,$suffix) = fileparse($trclass,qr/\.\w*/); # suf-'.trclass' or suf = qr/\.\w*/
  my $aaqualf = "$path/inputset/$name.aa.qual"; 
  $aaqualf = "$path/publicset/$name.aa.qual" unless(-f $aaqualf);

  ## this is a mess, should update aa.qual fields in .trclass, others, to include nnnX gaps and reduce? aasize
  if(-f $aaqualf) {
    my($ok,$hin)= openRead($aaqualf);  
    while(<$hin>) {
      next unless(/^\w/);
      my($od,$awnogap,$ngap,$aaq,$clen)=split;
      my($aw)= $aaq=~/(\d+)/;
      if($awnogap < $aw) {
        $aaq =~ s/$aw/$awnogap/; $aaq.=",${ngap}X";
          # this isn't effective for publicset/ ids when trinfo has orig ids.. fixme
        if(my $tri= $trinfo->{$od}) { ## change trinfo.aq value?
          $tri =~ s/^\S+/$aaq/; $trinfo->{$od}= $tri;
        }
      }
      $aaqual{$od}=$aaq; #? all or just awnogap?
    } close($hin);
  }  else {
    my $aaseq = "$path/publicset/$name.aa_pub.fa"; 
    $aaseq .= ".gz" unless(-f $aaseq); 
	  if(-f $aaseq) {
	    my $fasizes= fasizes_nogap($aaseq,'amino');
	    foreach my $od (sort keys %$fasizes) {
        my($awnogap,$aw,$ngap)=split "\t",$fasizes->{$od};
        if($awnogap < $aw) {
          # this isn't effective for publicset/ ids when trinfo has orig ids.. fixme
          if(my $tri= $trinfo->{$od}) { ## change trinfo.aq value?
            $tri =~ s/^\d+/$awnogap/; 
            $trinfo->{$od}= $tri;
          }
        }
        $aaqual{$od}=$awnogap; #? all or just awnogap?
	    }
	  }
	}
	return \%aaqual;
}

=item revisePublicset

	need logic from evgmrna2tsa to do in compatible way (several tables to update) ..
	revise publicset/$trname.{pubids,mrna_pub.fa,cds_pub.fa,aa_pub.fa,ann.txt} ..
	also revise .trclass okay/drop 
	# my @keepids= sort keys %{$keepdrop{'keep'}}; # new pd, od, old pd ..
	# my @dropids= sort keys %{$keepdrop{'drop'}}; # new pd, od, old pd ..

	## evgmrna2tsa2.pl	
	my($pubmrna,$npm)	= make_pubseq($cdnaseq,'mRNA',\%annothash);
	my($pubaa,$npa) 	= make_pubseq(makename($cdnaseq,'.aa'),'protein',\%annothash);
	my($pubcds,$npc)	= make_pubseq(makename($cdnaseq,'.cds'),'CDS',\%annothash);

=cut

sub revisePublicset 
{

# 	my($pubd,@ft)= getFileset('publicset','mrna_pub.fa|cds_pub.fa|aa_pub.fa');   
# 	foreach my $inseq (@ft) {
# 		my($drop,$nin,$fkeep,$fchange,$fdrop)=(0,0);
# 		my($ok,$inh)= openRead($inseq);
# 		(my $outfa=$inseq) =~ s/_pub/_pubrealt/;  #realt_pubids
# 		$ok= open($outh,'>',$outfa) if($ok);
# 		while(<$inh>) { 
# 			if(/^>(\S+)/) {  my $d=$1; # old pubid
# 				$drop= ($keepdrop{'drop'}{$d}) ? 1 : 0;
# 				my $nd= $newid{$d}; if($nd and $nd ne $d and not $drop) { s/>$d/>$nd/; $fchange++;}
# 				$nin++; $fkeep++ unless($drop); $fdrop++ if($drop);
# 			} 
# 			print $outh $_ unless($drop);
# 		}
#   	close($inh); close($outh);
# 		#loggit(0,"revised $inseq to $outfa: nin=$nin, nkeep=$fkeep, nchangeid=$fchange"); 
# 	}
	
# 	my($annotab, $ntra1, $ncdsa1) 
# 		= make_annotab($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
# 	loggit(0,"revised publicset: ",$pubmrna,$pubaa,$pubcds,$annotab); 

}


sub reclassAlts {
	my($pubids,$trinfo,$aaqual,$mapinfo)= @_;
	my($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene)=(0) x 10;
	my($ok,$pubidh)= openRead($pubids); die "ERR: reading $pubids" unless($ok);
	while(<$pubidh>) {
		unless(/^\w/) {  # BUT print header , adding new col names
		#pubids.hdr #Public_mRNA_ID originalID      PublicGeneID    AltNum
		if(/^#Pub/) { s/$/\tClass\tAAqual\tpIdAln\tNotes/; print $outh $_; }
		next; } 
		chomp; 
		my($pd,$od,$gd,$ti)=split"\t"; $nin++; # ,$aq,$pia,$cla
		
		unless($gd eq $lgd) { 
			if($lgd) { my($aout,$gdiff,$anum,$akeep,$adrop)= putgene(); 
				$ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; $ndrop+=$adrop; $nkeep+=$akeep; }
			%trscore=%trv=%trline= ();
		}
		
		my $trinfo= $trinfo->{$od}||"0";
		my($aq,$pia,$cla,$aaref)= split"\t", $trinfo;
		## fixme: add $aaref from trinfo, often 0/empty
		my($aw,$ap,$ac)=split",",$aq;
		my($pi,$pa)=split"/",$pia; 

		my $mapinfo= $mapinfo->{$pd}||$mapinfo->{$od}||"0"; # antisense info; which id here? pd or new realt pd or od ?
				## mapinfo == "cov=$cov,nexon=$nx,splicex=$nsx,sense=$asense" : keep any more than antisense?
		my $antisense=($mapinfo =~ /antisense/)?-1:0; # or -1:0 or "antisense" ?
		#	 $trline .= ",antisense" if($antisense<0); # ?? trinfo == "$aq\t$pia\t$cl" >> class,anti ? or aaqual,anti ?
	  if($antisense<0) {
	  	$trinfo=join("\t","$aq,antisense",$pia,$cla,$aaref); # more sensible to add to aaqual field
	  }
		my $trline="$_\t$trinfo"; chomp($trline); # append trinfo !!
		
		my $samecds= ($pi>=$SAME_PI and $pa>=$SAME_PA)?1:0;
		my $diffcds= ($pa<=$DIFF_PA)?1:0; # skip ident score here? alt is diff enough when align score low
		   $samecds=0 if($diffcds); # shouldnt need..
		   
		my $icla= ($ac =~ /complete/)?3:($ac =~ /partial[35]/)?2:1; 
		$icla-- if($cla =~ /frag/);
		
## FIXME: *** $aw aa-size from trclass is wrong, no gap adjust; 
##     original trclass/pubids sorting is right from using aa.qual table **
## FIXME2: add $aaref : must-keep if $aaref =~ /ref(best|good|ok)/, maybe-keep if $aaref=~ bitscore,refid
##     add to trv{pd} BUT, dont change trscore sort order for must-keep.

		#BAD# my $trscore= $aw;
    my $aaqual= $aaqual->{$pd} || $aaqual->{$od} || "";
    my ($aasize)= $aaqual=~m/^(\d+)/;
    if($aasize>0 and $aasize<$aw) { 
      $aw=$aasize;
    } # fixme: update trinfo/trline ??
    my $trscore= $aasize || 99999 - $ti;

		$trscore = 1 + int($trscore*0.75) if($antisense<0); # antisense reduces score to lower value in putgene
		$trscore{$pd}= $trscore; ## add icla to score? or $aw * (97+$icla)/100
		
		$trv{$pd}=join "\t", $aw,$cla,$icla,$samecds,$diffcds,$antisense,$aaref; ## add ,cla ***
		$trline{$pd}= $trline;
		$lgd=$gd;  
	} close($pubidh);
	
	if($lgd) { my($aout,$gdiff,$anum,$akeep,$adrop)= putgene(); 
		$ngene++; $ngdiff++ if($gdiff); $nrenum+=$anum; $nout+=$aout; $ndrop+=$adrop; $nkeep+=$akeep; }
	## add summary counts:  ndrop, nrenum, ngenechanged, nin, nout
	return($nin,$nout,$ndrop,$nkeep,$nrenum,$ngdiff,$ngene);
}

sub putgene {
  my ($changed,$renum)=(0,0); 
  my @tr= sort{ $trscore{$b} <=> $trscore{$a} or $a cmp $b } keys %trscore;
  my $t1= shift @tr or return 0; # main/longest tr

	#	$trv{$pd}=join"\t",$aw,$cla,$icla,$samecds,$diffcds,$antisense,$aaref; ## add ,cla ***
  my($t1aw,$t1cla,$t1icl,$ti1same,$t1diff,$t1anti,$t1aaref)= split"\t", $trv{$t1};
  ## if t1 is antisense, need to change below drops, add antisense -trscore ??
  $t1aw= $trscore{$t1} if($t1anti<0); ## reduce t1aw, use trscore?
  
  my (@drops,@keeps);
  my $ialt=1; $t1aw||=1;  my ($ltaw)=($t1aw);

	my $AABITS_MIN_NOCL = 69; # or what?
	my $AABITS_MIN_ALT  = 99; # or what?
	
  foreach my $tr ($t1,@tr) {
    my ($drop,$keep)=(0,0); 
    
    if($tr eq $t1) { # check for noclass drop option
    	## NO: t1icl == number not trclass
      if($noclasscut>0 and $t1cla =~ /noclass/) {  # check for noclass drop option
        $drop=($t1aw < $noclasscut)?1:0;
        if($t1aaref=~/ref(best|good|ok)/) { $keep=1; $drop=0; } ##?? PROBLEM FAILING TEST ** why? trv:\t,
        elsif($t1aaref=~/^(\d+)/) { my $bits=$1; $drop=0 if($bits>$AABITS_MIN_NOCL); } # what?
        
      } else { # should be main, but maybe not ..
        $keep=1;
      }
      
    } else { # if($tr ne $t1)
      my ($taw,$tcla,$ticl,$tisame,$tidiff,$tianti,$tiaaref)= split"\t", $trv{$tr};
   		$taw= $trscore{$tr} if($tianti<0); ## reduce taw/use trscore?
      my $paw=$taw/$t1aw;
      
      $drop=(($ticl< $t1icl and $paw<$DROPSHORT_AAPART and not $tidiff) 
      	  or ($ticl>=$t1icl and $tisame and $paw < $DROPSHORT_AAFULL)) ? 1 :0;
      ## ^^ need to check if dropping all slightly short partials is good idea; 

			$drop=1 if($tianti < 0 and $paw < $DROPSHORT_ANTISENSE); ## antisense check ..
			
			## aablastp ref score keeps some drops, doesn't change alt order
			## HOWEVER, may have several alts same aaref score.. modify maxsame check? add trhoscore?
			if($tiaaref) { ## and $drop < NOT this
			 if($tiaaref=~/ref(best|good|ok)/) { $keep=1; $drop=0; }
			 elsif($tiaaref=~/^(\d+)/) { my $bits=$1; $drop=0 if($bits>$AABITS_MIN_ALT); }
      }
      
      if($ialt > $maxsame and not $drop and not $keep) { $drop=1 if($taw > $ltaw-12); } # $taw == $ltaw or 
      ## ^ this classing need checking. may be dropping valid alts
      
      $ltaw=$taw; # unless drop ??
    }
    
    $drop=0 if $keep; # fix what?
    ## new col for aaref? missing often; append to Notes col
    my $trline=$trline{$tr}; 
    if($drop) {
      push @drops, $trline;
    } else { 
    	push @keeps, $trline;
    	$ialt++; # need to count good alts here for maxsame test. reset for output
    }
    
#     my($pd,$od,$gd,$ti,$aq,$pia,$cla,$aaref)=split"\t",$trline; ## ,antisense now on $cla : move to add attr?    
#     my $notes=($aaref)?"aaref:$aaref,":"";
#     if($drop) {
#       push @drops, join("\t",$pd,$od,$gd,$ti,$aq,$pia,'drop'.$cla,$notes)."\n";  
#     } else { 
#     	# collect @keeps, as @drops, recheck for dupl prots >> maxsame, then print
#     	push @keeps, join("\t",$pd,$od,$gd,$ti,$aq,$pia,$cla,$notes)."\n";  
#     	$ialt++; # need to count good alts here for maxsame test. reset for output
#     }
    
  }
	
	# maybe recheck @keeps here for excess same alts
	# save drop/keep pd,od,oldpd in hash for other uses, ie rewrite publicset files

	## UPDATE unless($arenum) put newID in notes column, and maybe put ialt into ti column?
	## FIXME: change classes: 1st class to 'main', change 'dropmain' to 'dropalt', main>alt if not 1st
	##  .. and preserve oldclass:$cla in Notes 
  ## BAD Output :  ^,newid:PitaTv1R000040t115  n=110585; CHOMP bug ??  on trline?? 
	
	$ialt=0; # reset for output.
  foreach my $trline (@keeps) {
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes)=split"\t",$trline; 
   	$ialt++; my $oldid=$pd; 
   	if($notes) { $notes="aaref:$notes" if($notes =~ /^\d/); $notes.=","; }
   	else { $notes=""; } # not "0"
   	
   	## if( $aq =~ s/,antisense// ) { $notes.="antisense,"; } ## leave on aaqual or class field
   	if($ti != $ialt) { my $newd=$gd.'t'.$ialt; 
  		if($arenum) { $notes.="oldid:$pd,"; $pd=$newd; $ti=$ialt; $renum++; $changed++; }
   		else { $notes.="newid:$newd,";  $ti=$ialt; } # change ti or not?
   	}
    print $outh join("\t",$pd,$od,$gd,$ti,$cla,$aq,$pia,$notes||'.')."\n"; 
    my @ids=($pd,$od); if($oldid ne $pd) { push @ids, $oldid; $newid{$oldid}=$pd; } $keepdrop{'keep'}{@ids}=$trline; #?
  }
  
  foreach my $trline (@drops) {
    my($pd,$od,$gd,$ti,$aq,$pia,$cla,$notes)=split"\t",$trline; 
    $ialt++; my $oldid=$pd;
   	if($notes) { $notes="aaref:$notes" if($notes =~ /^\d/); $notes.=","; }
   	else { $notes=""; } # not "0"
   	## if( $aq =~ s/,antisense// ) { $notes.="antisense,"; } ## leave on aaqual or class field
   	if($ti != $ialt) { my $newd=$gd.'t'.$ialt; 
  		if($arenum) { $notes.="oldid:$pd,"; $pd=$newd;  $ti=$ialt; $renum++; }
   		else { $notes.="newid:$newd,";  $ti=$ialt; } # change ti or not?
   	}
    print $outh join("\t",$pd,$od,$gd,$ti,'drop'.$cla,$aq,$pia,$notes||'.')."\n"; 
    my @ids=($pd,$od); if($oldid ne $pd) { push @ids, $oldid; $newid{$oldid}=$pd; } $keepdrop{'drop'}{@ids}=$trline; #?
    $changed++;
  }
  return ($ialt,$changed,$renum, scalar(@keeps),scalar(@drops));
}

__END__

evigene/scripts/rnaseq/asmrna_altreclass.pl -maxalt=19 -altrenum  -trclass kfish2evg367mixx.trclass -out 
#altreclass: nin/out=499047/499047, ndrop/keep=141000/358047, nrenum=297357, ngenediff=35964/149901

## adding gmapsense..
pt=kfish2evg367mixx
$evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $pt.trclass -mapsensetab publicset/$pt.gmapsense.tab \
  -out $pt.asenrealt2.idtab -maxalt=19 -altrenum -debug
#altreclass: nin/out=499047/499047, ndrop/keep=110027/389020, nrenum=299194, ngenediff=34890/149901
  ^^ more are kept here, why?
cat kfish2evg367mixx.asenrealt2.idtab | grep antisens | wc -l =    7552
	drop  = 4502; keep = 3050
	dropmain = 340 << these need to be checked w/ blastp for bad antisense calls.
	Funhe2Exx3m006258t29    Fungr1EG3m002911t1      Funhe2Exx3m006258       29      dropmain        747,78%,complete,antisense      99/100  oldid:Funhe2Exx3m006258t1
	Funhe2Exx3m010025t7     Fungr1EG3m004981t1      Funhe2Exx3m010025       7       dropmain        573,76%,complete,antisense      99/100  oldid:Funhe2Exx3m010025t1
		^^ antisense yes, but revaa/tr may be best model, should recalc fwd aa,cds, not drop
		
cat kfish2evg367mixx.asenrealt2.idtab | grep antisens | grep -v drop | head
Funhe2Exx3m001262t1     Funhe2Emap3m000720t1    Funhe2Exx3m001262       1       main    1502,89%,complete,antisense     99/99   .
Funhe2Exx3m006254t1     Funhe2Emap3m004077t1    Funhe2Exx3m006254       1       main    748,21%,complete-utrbad,antisense       99/100  .
Funhe2Exx3m006448t3     Fungr1EG3m003035t1      Funhe2Exx3m006448       3       main    734,62%,complete,antisense      98/97   oldid:Funhe2Exx3m006448t1
Funhe2Exx3m006448t4     Funhe2Eq7m051374t1      Funhe2Exx3m006448       4       altmid  700,63%,complete,antisense      98/97   .
		^^ 448t3,4 : was t1,t2; old t3,4 became new t1,t2 due to antisense.
Funhe2Exx3m007789t3     Funhe2Emap3m005117t1    Funhe2Exx3m007789       3       main    662,51%,complete-utrpoor,antisense      99/100  oldid:Funhe2Exx3m007789t1

grep Funhe2Exx3m006448t kfish2evg367mixx.asenrealt2.idtab
** FIXME: main reset here, need change class output..
## human:UniRef50_Q9H0H3 Ectoderm-neural cortex protein 2, 589aa matches new 448t2, t1 may be bogus partial5.
Funhe2Exx3m006448t1     Funhe2E6bm006723t2      Funhe2Exx3m006448       1       altmid  654,84%,partial5        98/91/Funhe2Eq7m051374t1        oldid:Funhe2Exx3m006448t3
Funhe2Exx3m006448t2     Funhe2E6bm006723t3      Funhe2Exx3m006448       2       althi1  589,82%,complete        98/100/Funhe2Eq7m051374t1       .
Funhe2Exx3m006448t3     Fungr1EG3m003035t1      Funhe2Exx3m006448       3       main    734,62%,complete,antisense      98/97   oldid:Funhe2Exx3m006448t1
Funhe2Exx3m006448t4     Funhe2Eq7m051374t1      Funhe2Exx3m006448       4       altmid  700,63%,complete,antisense      98/97   .


eg. keep partial?  t25 part3 here is bigger than t6+, 99/92 says 8% differs from old main
.. map view says likely good form, almost complete
Funhe2Exx3m001853t25    Funhe2Eq7m065601t1      Funhe2Exx3m001853       25      1264,95%,partial3       99/92   dropalthi1      oldid:Funhe2Exx3m001853t10
..
Funhe2Exx3m001853t1     Funhe2E6bm001704t3      Funhe2Exx3m001853       1       1513,95%,complete       99/87   althi   oldid:Funhe2Exx3m001853t3
Funhe2Exx3m001853t2     Funhe2E6bm001704t6      Funhe2Exx3m001853       2       1286,89%,complete       99/90   althi   oldid:Funhe2Exx3m001853t6
Funhe2Exx3m001853t3     Funhe2E6bm001704t37     Funhe2Exx3m001853       3       1279,90%,complete       99/90   althi1  oldid:Funhe2Exx3m001853t48
Funhe2Exx3m001853t4     Funhe2E6bm001704t40     Funhe2Exx3m001853       4       1277,91%,complete       99/91   althi1  oldid:Funhe2Exx3m001853t41
Funhe2Exx3m001853t5     Funhe2E6bm001704t35     Funhe2Exx3m001853       5       1270,89%,complete       99/92   althi1  oldid:Funhe2Exx3m001853t42
Funhe2Exx3m001853t6     Funhe2E6bm001704t14     Funhe2Exx3m001853       6       1116,83%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t25
Funhe2Exx3m001853t7     Funhe2E6bm001704t28     Funhe2Exx3m001853       7       1116,90%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t37
Funhe2Exx3m001853t8     Funhe2E6bm001704t38     Funhe2Exx3m001853       8       1116,87%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t45
Funhe2Exx3m001853t9     Funhe2E6bm001704t29     Funhe2Exx3m001853       9       1108,93%,complete       99/94   althi1  oldid:Funhe2Exx3m001853t31
Funhe2Exx3m001853t10    Funhe2E6bm001704t25     Funhe2Exx3m001853       10      997,79%,complete        99/94   althi1  oldid:Funhe2Exx3m001853t26
...
Funhe2Exx3m001853t20    Funhe2E6bm001704t1      Funhe2Exx3m001853       20      1292,99%,partial        99/91   dropmaina2      oldid:Funhe2Exx3m001853t1
 > case of dropped main .. is this right? partial but 2nd longest. prob same as Funhe2Exx3m001853t6 = 1286aa complete + 6aa; maybe bad orf call


/bio/bio-grid/kfish2/rnas/kf2evgr/trevg367mixx
poor alt removal addition to trclassing 
may need to do after making pubids ?
separate altpoor classifier?

#......
gmapz/altpoordrop.info

## redo / add to trclass to remove many poor/partial alts,
   use pubids + trclass info: aaqual, class, maybe also p-id/aln score?
   sort per gene by aasize, qual?
   per gene remove (to other file) cases where alt.aasize << main.aasize, alt class=frag, alt.aaqual = partial
   especially for common case of % id/aln ~ 99/100 (ie subset alt)
 
# redo, add $pia to pubidx 
grep okay *.trclass | cat - publicset/kfish2evg367mixx.pubids | perl -ne\
'if(/\tokay/) { @v=split"\t"; ($od,$cl,$pia,$aq)=@v[0,2,4,5]; \
($piav)= $pia =~ m,^(\d+/\d+),; $oda{$od}="$aq\t$piav\t$cl"; } \
elsif(/^Fun/) { ($pd,$od,$td,$ti)=split; $oda=$oda{$od}||0; s/$/\t$oda/; print; }' |\
sort -k3,3 -k5,5nr -k1,1 > publicset/kfish2evg367mixx.pubidx2

# alt-reclassifier: require sorted by gene at least, also by aa-size? or do that here?
# .. keep orig pubid or renum alts? need option
cat publicset/kfish2evg367mixx.pubidx2 | env altrenum=0 perl altreclass.pl

