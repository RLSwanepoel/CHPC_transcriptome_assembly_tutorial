#!/bin/bash
## env bam=xxx trfasta=xxx datad=`pwd` qsub -q shared rnaextest.sh
#PBS -N rnax
#PBS -A ind114
#PBS -l nodes=1:ppn=4,walltime=9:55:00
#PBS -o rnax.$$.out
#PBS -e rnax.$$.err
#PBS -V

ncpu=4
## .. rx-daphmag5xall13f_okall-Dman_84_gsnap took 6hr (4cpu, 178k tr x 33 Gb bam)
## Notes: gsnap reads x trasm is efficient, ~8h x 32cpu per read pair
##  ...   bowtie2 is failing all over, using --all map as needed, outamem and outatime (--fast)
##  ...  ** Should use mRNA.tr oriented by CDS 5-3 ends as express calcs 5/3 end biases **
##
## rnax -p --num_threads is hidden option; dont know if effective
## rnax man:  2 cores always used; not this: ** Maximum of 3 free processor cores required **
## bam list: <lib_1.sam,lib_2.sam,...,lib_N.sam>     A comma-separated list of filenames
## bam must be NAME sorted: samtools sort -n hits.bam 
## .. i have locus-sorted hits.bam ..
## readids: suffix identifiers attached to the end ('/1' and '/2') MUST BE stripped (v1.3.0; 131 drops this)

biobin=$HOME/bio/bin
#bampre='daphmag5xall13f_okall-'
bampre='-daphmag5xall13f'

if [ "X" = "X$datad" ]; then echo "ERR: missing env datad=/path/to/fasta"; exit -1; fi
if [ "X" = "X$bam" ]; then echo "ERR: missing env bam=mydata.bam"; exit -1; fi
if [ "X" = "X$trfasta" ]; then echo "ERR: missing env trfasta=mydata.tr"; exit -1; fi

cd $datad
tnam=`basename $trfasta .tr | sed 's/\.fasta//; s/\.fa$//;'`
bnam=`basename $bam .bam | sed "s/\..*//; s/_1$//; s/$bampre//;"`
## bam=../gsofq4x/Dman_84_merge-daphmag5xall13f.nsort.bam
## bam=../bowtrfq4x1dn/daphmag5xall13f_okall-Dman_70_GCCAAT_L008_1.nsort.bam
## .. or symlink-rename these bams ; esp for list input ..
outrx=rx-$tnam-$bnam
mkdir $outrx

echo "START rnax: `date`"

## add -p $ncpu
echo $biobin/rnaexpress -p $ncpu -o $outrx --no-update-check $trfasta $bam
$biobin/rnaexpress -p $ncpu -o $outrx --no-update-check $trfasta $bam

echo "DONE rnax: `date`"

