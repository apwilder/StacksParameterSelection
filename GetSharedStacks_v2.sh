#!/bin/bash
# This script takes output from StacksDeNovo.sh run in parameter exploration mode and creates summary files that can be plotted () to asses optimal parameter values
 
POPMAP=$1 # path to list of sample IDs and population membership, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/PopulationMap_excl.txt
WKDIR=$2 #working directory with stack output folders, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/stacks2_Outfiles/
SPECIES=$3 #species name to append to output (e.g., SWR, Waterbuck, BFI)

cd $WKDIR
DIRS=`ls -d */` #list directories in main Stacks output directory

###Get stats for varying m #####################
# create header for output file
printf "%s\t" m M n Loci_50pctShared PolymLoci_50pctShared SNPs_50pctShared Loci_80pctShared PolymLoci_80pctShared SNPs_80pctShared > $SPECIES'_SharedPolymorphisms.txt'
printf '\n' >> $SPECIES'_SharedPolymorphisms.txt'

for DIR in $DIRS; do #for each combo of m and M
cd $DIR
m=`echo $DIR | cut -c2` #get m
M=`echo $DIR | cut -c4` #get M
cp ustack.stats ../$SPECIES'_m'$m'M'$M'_ustack.stats' #rename ustack stats files and move to main directory
SBDIRS=`ls -d */` 
for SBDIR in $SBDIRS; do #for each n
cd $SBDIR
`echo pwd`
m=`echo $SBDIR | cut -c2`
M=`echo $SBDIR | cut -c4`
n=`echo $SBDIR | cut -c6- | sed 's"/""g'`

#Output the number of matches to the catalog for each sample
echo SNPs > $SPECIES'_SNPsPerIndiv.txt'
for SAMPLE in `cat $POPMAP | cut -f1`; do
ID=$(echo $SAMPLE | cut -f 1 -d '.')
zcat $ID'.matches.tsv.gz' | cut -f 4,5 | tr '\t' '_' | sort | uniq -c | grep "2 " | wc -l >> $SPECIES'_SNPsPerIndiv.txt'
done

#Get the number of shared loci and SNPs for 50% and 80% of the population
POLYMLOC50=`zcat populations_r0.5.hapstats.tsv.gz | awk '$6>1 {print}' | wc -l`
LOC50=`zcat populations_r0.5.hapstats.tsv.gz | wc -l`
SNPS50=`wc -l populations_r0.5.snps.vcf | cut -f 1 -d' '`
POLYMLOC80=`zcat populations_r0.8.hapstats.tsv.gz | awk '$6>1 {print}' | wc -l`
LOC80=`zcat populations_r0.8.hapstats.tsv.gz | wc -l`
SNPS80=`wc -l populations_r0.8.snps.vcf | cut -f 1 -d' '`

printf "%s\t" $m $M $n $LOC50 $POLYMLOC50 $SNPS50 $LOC80 $POLYMLOC80 $SNPS80 >> $WKDIR$SPECIES'_SharedPolymorphisms.txt'
printf '\n' >> $WKDIR$SPECIES'_SharedPolymorphisms.txt'
cd ..
done
cd ..
done


