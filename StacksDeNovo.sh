#!/bin/bash

#This program builds de novo stacks from RADseq fastq files, and can be used to find optimal parameter settings for m, M and n when run in exploratory mode

#There are two ways to run this program: 
# 1. explore a range of each parameter value of m, M and n, holding other parameter values constant 
# 2. specify m, M and n and the program will run only those values

#Define variables
R1DIR=$1 #directory of R1 reads, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/ProcessRADtagsOutput/ConcatenatedR1/
R2DIR=$2 #directory of R2 reads, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/ProcessRADtagsOutput/ConcatenatedR2/
OUTPUTMAIN=$3 #directory for output, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/stacks2_Outfiles/
POPMAP=$4 #path to popmap file, e.g. /media/jamie/4TB/TAM_run01/Waterbuck_January2018/PopulationMap_excl.txt (tab-delimited table of sample ID, Population ID)
m1=$5 #Range of m (lower bound) [or m if running specific parameter values]
m2=$6 #Range of m (upper bound) [or M if running specific parameter values]
M1=$7 #Range of M (lower bound) [or n if running specific parameter values]
M2=$8 #Range of M (upper bound) [or nothing if running specific parameter values]
n1=$9 #Range of n (lower bound) [or nothing if running specific parameter values]
n2=${10} #Range of n (upper bound) [or nothing if running specific parameter values]
HOLDm=${11} #The default value to hold m when varying other parameters [or nothing if running specific parameter values]
HOLDM=${12} #The default value to hold M when varying other parameters [or nothing if running specific parameter values]
HOLDn=${13} #The default value to hold n when varying other parameters [or nothing if running specific parameter values]

##### Example of call to explore range of parameter values:
# nohup ./Stacksv2.sh /media/jamie/4TB/TAM_run01/SWRhino_November2017/ProcessRADtagsOutput/R1/R1cat/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/ProcessRADtagsOutput/R2/R2cat/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/stacks_denovo_output/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/PopulationMap.txt \
# 1 3 1 3 1 10 3 2 0 >& Stacks2_SWR_m1-3_M1-3_n1-10_m3M2n0.nohup &

##### Example of call to run specified parameter values:
# nohup ./Stacksv2.sh /media/jamie/4TB/TAM_run01/SWRhino_November2017/ProcessRADtagsOutput/R1/R1cat/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/ProcessRADtagsOutput/R2/R2cat/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/stacks_denovo_output/ \
# /media/jamie/4TB/TAM_run01/SWRhino_November2017/PopulationMap.txt \
# 1 3 1 3 1 10 3 2 0 >& Stacks2_SWR_m1-3_M1-3_n1-10_m3M2n0.nohup &

##########################MAIN FUNCTION####################################################
StacksPipe () {
OUTPUTDIR=$OUTPUTMAIN'm'$m'M'$M'/' #define output directory with this combo of M and m to put in main output folder
if [ -d $OUTPUTDIR ] #if directory already exists, then skip ustacks for this combo of M and m (these parameters have already been done)
	then
	echo 'Directory '$OUTPUTDIR' exists!\n'
	cd $OUTPUTDIR

else #otherwise create the directory
	mkdir $OUTPUTDIR

	###Align reads into stacks of identical sequences using ustacks#############
	cd $OUTPUTDIR
	###create header for summary output file for ustacks (compiles info for all samples)
	echo -e 'sampleID\treads\tstacks\treads_in_stacks\t2ndry_stacks\treads_in_2ndry_stacks\tstack_cvg_mean\tstack_cvg_stdev\tstack_coverage_max\trawreads_in_stacks\trpt_stacks\tcvg_after_rpt_rmvl_mean\tcvg_after_rpt_rmvl_stdev\tcvg_after_rpt_rmvl_max\treads_in_stacks_after_rpt_rmvl\tnon-rpt_stacks\tassmbld_stacks\tnot_assmbld_stacks\tcvg_assmbld_mean\tcvg_assmbld_stdev\tcvg_assmbld_max\treads_in_stacks_assmbld\tfinal_cvg_mean\tfinal_cvg_stdev\tfinal_cvg_max\tfinal_reads_in_stacks' > $OUTPUTDIR'ustack.stats'
	#loop through all samples in PopMap
	for SAMPLE in `cat $POPMAP | cut -f1`; do 
	ID=$(echo $SAMPLE | cut -f 1 -d '.') #define sample ID (if the sample is the first entry of a dot-delimited string)

	###disable use of secondary reads if m>6, as recommended by Paris et al. 2017
	if [ $m -gt 6 ] #if >m6 disable use of secondary reads, should be > for linux!! -gt for mac
	then
	ustacks -f $R1DIR$SAMPLE'.fq' -i $ID -o $OUTPUTDIR -M $M -m $m -p 12 -t fastq -d -H --model_type snp >& $OUTPUTDIR$ID'_ustack.out'
	else #otherwise use secondary reads
	ustacks -f $R1DIR$SAMPLE$'.fq' -i $ID -o $OUTPUTDIR -M $M -m $m -p 12 -t fastq -d --model_type snp >& $OUTPUTDIR$ID'_ustack.out'
	fi
	###get info from *_ustack.out and write it to a summary file
	printf "%s\t" $ID $(grep " reads; formed:" $ID'_ustack.out' | cut -f2 -d " ") $(grep "primary reads" $ID'_ustack.out' | tr -s ' ' | cut -f 2,5 -d' ' | tr ' ' '\t') \
	$(grep "secondary stacks representing " $ID'_ustack.out' | tr -s ' ' | cut -f 2,6 -d' ' | tr ' ' '\t') \
	$(grep "Stack coverage: " $ID'_ustack.out' | tr '=;(' '\t' | cut -f2,4,6,8) $(grep Blacklisted $ID'_ustack.out' | cut -d' ' -f2) \
	$(grep "repeat removal: " $ID'_ustack.out' | tr '=;(' '\t' | cut -f 2,4,6,8) $(grep Assembled $ID'_ustack.out' | tr '; ' '\t' | cut -f2,5,8) \
	$(grep "after assembling stacks:" $ID'_ustack.out' | tr '=;(' '\t' | cut -f 2,4,6,8) $(grep "Final coverage:" $ID'_ustack.out' | tr '=;(' '\t' | cut -f 2,4,6,8) >> $OUTPUTDIR'ustack.stats'
	printf '\n' >> $OUTPUTDIR'ustack.stats'
	done
fi

###Function to run through n-dependent part of stacks ######################
nPipe () {
if [ -d $OUTPUTDIR'm'$m'M'$M'n'$n ] #if directory already exists, skip this step (these parameters have already been done)
then
echo 'Directory '$OUTPUTDIR'm'$m'M'$M'n'$n' exists!\n'
cd $OUTPUTDIR
else #otherwise create the directory

###Build a catalog from consensus loci using cstacks#############
cstacks -P $OUTPUTDIR -M $POPMAP -b $n -n $n -p 12 >& $OUTPUTDIR'Cstack.out'

#compile cstacks statsfile showing how many SNPs are added to the catalog per sample
echo NA > $OUTPUTDIR'cstack_1.stats'
grep -w "loci were matched to a catalog locus." Cstack.out | tr -s ' ' | cut -f 2 -d' ' >> $OUTPUTDIR'cstack_1.stats'
grep -w "loci were newly added to the catalog." Cstack.out | tr -s ' ' | cut -f 2 -d' ' > $OUTPUTDIR'cstack_2.stats'
echo NA > $OUTPUTDIR'cstack_3.stats'
grep -w "loci matched more than one catalog locus and were excluded." Cstack.out | tr -s ' ' | cut -f 2 -d' ' >> $OUTPUTDIR'cstack_3.stats'
paste cstack_*.stats > $OUTPUTDIR'cstack.stats'
sed -i '1 i\Loci_matched\tLoci_added\tLoci_excluded' $OUTPUTDIR'cstack.stats'
rm $OUTPUTDIR'cstack_'*'.stats'

###Match all samples against the catalog using sstacks#############
sstacks -P $OUTPUTDIR -M $POPMAP -p 12 -b $n >& $OUTPUTDIR'Sstack.out' 

tsv2bam -P $OUTPUTDIR -M $POPMAP -t 12 -R $R2DIR -b $n #convert to tsv files to bam

###Merge paired-end reads to contigs, align reads and call SNPs #############
gstacks -P $OUTPUTDIR -M $POPMAP -O $OUTPUTDIR -t 12

#### COMPUTE POPULATION STATISTICS USING POPULATIONS ####
populations -P $OUTPUTDIR -M $POPMAP -t 12 --write_single_snp --genepop --vcf -r 0.50 --min_maf 0.02 #run populations at r=0.5 (output SNPs present in 50% of samples)
sed -i 's/# Pop ID/Pop_ID/g' populations.sumstats_summary.tsv #edit an output file so it's easier to use downstream
rename 's/populations./populations_r0.5./' populations.* #change the names of output files so they're unique to the r value

populations -P $OUTPUTDIR -M $POPMAP -t 12 --write_single_snp --genepop --vcf -r 0.80 --min_maf 0.02 #run populations at r=0.8 (output SNPs present in 80% of samples)
sed -i 's/# Pop ID/Pop_ID/g' populations.sumstats_summary.tsv #edit an output file so it's easier to use downstream
rename 's/populations./populations_r0.8./' populations.* #change the names of output files so they're unique to the r value

#make subdirectory for n-specific files and move all n-specific files into it
mkdir 'm'$m'M'$M'n'$n 
mv cstack.stats 'm'$m'M'$M'n'$n
mv Cstack.out 'm'$m'M'$M'n'$n
mv 'batch_'$n'.catalog.'*'.tsv' 'm'$m'M'$M'n'$n
mv *'.matches.tsv' 'm'$m'M'$M'n'$n
mv Sstack.out 'm'$m'M'$M'n'$n
mv *'.matches.bam' 'm'$m'M'$M'n'$n
mv tsv2bam.log 'm'$m'M'$M'n'$n
mv gstacks.* 'm'$m'M'$M'n'$n
mv populations* 'm'$m'M'$M'n'$n

#compress files (may want to just delete .tsv files)
gzip 'm'$m'M'$M'n'$n'/'*'.tsv' &
fi
}

##Run through range of n when M and m equal the defaults specified above, else hold n=default
if [ $M -eq $HOLDM ] && [ $m -eq $HOLDm ]
then

for n in `seq $n1 $n2`; do
nPipe
done

else
n=$HOLDn
nPipe
fi

#compress files (may want to just delete .tsv files)
gzip *.tsv &
}
########################################################################################################

###Launch StacksPipe function############
if [ $# -eq 13 ] #if we're running through a range of parameter values to explore optimal settings
then
for m in `seq $m1 $m2`; do #vary m between upper and lower bounds
M=$HOLDM #hold M
n=$HOLDn #hold n
StacksPipe #run through the main function
done

for M in `seq $M1 $M2`; do #vary M between upper and lower bounds
m=$HOLDm #hold m
n=$HOLDn #hold n
StacksPipe #run through the main function
done

else #otherwise we're running specific parameter values
m=$5 #m
M=$6 #M
HOLDn=$7 #n
StacksPipe #run through the main function
fi

