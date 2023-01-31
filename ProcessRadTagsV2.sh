#!/bin/bash
# Script written by Aryn Wilder, Asako Chaille and Hayley Nuetzel

## specify which version of STACKS to run
cd ' /home/Jamie/programs/stacks-2.41'

## here the path directories are specified for three files
## input fastq files, output files to go in, population map file
## in order so they are recognized as you drag them into terminal, in order 
FASTQDIR=$1 # path to directory of input fastq files
OUTDIR=$2 # path to output directory
POPMAP=$3 # path to POPMAP

## specify path to a directory of files (-p), data is paired-end reads (-P), input file type is gzfastq (-i), that the data is clean (-c), discard seq. with low QC (-q), disable rad check since Texas A&M already stripped RE sites (-disable_rad_check), and specify path for output (-o)
process_radtags -p $FASTQDIR -P -i gzfastq -c -q --disable_rad_check -o $OUTDIR

## reorganize output files by changing directory for output data
## rename output folders by removing the projectID and other unnecessary identifiers so just studbook number ID
cd $OUTDIR 
for file in *.fq.gz; do
mv $file "${file/19157Ivy_BP/}"
done

##make directories for the four datafile types
mkdir R1
mkdir R2
mkdir R1_rem
mkdir R2_rem


##move the four datafile types into their respective directories
mv *.rem.1.fq.gz R1_rem
mv *.rem.2.fq.gz R2_rem
mv *.1.fq.gz R1
mv *.2.fq.gz R2

cd R1
for i in *.1.fq.gz
do mv "$i" "$(echo $i | awk -F"_" '{print $1}').fq.gz"
done

cd ../R2
for i in *.2.fq.gz
do mv "$i" "$(echo $i | awk -F"_" '{print $1}').2.fq.gz"
done



####example code if you want to just remove .1 in 024.1.fq.gz to make file name 024.fq.gz
####for f in *.1.fq.gz; do
####mv -- "$f" "${f%.1.fq.gz}.fq.gz"
####done
