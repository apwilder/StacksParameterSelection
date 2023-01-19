# StacksParameterSelection
This program selects Stacks de novo parameter (M, m and n) values that maximize the number of SNPs called in 80% of samples

Steps for running Stacks parameter optimization pipeline

1.	After processing RADtags, run StacksDeNovo.sh (in the directory /media/jamie/4TB/TAM_run01/StacksPipeline/) in exploratory mode:
a.	This program builds de novo stacks from RADseq fastq files, and can be used to find optimal parameter settings for m, M and n when run in exploratory mode
b.	There are two ways to run this program: 
i.	explore a range of each parameter value of m, M and n, holding other parameters constant at the default values
ii.	specify m, M and n and the program will run only those values
c.	The script is set up so that there is a folder with R1 reads and a folder with R2 reads.
i.	The script assumes files are named something like SampleID.OptionalOtherinformation.fq (a dot-delimited string where text before the first dot is the sample ID, and the text after the last dot is .fq. This can be changed in the script on lines 55, 60 and 62.
d.	The PopulationMap.txt file is a tab-delimited list of sample ID’s and populations
i.	Sample ID’s should be the file name before the .fq
e.	To run exploratory mode use a command like:

nohup ./StacksDeNovo.sh /media/jamie/4TB/TAM_run01/R1/ \
/media/jamie/4TB/TAM_run01/R2/ \
/media/jamie/4TB/TAM_run01/stacks_denovo_output/ \
/media/jamie/4TB/TAM_run01/PopulationMap.txt \	
1	3 1 3 1 10 3 2 0 >& outfile.nohup &

f.	The nohup … >& outfile.nohup & parts of the command will run the script in no hangup mode, so you can close your connection and the program will keep running in the background. Any output that would normally print to the screen will print to a file called outfile.nohup
g.	The other parameters are:
a.	R1DIR: the R1 directory
b.	R2DIR: the R2 direcory
c.	OUTPUTDIR: the Output directory
d.	POPMAP: the population map
e.	m1: #Range of m (lower bound) [or m if running specific parameter values]
f.	m2: #Range of m (upper bound) [or M if running specific parameter values]
g.	M1: #Range of M (lower bound) [or n if running specific parameter values]
h.	M2: #Range of M (upper bound) [or nothing if running specific parameter values]
i.	n1: #Range of n (lower bound) [or nothing if running specific parameter values]
j.	n2: #Range of n (upper bound) [or nothing if running specific parameter values]
k.	HOLDm: #The default value to hold m when varying other parameters [or nothing if running specific parameter values]
l.	HOLDM: #The default value to hold M when varying other parameters [or nothing if running specific parameter values]
m.	HOLDn: #The default value to hold n when varying other parameters [or nothing if running specific parameter values]


2.	After StacksDeNovo.sh is finished, run GetSharedStacks_v2.sh with the command:

nohup ./GetSharedStacks_v2.sh /media/jamie/4TB/TAM_run01/PopulationMap.txt \
/media/jamie/4TB/TAM_run01/stacks_denovo_output/ \
Speciesname >& outfile2.nohup &

This will create summary tables from the results of the stacks runs that can be plotted in R

3.	Finally, to plot a summary, run PlotStacksSummary.R with the command:

Rscript PlotStacksSummary.R SpeciesName m1 m2 M1 M2 n1 n2
(where m1 and m2 are the lower and upper bounds of the parameters you explored for m, and so on.)
This will plot the number of shared SNPs across samples at various parameter values.

4.	Once you’ve found your optimal parameters, run them again in StacksDeNovo.sh in set parameter mode with the command (assuming m=3, M=2 and n=0 are optimal):

nohup ./StacksDeNovo.sh /media/jamie/4TB/TAM_run01/R1/ \
/media/jamie/4TB/TAM_run01/R2/ \
/media/jamie/4TB/TAM_run01/stacks_denovo_output/ \
/media/jamie/4TB/TAM_run01/PopulationMap.txt \
3 2 0 >& outfile.nohup &

