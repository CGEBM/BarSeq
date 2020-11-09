#!/bin/bash
############################################################################################################################################################
# barcodeDetectionCodeDispatcher - A "dispatcher" component to submit jobs using its counterpart "generalBBDUK4BarSeq.pl" wrapper script over the BBDuk tool for NGS Bar-seq analysis.
# Copyright (C) 2015-present Antonio Claudio Bello Ribeiro and contributors
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#	Contact details: antonio.ribeiro@abdn.ac.uk or at Centre for Genome-Enabled Biology and Medicine (CGEBM) 23 St Machar Drive - Aberdeen AB24 3RY - United Kingdom
############################################################################################################################################################

baseDir=$1 # The ABSOLUTE PATH (just like when obtained via pwd command; i.e. without trailing slashes) of the base directory that contains the samples' subfolders (usually these latter are the ones that contain the *fastq.gz files)
mismatchRate=$2 # The aimed mismatch rate being used in the experiment (for usage in the bbduk third-party script later) (user must choose between the values 0 to 3)
list_file=$3 # The ABSOLUTE PATH to the FASTA-like file with the list of barcode sequences being searched in the experiment (for usage in the bbduk third-party script later)
time=$4 # The INT time period (in hours) that should be allocated for each job to run in the cluster (user must choose an integer value of at least 1)
readType=$5 # The type of read file to be analysed by the script (user must choose between the value 1 for analysing R1 reads or 2 for analysing R2 reads)

# Checking for the correct number of args
if [ $# -ne 5 ]
then
	echo "Error in $0 - Invalid argument count"
	echo "Syntax: $0 <absolute path to the base directory that contains the samples subfolders (as when obtained via pwd command; i.e. without trailing slashes)> <mismatch rate allowed when looking for the barcode sequence (options: 0 to 3)> <absolute path to the list of barcodes FASTA file> <INT value representing the time (in hours) that should be allocated for each job to run> <R1 or R2 read type to be analysed (options: 1 for R1 reads or 2 for R2 reads)>"
	exit
fi
#if [ $2 -lt 0 ] || [ $2 -gt 3 ]
if (($2 < 0 || $2 > 3))
then
	echo "Sorry, I'm not prepared to deal with this mismatch rate input... Please, inform an INTEGER number between 0 and 3."
	exit
fi
#if [ $4 -gt 1 ]
if (($4 <= 0))
then
	echo "Sorry, I'm not prepared to deal with this time input... Please, inform an INTEGER number of at least 1."
	exit
fi
#if [ $5 -ne 1 ] || [ $5 -ne 2 ]
if (($5 < 1 || $5 > 2))
then
	echo "Sorry, I'm not prepared to deal with this read type input... Please, opt for 1 for analysing R1 reads or 2 for analysing R2 reads."
	exit
fi

### Strings and formatting variables
underscore="_"
slash="/"
reportString="report.txt"

printf " \n"
echo "Heating up..."
printf " \n"

######################################################################################################################################################################
### Obs.: Please, bear in mind that this script and its associate "generalBBDUK4BarSeq.pl" counterpart have been developed/tested using openjdk-8.0.152 (to enable the Java environment for BBDUK). So, when running the script, the following SLURM module loading instruction is also issued. Depending on your system/environment, you might need to adapt this command to the appropriate version available (or event to omit/comment it out).
######################################################################################################################################################################
module load openjdk-8.0.152
######################################################################################################################################################################

# Listing all of the directories in this dir and extracting all of their names to an array
dirs=`ls -l | grep "^d" | awk '{ print $11 }'` # Please, bear in mind that this instruction might need modification depending on the Unix system: typically, field $9 is another possibility...

# For each directory
for i in $dirs
do
	echo "Hopping to directory $i"
	printf " \n"
	cd $i
	# Filenames' convention examples that could be present in each directory
	#98_PXXXA_1-ds.541e63d043484d54b1be4ecbf42e9d90/98-PXXXA-1_S127_L001-4_R1_001.fastq.gz
	#98_PXXXA_1-ds.541e63d043484d54b1be4ecbf42e9d90/98-PXXXA-1_S127_L001-4_R2_001.fastq.gz
	#99_PXXXA_1-ds.b15273753c5b4ff6bf49f9afd14dc037/99-PXXXA-1_S126_L001-4_R1_001.fastq.gz
	#99_PXXXA_1-ds.b15273753c5b4ff6bf49f9afd14dc037/99-PXXXA-1_S126_L001-4_R2_001.fastq.gz
	#9_PXXXA_1-ds.1174eadeb5bc4d65ba30aab8546f85ab/9-PXXXA-1_S55_L001-4_R1_001.fastq.gz
	#9_PXXXA_1-ds.1174eadeb5bc4d65ba30aab8546f85ab/9-PXXXA-1_S55_L001-4_R2_001.fastq.gz

	# Getting some substrings from the filenames
	echo "Getting some substrings from the filenames which will be needed later..."
	printf " \n"
	echo "I'll work with R${readType} reads..."
	printf " \n"
	file2Work=`ls -l *_L001_R${readType}_001.fastq.gz | awk '{ print $11 }'` # Please, bear in mind that this instruction might need modification depending on the Unix system: typically, field $9 is another possibility...
	echo "Echoing the names of R${readType} files..."
	printf " \n"
	# Iterate over these
	for m in $file2Work
	do
		echo $m
	done
	testFilename2Work=`echo $m`
	# Extracting filename without extension
	testFilename2WorkNoExt=`echo $testFilename2Work | basename $testFilename2Work .fastq.gz`
	printf " \n"
	echo "The name of the file that I'll use as a test is: $testFilename2Work"
	printf " \n"
	echo "The basename of this file is: $testFilename2WorkNoExt"
	printf " \n"
	numUnds=`echo $testFilename2WorkNoExt | awk -F_ '{print NF-1}'`
	echo "I am dealing with $numUnds underscores in the basename $testFilename2WorkNoExt"
	numUndsMinus1=`expr $numUnds - 1`; echo -e "\t numUndsMinus1= $numUndsMinus1"
	numUndsMinus2=`expr $numUnds - 2`; echo -e "\t numUndsMinus2= $numUndsMinus2"
	numUndsMinus3=`expr $numUnds - 3`; echo -e "\t numUndsMinus3= $numUndsMinus3"
	numUndsMinus4=`expr $numUnds - 4`; echo -e "\t numUndsMinus4= $numUndsMinus4"
	printf " \n"
	sampleName=`awk -F '[_.]' -v xx=$numUnds '{print $(NF-xx)}' <<< "$testFilename2WorkNoExt"`; echo -e "\t sampleName= $sampleName"
	sampleNumber=`awk -F '[_.]' -v xx=$numUndsMinus1 '{print $(NF-xx)}' <<< "$testFilename2WorkNoExt"`; echo -e "\t sampleNumber= $sampleNumber"
	lane=`awk -F '[_.]' -v xx=$numUndsMinus2 '{print $(NF-xx)}' <<< "$testFilename2WorkNoExt"`; echo -e "\t lane= $lane"
	read=`awk -F '[_.]' -v xx=$numUndsMinus3 '{print $(NF-xx)}' <<< "$testFilename2WorkNoExt"`; echo -e "\t read= $read"
	flowCellIndex=`awk -F '[_.]' -v xx=$numUndsMinus4 '{print $(NF-xx)}' <<< "$testFilename2WorkNoExt"`; echo -e "\t flowCellIndex= $flowCellIndex"
	printf " \n"
	
	# Submitting the barcode detection code
	echo "Running the barcode detection code over the $file2Work reads present in the directory..."
	printf " \n"
	#######################################################
	# Running the command without using a scheduler system:
	#######################################################
	#perl ~/CGEBM/BarSeq/generalBBDUK4BarSeq.pl $baseDir$slash$i$slash$files$file2Work $mismatchRate $list_file $sampleName$underscore$sampleNumber$underscore$read$underscore$reportString
	#####################################################################################################
	# Original submission command for Sun Grid Engine scheduler by the time of the prototype elaboration:
	#####################################################################################################
	#qsub -cwd -V -l h_vmem=10G -l h_rt=$time:0:0 -b y perl ~/CGEBM/BarSeq/generalBBDUK4BarSeq.pl $baseDir$slash$i$slash$files$file2Work $mismatchRate $list_file $sampleName$underscore$sampleNumber$underscore$read$underscore$reportString
	###############################################################################
	# Adaptation of the command to the SLURM scheduler (development/test scenario):
	###############################################################################
	###############################################################################################################
	### Just echoing the final command in the development/test scenario while making different tests with the code:
	###############################################################################################################
	#echo "sbatch --mem=10G --time=$time:00:00 --wrap=\"perl ~/CGEBM-open/NGS-BarSeq/generalBBDUK4BarSeq.pl $baseDir$slash$i$slash$file2Work $mismatchRate $list_file $sampleName$underscore$sampleNumber$underscore$read$underscore$reportString\""
	sbatch --mem=10G --time=$time:00:00 --wrap="perl ~/CGEBM/BarSeq/generalBBDUK4BarSeq.pl $baseDir$slash$i$slash$file2Work $mismatchRate $list_file $sampleName$underscore$sampleNumber$underscore$read$underscore$reportString"
	printf " \n"
	echo "I'm done with this directory... The barcode detection code should have been submitted over the $file2Work file of the folder in question... Hopping back to the upper level folder..."
	printf " \n"
	cd $baseDir
done
echo "JOB DONE!!!"
printf " \n"
