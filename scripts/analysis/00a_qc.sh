#!/bin/bash

# created by: Elizabeth Brooks

# script to perform fastqc quality control of paired end reads
# usage: bash 00a_qc.sh analysisType
# usage ex: bash 00a_qc.sh raw
# usage ex: bash 00a_qc.sh 01_merged
# usage ex: bash 00a_qc.sh 02_trimmed
# usage ex: bash 00a_qc.sh 03_filtered

# load software module for remote servers
module load bio/2.0

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# retrieve raw paired reads absolute path for alignment
	readPath=$(grep "pairedReads:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/pairedReads://g")
else 
	# set reads path
	readPath=$outputsPath"/"$analysisType
fi

#Make a new directory for analysis
outputsPath=$outputsPath"/00a_qc"
mkdir $outputsPath

#Make a new directory for analysis
outputsPath=$outputsPath"/qc_"$analysisType
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Move to the new directory
cd $outputsPath

# status message
echo "Beginning analysis..."

# perform QC
fastqc $readPath"/"*\.f*q* -o $outputsPath

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# run multiqc
	multiqc $outputsPath -o $outputsPath -n "qc_raw"
elif [[ $analysisType == "trimmed"* ]]; then
	# run multiqc on all
	multiqc $outputsPath -o $outputsPath -n "qc_"$analysisType"_all"
	# run multiqc on paired data
	#multiqc $outputsPath"/"*_p* -o $outputsPath -n "qc_"$analysisType"_paired"
else
	# run multiqc
	multiqc $outputsPath -o $outputsPath -n "qc_"$analysisType
fi

#Print status message
echo "Analysis complete!"
