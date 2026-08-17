#!/bin/bash

# created by: Elizabeth Brooks

# script to perform trimmomatic trimming of paired end reads
# usage: bash 01_merge.sh 

# DNA library: CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG(Nx40)CGCTGTCCGT -> 81bp total
# raw data: all reads are 81 bp

# retrieve paired reads absolute path for alignment
inputsPath=$(grep "pairedReads:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/pairedReads://g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve software absolute path
softwarePath=$(grep "software_FLASH:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_FLASH://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# create outputs directory, if needed
mkdir $outputsPath

# make a new directory for analysis
outputsPath=$outputsPath"/01_merged"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the software directory
cd $softwarePath

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*_R1_001\.fastq\.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R1_001\.fastq\.gz//')
	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R1_001\.fastq\.gz//')
	# status message
	echo "Processing $sampleTag"
	# perform merging of paired reads
	#./flash -t 4 -o $sampleTag -d $outputsPath -r 81 $f1 $f2
	./flash -t 4 -o $sampleTag -d $outputsPath -M 81 $f1 $f2
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
