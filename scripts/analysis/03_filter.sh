#!/bin/bash

# created by: Elizabeth Brooks

# script to filter reads and keep sequences with matching up- and down-stream sequences
# usage: bash 03_filter.sh

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# DNA library: CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG(Nx40)CGCTGTCCGT -> 81bp total

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set inputs path
inputsPath=$outputsPath"/02_trimmed"

# make a new directory for analysis
outputsPath=$outputsPath"/03_filtered"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# status message
echo "Beginning analysis..."

# unzip any gz read files
gunzip -v $inputsPath"/"*\.gz

# loop through all samples
for f1 in $inputsPath"/"*\.fq; do
	# status message
	echo "Processing $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_combined\.fq/_filtered\.fq/')
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | grep -Ex -B1 -A2 '.*CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG.{40}CGCTGTCCGT.*' | grep -v "^--$" > $outputsPath"/"$newName
done

# status message
echo "Analysis complete!"
