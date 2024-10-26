#!/bin/bash

# script to clean reads and keep only the variable 40bp region
# usage: bash 04_clean.sh

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/filtered"

# make a new directory for analysis
outputsPath=$outputsPath"/cleaned"
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

# loop through all samples
for f1 in $inputsPath"/"*\.fq; do
	# status message
	echo "Processing $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/\.fq/\.fa/')
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | sed "s/^.*CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG/START/g" | sed "s/CGCTGTCCGT.*$/END/g" | grep -Ex -B1 'START.{40}END' | grep -v "^--$" | sed "s/START//g" | sed "s/END//g" > $outputsPath"/"$newName
done

# status message
echo "Analysis complete!"
