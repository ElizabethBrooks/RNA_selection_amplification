#!/bin/bash

# script to filter reads and keep sequences with matching up- and down-stream sequences
# usage: bash 05_filter.sh analysisType

# retrieve input analysis type
analysisType=$1

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/"$analysisType

# make a new directory for analysis
filterOut=$outputsPath"/filtered_"$analysisType
mkdir $filterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $filterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $filterOut

# status message
echo "Beginning analysis..."

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $inputsPath"/"*\.fq; do
	# status message
	echo "Processing $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_combined\.fq/_filtered\.fq/')
	# filter to keep sequences with matching up- and down-stream sequences
	#cat $f1 | awk '/GGACAGCG.{40}CGCTGTCC/{if (a && a !~ /GGACAGCG.{40}CGCTGTCC/) print a; print} {a=$0}' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' > $filterOut"/"$newName
	#cat $f1 | grep -Ex -B1 '.*GGACAGCG.{40}CGCTGTCC.*' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' | grep -v "^--$" > $filterOut"/"$newName
	cat $f1 | grep -Ex -B1 -A2 '.*GGACAGCG.{40}CGCTGTCC.*' | grep -v "^--$" > $filterOut"/"$newName
done

# status message
echo "Analysis complete!"
