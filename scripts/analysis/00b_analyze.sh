#!/bin/bash

# created by: Elizabeth Brooks

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 00b_analyze.sh analysisType
# usage ex: bash 00b_analyze.sh 01_merged
# usage ex: bash 00b_analyze.sh 02_trimmed
# usage ex: bash 00b_analyze.sh 03_filtered
# usage ex: bash 00b_analyze.sh 04_cleaned
# usage ex: bash 00b_analyze.sh 05_combined
# usage ex: bash 00b_analyze.sh 06_formatted

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
outputsPath=$outputsPath"/00b_analyzed"
mkdir $outputsPath

#Make a new directory for analysis
outputsPath=$outputsPath"/analyzed_"$analysisType
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

# unzip any gz read files
gunzip -v $readPath"/"*

# name output files
namesOut=$outputsPath"/names_count.txt"
uniqueOut=$outputsPath"/unique_sequences_count.txt"
lengthsOut=$outputsPath"/read_lengths.txt"
countsOut=$outputsPath"/read_counts.txt"

# add headers
echo "run,length" > $lengthsOut
echo "run,count,sequence" > $countsOut
#echo "run,total,quality,unique,diversity"

# check analysis type
if [[ $analysisType == "cleaned" || $analysisType == "combined" || $analysisType == "formatted" ]]; then
	# count read names
	echo "Calculating read name counts..."
	for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%2==1' | wc -l; done > $namesOut
	# count unique read sequences
	echo "Calculating unique read sequence counts..."
	#for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%2==0' | sort -u | wc -l; done > $uniqueOut
	for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%2==0' | awk '!seen[$0]++' | wc -l; done > $uniqueOut
else
	# count read names
	echo "Calculating read name counts..."
	for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%4==1' | wc -l; done > $namesOut
	# count unique read sequences
	echo "Calculating unique read sequence counts..."
	#for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%4==2' | sort -u | wc -l; done > $uniqueOut
	for i in $readPath"/"*; do echo $i; cat $i | awk 'NR%4==2' | awk '!seen[$0]++' | wc -l; done > $uniqueOut
fi

# loop over each file
for f1 in $readPath"/"*; do
	# status message
	echo "Processing file: $f1"
	# check analysis type
	if [[ $analysisType == "cleaned" || $analysisType == "combined" || $analysisType == "formatted" ]]; then
		# print read lengths
		#cat $f1 | awk 'NR%2==0 {print length($0)}' | sort -n | uniq -c > $lengthsOut".tmp.txt"
		cat $f1 | awk 'NR%2==0 {print length($0)}' | awk '!seen[$0]++' > $lengthsOut".tmp.txt"
		# print read counts
		cat $f1 | awk 'NR%2==0' | sort | uniq -c | sort -nrk1 > $countsOut".tmp.txt"
		#cat $f1 | awk 'NR%2==0' | awk '!seen[$0]++' | sort -rk1 > $countsOut".tmp.txt"
	else
		# print read lengths
		#cat $f1 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c > $lengthsOut".tmp.txt"
		cat $f1 | awk 'NR%4==2 {print length($0)}' | awk '!seen[$0]++' > $lengthsOut".tmp.txt"
		# print read counts
		cat $f1 | awk 'NR%4==2' | sort | uniq -c | sort -nrk1 > $countsOut".tmp.txt"
		#cat $f1 | awk 'NR%4==2' | awk '!seen[$0]++' | sort -rk1 > $countsOut".tmp.txt"
	fi
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $lengthsOut".tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./$runTag/g" > $lengthsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $lengthsOut".run.tmp.txt" $lengthsOut".tmp.txt" | tr -s ' ' | sed "s/ /,/g" >> $lengthsOut
	# get length of output file
	outLength=$(wc -l $countsOut".tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./$runTag/g" > $countsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $countsOut".run.tmp.txt" $countsOut".tmp.txt" | tr -s ' ' | sed "s/ /,/g" >> $countsOut
	# clean up
	rm $lengthsOut".tmp.txt"
	rm $lengthsOut".run.tmp.txt"
	rm $countsOut".tmp.txt"
	rm $countsOut".run.tmp.txt"
done

# status message
echo "Analysis complete!"
