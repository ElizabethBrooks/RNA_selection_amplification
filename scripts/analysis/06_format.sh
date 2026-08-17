#!/bin/bash

# created by: Elizabeth Brooks

# script to subset sequences and format headers
# usage: bash 06_format.sh

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/05_combined"
#inputsPath=$outputsPath"/14_randomized_sequences"

# make a new directory for analysis
outputsPath=$outputsPath"/06_formatted"
#outputsPath=$inputsPath"/06_formatted"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd $outputsPath

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*_combined\.RC\.fa; do
	# status message
	echo "Processing file: $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_combined\.RC\.fa/_formatted/')
	# print read counts
	# for fasta files
	cat $f1 | awk 'NR%2==0' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName"_counts.tmp.txt"
	# for fastq files
	#cat $f1 | awk 'NR%4==2' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName"_counts.tmp.txt"
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $outputsPath"/"$newName"_counts.tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./>$runTag/g" > $outputsPath"/"$newName"_run.tmp.txt"
	# make a file with a sequence of numbers to use as unique sequence IDs
	seq -f %1.0f 1 $outLength > $outputsPath"/"$newName"_ID.tmp.txt"
	# combine output counts and sequence file with run tags and convert to csv
	paste -d" " $outputsPath"/"$newName"_run.tmp.txt" $outputsPath"/"$newName"_ID.tmp.txt" $outputsPath"/"$newName"_counts.tmp.txt" | tr -s ' ' | sed "s/ /,/g" > $outputsPath"/"$newName".tmp.txt"
	# filter to remove sequences with less than 3 reads
	awk -F ',' '($3 > 2)' $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_above2.tmp.txt"
	# cut out the header data
	cut -d"," -f1-3 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_header.tmp.txt"
	cut -d"," -f1-3 $outputsPath"/"$newName"_above2.tmp.txt" > $outputsPath"/"$newName"_above2_header.tmp.txt"
	# cut out the seqeunce data
	cut -d"," -f4 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_seq.tmp.txt"
	cut -d"," -f4 $outputsPath"/"$newName"_above2.tmp.txt" > $outputsPath"/"$newName"_above2_seq.tmp.txt"
	# interleave the seqeunce data with the headers
	paste -d'\n' $outputsPath"/"$newName"_header.tmp.txt" $outputsPath"/"$newName"_seq.tmp.txt" > $outputsPath"/"$newName".fa"
	paste -d'\n' $outputsPath"/"$newName"_above2_header.tmp.txt" $outputsPath"/"$newName"_above2_seq.tmp.txt" > $outputsPath"/"$newName"_above2.fa"
	# clean up
	rm $outputsPath"/"$newName*".tmp."*
done

# create combined data file for plotting
#cat /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/*_L001_formatted.fa | sed "s/$/,/g" | tr -d '\n' | sed "s/>/\n>/g" | sed "s/,$//g"  > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa  
#tail -n+2 /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa.tmp                                       
#echo "" >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa.tmp
#echo "run_name,sequence_ID,counts,sequence" > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa    
#cat /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa.tmp | sed "s/>r//g" >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/combined_formatted_no_doped.fa 

# status message
echo "Analysis complete!"
