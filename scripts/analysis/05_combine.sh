#!/bin/bash

# created by: Elizabeth Brooks

# script to combine files of trimmed forward paired reads with trimmed forward unpaired reads with merged reads
# usage: bash 05_combine.sh

# retrieve analysis tag
analysisTag=$(echo $1 | sed "s/trimmed_//g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/04_cleaned"

# make a new directory for analysis
outputsPath=$outputsPath"/05_combined"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsCombined

# loop through all samples
for f1 in $inputsPath"/"*_trimmed\.fa; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_trimmed\.fa//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/"$sampleTag"_pForward.fa" $inputsPath"/"$sampleTag"_uForward.fa" | perl -lpe'tr/ATGCatgc/TACGtacg/; $_ = reverse $_ unless $.%2;' >> $outputsPath"/"$sampleTag"_combined.RC.fa"
done

# double check that there are no duplicate reads
#for i in $outputsCombined"/"*_combined\.fa; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
