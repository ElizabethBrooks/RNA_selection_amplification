#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 02_combine.sh analysisType
# usage: bash 02_combine.sh trimmed_s4q20
# usage: bash 02_combine.sh trimmed_merged

# retrieve input analysis type
analysisType=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set the inputs directory
inputsPath=$outputsPath"/"$analysisType

# make a new directory for analysis
outputsCombined=$outputsPath"/combined_"$analysisTag
mkdir $outputsCombined
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsCombined directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsCombined

# status message
echo "Analyzing combined data..."

# unzip any gz read files
gunzip -v $inputsPath"/"*\.gz

# loop through all samples
for f1 in $inputsPath"/"*_stiched_reads\.fq; do
	# trim file extension
	sampleFile=$(echo $f1 | sed 's/\.fq//')
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.fq//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged, failed merged, and unpaired trimmed reads
	cat $f1 $sampleFile"_failed.fq_1.fastq" $sampleFile"_failed.fq_2.fastq" >> $outputsCombined"/"$sampleTag"_combined.fq"
	#cat $f1 $inputsPath"/trimmed_"$analysisTag"/"$sampleTag*\.fq >> $outputsCombined"/"$sampleTag"_combined.fq"
done

# double check that there are no duplicate reads
#for i in $outputsCombined"/"*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
