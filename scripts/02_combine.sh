#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 02_combine.sh

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set inputs path
inputsPath=$outputsPath"/trimmed_"$analysisTag

# make a new directory for analysis
outputsPath=$outputsPath"/combined_"$analysisTag
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# status message
echo "Analyzing combined data..."

# unzip any gz read files
gunzip -v $inputsPath"/"*\.gz

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/"*_pForward\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_pForward\.fq//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/"$sampleTag"_pReverse.fq" $inputsPath"/"$sampleTag"_uForward.fq" $inputsPath"/"$sampleTag"_uReverse.fq" > $outputsPath"/"$sampleTag"_combined.fq"
done

# loop through all samples
#for f1 in $inputsPath"/"*_stiched_reads\.fq; do
	# trim file extension
	#sampleFile=$(echo $f1 | sed 's/\.fq//')
	# trim to sample tag
	#sampleTag=$(basename $f1 | sed 's/\.fq//')
	# status message
	#echo "Processing $sampleTag ..."
	# combine un-filtered merged, failed merged, and unpaired trimmed reads
	#cat $f1 $sampleFile"_failed.fq_1.fastq" $sampleFile"_failed.fq_2.fastq" > $outputsPath"/"$sampleTag"_combined.fq"
#done

# double check that there are no duplicate reads
#for i in $outputsPath"/"*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
