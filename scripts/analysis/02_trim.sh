#!/bin/bash

# created by: Elizabeth Brooks

# script to perform trimmomatic trimming of paired end reads
# usage: bash 02_trim.sh 

# DNA library: CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG(Nx40)CGCTGTCCGT -> 81bp total

# load the software module
module load bio/2.0

# set phred score for trimming
# https://www.drive5.com/usearch/manual/quality_score.html
score=33

# retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter_TruSeq:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/adapter_TruSeq://g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set inputs path
inputsPath=$outputsPath"/01_merged"

# make a new directory for analysis
outputsPath=$outputsPath"/02_trimmed"
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

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/"*\.extendedFrags.fastq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.extendedFrags\.fastq//')
	# status message
	echo "Processing $sampleTag"
	# perform adapter trimming on paired reads
	# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
	#trimmomatic PE -threads 4 -phred"$score" $inputsPath"/"$sampleTag".notCombined_1.fastq" $inputsPath"/"$sampleTag".notCombined_2.fastq" $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" MAXINFO:56:0.8 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	#trimmomatic SE -threads 4 -phred"$score" $f1 $sampleTag"_trimmed.fq.gz" MAXINFO:56:0.8 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	#trimmomatic PE -threads 4 -phred"$score" $inputsPath"/"$sampleTag".notCombined_1.fastq" $inputsPath"/"$sampleTag".notCombined_2.fastq" $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" SLIDINGWINDOW:4:20 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	#trimmomatic SE -threads 4 -phred"$score" $f1 $sampleTag"_trimmed.fq.gz" SLIDINGWINDOW:4:20 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	trimmomatic PE -threads 4 -phred"$score" $inputsPath"/"$sampleTag".notCombined_1.fastq" $inputsPath"/"$sampleTag".notCombined_2.fastq" $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" AVGQUAL:30 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	trimmomatic SE -threads 4 -phred"$score" $f1 $sampleTag"_trimmed.fq.gz" AVGQUAL:30 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
