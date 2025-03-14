#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_merge_NGmerge_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to perform merging of paired end reads into single reads
# usage: qsub 01b_trim.sh
## job 870256

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve paired reads absolute path for alignment
inputsPath=$(grep "pairedReads:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/pairedReads://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve software absolute path
softwarePath=$(grep "software_NGmerge:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_NGmerge://g")

# make a new directory for analysis
outputsPath=$outputsPath"/trimmed_"$analysisTag
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# create output logs directory
mkdir $outputsPath"/logs"

# move to the software directory
cd $softwarePath

# status message
echo "Processing..."

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
	./NGmerge -v -n 8 -m 20 -p 0.1 -1 $f1 -2 $f2 -o $outputsPath"/"$sampleTag"_stiched_reads.fq" -l $outputsPath"/logs/"$sampleTag"_log_stitching_results.txt" -f $outputsPath"/"$sampleTag"_stiched_reads_failed.fq" -j $outputsPath"/logs/"$sampleTag"_log_formatted_alignments.txt"
	#./bbmerge.sh in1=/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/trimmed_q20/r1_S1_L001_pForward.fq in2=/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/trimmed_q20/r1_S1_L001_pReverse.fq out=/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tests/merged_bbmerge/r1_merged.fq outu=/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tests/merged_bbmerge/r1_unmerged.fq ihist=/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tests/merged_bbmerge/r1_ihist.txt
	# status message
	echo "$sampleTag processed!"
done

#Print status message
echo "Analysis complete!"
