#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_format_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to subset sequences and format headers
# usage: qsub 05_format.sh inputFile
# usage: qsub 05_format.sh cleaned_merged

# retrieve input analysis type
analysisType=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/"$analysisType

# make a new directory for analysis
outputsPath=$outputsPath"/formatted_"$analysisTag
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
for f1 in $inputsPath"/"*\.fa; do
	# status message
	echo "Processing file: $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_cleaned\.fq/_formatted\.fq/')
	# print read counts
	# for fasta files
	cat $f1 | awk 'NR%2==0' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName".counts.tmp.txt"
	# for fastq files
	#cat $f1 | awk 'NR%4==2' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName".counts.tmp.txt"
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $outputsPath"/"$newName".counts.tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./>$runTag/g" > $outputsPath"/"$newName".run.tmp.txt"
	# combine output counts and sequence file with run tags and convert to csv
	paste -d" " $outputsPath"/"$newName".run.tmp.txt" $outputsPath"/"$newName".counts.tmp.txt" | tr -s ' ' | sed "s/ /,/g" > $outputsPath"/"$newName".tmp.txt"
	# cut out the header data
	cut -d"," -f1,2 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName".header.tmp.txt"
	# cut out the seqeunce data
	cut -d"," -f3 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName".seq.tmp.txt"
	# interleave the seqeunce data with the headers
	paste -d'\n' $outputsPath"/"$newName".header.tmp.txt" $outputsPath"/"$newName".seq.tmp.txt" > $outputsPath"/"$newName
	# clean up
	rm $outputsPath"/"$newName*".tmp.txt"
done

# status message
echo "Analysis complete!"
