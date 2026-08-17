#!/bin/bash
#$ -r n
#$ -N RNA_randomize_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to create files with randomized nucelotide sequences
# usage: qsub 14_randomize.sh numNum
# usage example: qsub 14_randomize.sh 1
## 1 to 4
## jobs 1735273 to 1735276

# retrieve input run num
runNum=$1

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/14_randomized_sequences_2M"

# create outputs directory
mkdir $outDir

# set of nucleotides 
nucSet="ACTG"

# length of seqeuences
seqLength=40

# length of file
numSeqs=2000000

# name out output file
outFile=$outDir"/random"$runNum"_sequences_combined.RC.fa"

# add header
echo "run_name,sequence_ID,read_counts,sequence,counts,counts_run_name" > $outFile

# loop over each line of the file
for seqNum in $(seq 1 $numSeqs); do
	# status message
	echo "Processing sequence number $seqNum"
	# reset the current sequence
	currSeq=""
	# loop over each nucleotide
	for nucNum in $(seq 1 $seqLength); do
		# generate a random number
		nucIndex=$((RANDOM % 4))
		#nucIndex=$(shuf -i 1-4 -n 1)
		# select a random nucleotide
		randomNuc="${nucSet:$nucIndex:1}"
		# add the random nucelotide to the sequence
		currSeq="$currSeq$randomNuc"
	done
	# add the randomized sequence to the file
	echo "random$runNum,$seqNum,1,$currSeq,1,$runNum" >> $outFile
done

# status message
echo "Analysis complete!"
