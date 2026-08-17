#!/bin/bash
#$ -r n
#$ -N RNA_conserved_r0_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to run job scripts that identify conserved regions
# usage: qsub 13e_conservation_r0.sh runNum
# usage ex: qsub 13e_conservation_r0.sh 1

# load the software
module load bio/0724

# retrieve input run number
runNum=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13e_conservation_r0_run"$runNum

# create outputs directory
mkdir $outDir

# status message
echo "Processing random$runNum sequences..."

# read in sequence count data for the specified round
seqsFile="/scratch365/ebrooks5/RNA_evolution/outputs/14_randomized_sequences/random"$runNum"_sequences_combined.RC.fa"

# run the analysis
Rscript 13_sequence_conservation.R $outDir $seqsFile

# characterize and summarize the results
Rscript 13f_characterize_r0.R $outDir $outDir"/overhang_data_wobble.csv"

# status message
echo "Analysis complete!"
