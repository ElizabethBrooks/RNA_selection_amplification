#!/bin/bash

# created by: Elizabeth Brooks

# script to run job scripts that identify conserved regions
# usage: bash 13b_conservation.sh
## test_23May2025
## jobs 1733990 to 1733997

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13b_overhang_conservation_above2"
#outDir=$outputsPath"/test_22May2025/13b_overhang_conservation_above2"
#outDir=$outputsPath"/test_23May2025/13b_overhang_conservation_above2"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
seqsFile=$outputsPath"/09b_quantified_above2/counts_plot_table_noDoped.csv"

# loop over each input run num
for runNum in {1..8}; do 
	# retrieve the input round number
	roundNum=$runNum
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# submit job script
	qsub 13_conserved.sh $roundNum $outDir $seqsFile	
done

# status message
echo "Analysis complete!"
