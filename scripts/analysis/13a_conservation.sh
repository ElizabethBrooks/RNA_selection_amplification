#!/bin/bash

# created by: Elizabeth Brooks

# script to run job scripts that identify conserved regions
# usage: bash 13a_conservation.sh

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13a_overhang_conservation_all"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
seqsFile=$outputsPath"/09a_quantified_all/counts_plot_table_noDoped.csv"

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
