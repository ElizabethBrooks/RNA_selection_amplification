#!/bin/bash

# created by: Elizabeth Brooks

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 11_identification.sh

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/11b_family_identification_above2"

# create outputs directory
mkdir $outDir

# set the cluster family sequence data
peaksFile=$outputsPath"/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"

# set the sequence count data for the specified round
seqsFile=$outputsPath"/09b_quantified_above2/counts_plot_table.csv"
countsFile=$outputsPath"/09b_quantified_above2/counts_plot_table_noDoped.csv"

# create the combined counts plotting data
head -1 $outputsPath"/09b_quantified_above2/r1_S1_L001_in_r1_S1_L001_counts_plot_table.csv" > $seqsFile
for i in $outputsPath"/09b_quantified_above2/"*_counts_plot_table.csv; do cat $i | tail -n+2 >> $seqsFile; done
cat $seqsFile | grep -v "doped" > $countsFile

# loop over each input run num
for runNum in {1..8}; do 
	# retrieve the input round number
	roundNum=$runNum
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# submit job script
	qsub 11a_identify.sh $roundNum $outDir $peaksFile $countsFile	
done

# status message
echo "Analysis complete!"
