#!/bin/bash

# created by: Elizabeth Brooks

# script to run job scripts that identify conserved regions
# usage: bash 13c_conservation_top10.sh

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13c_overhang_conservation_top10_above2"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
seqsFile=$outputsPath"/09d_quantified_top10_above2/counts_plot_table_noDoped.csv"

# run the analysis
Rscript 13_sequence_conservation.R $outDir $seqsFile

# status message
echo "Analysis complete!"
