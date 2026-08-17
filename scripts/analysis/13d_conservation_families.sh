#!/bin/bash

# created by: Elizabeth Brooks

# script to run job scripts that identify conserved regions
# usage: bash 13d_conservation_families.sh

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13d_overhang_conservation_families"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
seqsFile="/scratch365/ebrooks5/RNA_evolution/outputs/tables_and_figures/ST2_family_table/r8_family_count_data.csv"
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/ST2_family_table/r8_family_count_data.csv"

# run the analysis
Rscript 13_sequence_conservation.R $outDir $seqsFile

# status message
echo "Analysis complete!"
