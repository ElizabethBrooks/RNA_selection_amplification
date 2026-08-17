#!/bin/bash
#$ -r n
#$ -N RNA_conserved_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to run R scripts that count the number of sequences in sequence families
# usage: qsub 13_conserved.sh roundNum outDir seqsFile

# load the software
module load bio/0724

# retrieve the input round number
#roundNum="1"
roundNum=$1

# set outputs directory
#outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_all"
#outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_above2"
outDir=$2

# read in sequence count data for the specified round
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/counts_plot_table_noDoped.csv"
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table_noDoped.csv"
seqsFile=$3

# status message
echo "Beginning analysis of round $roundNum ..."

# run the R script
Rscript 13_round_conservation.R $roundNum $outDir $seqsFile

# status message
echo "Analysis complete!"
