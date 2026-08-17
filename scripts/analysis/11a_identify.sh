#!/bin/bash
#$ -r n
#$ -N RNA_identify_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to run R scripts that count the number of sequences in sequence families
# usage: qsub 11a_identify.sh roundNum outDir peaksFile seqsFile

# load the software
module load bio/0724

# retrieve the input round number
#roundNum="1"
roundNum=$1

# set outputs directory
#outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification_above2"
#outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification_all"
outDir=$2

# read in cluster family sequence data
#peaksFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"
peaksFile=$3

# read in sequence count data for the specified round
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv"
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/counts_plot_table.csv"
seqsFile=$4

# retrieve input split num
splitNum=$5

# status message
echo "Beginning analysis of round $roundNum ..."

# run the R script
Rscript 11b_family_identification.R	$roundNum $outDir $peaksFile $seqsFile $splitNum

# clean up
rm $seqsFile

# status message
echo "Analysis complete!"
