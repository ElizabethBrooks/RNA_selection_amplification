#!/bin/bash
#$ -r n
#$ -N RNA_summarize_jobOutput
#$ -q largemem

# script to run R scripts that count the number of sequences in sequence families
# usage: qsub 11a_identification.sh roundNum outDir peaksFile seqsFile
# usage ex: qsub 11a_identification.sh 1 /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv

# load the software
module load bio/2.0

# retrieve the input round number
#roundNum="1"
roundNum=$1

# set outputs directory
#outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification"
outDir=$2

# read in cluster family sequence data
#peaksFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"
peaksFile=$3

# read in sequence count data for the specified round
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv"
seqsFile=$4

# status message
echo "Beginning analysis of round $roundNum ..."

# run the R script
Rscript 11a_family_identification.R	$roundNum $outDir $peaksFile $seqsFile

# status message
echo "Analysis complete!"
