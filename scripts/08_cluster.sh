#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 32

# script to cluster sequences using clustalo
# usage: qsub 08_cluster.sh 

# load the software module
module load bio

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

# make a new directory for analysis
clusterOut=$outputsPath"/clustered"
mkdir $clusterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $clusterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis..."

# filter to keep sequences with matching up- and down-stream sequences
clustalo --threads=32 -v -i $formatOut"/combined.flt.fasta" -o $clusterOut"/clustered.fasta"

# status message
echo "Analysis complete!"