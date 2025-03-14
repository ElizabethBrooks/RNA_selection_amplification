#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_maxnumseq_clustalo_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo and --maxnumseq=1300000
# usage: qsub 06c_cluster.sh inputFile
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06c_cluster.sh "${fileList[$i]}"; done
## jobs

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_maxnumseq_1300000_"$analysisTag
mkdir $clusterOut

# make a new directory for analysis
clusterOut=$clusterOut"/"$nameTag
mkdir $clusterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $clusterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $nameTag ..."

# cluster sequences
clustalo --threads=$NSLOTS -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa" --maxnumseq=1300000

# status message
echo "Analysis complete!"
