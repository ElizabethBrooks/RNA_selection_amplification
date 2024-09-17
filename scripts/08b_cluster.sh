#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_500_clustalo_jobOutput
#$ -pe smp 8

# script to cluster sequences using clustalo and --cluster-size=500
# usage: qsub 08b_cluster.sh inputFile
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 08b_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 08b_cluster.sh "${fileList[$i]}"; done
## job 814062 -> FATAL: Memory allocation for distance matrix failed
## job 814063 -> FATAL: Memory allocation for distance matrix failed
## job 814064 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_size_500_r6_S6_L001_combined_fmt/clustered_r6_S6_L001_combined_fmt
## job 814065 
## job 814066 -> /opt/sge/crc/spool/d32cepyc217/job_scripts/814066: line 41: 87054 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500
## job 814067 -> /opt/sge/crc/spool/d32cepyc193/job_scripts/814067: line 41: 2721579 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500
## job 814068 -> /opt/sge/crc/spool/d32cepyc229/job_scripts/814068: line 41: 2345811 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500
## job 814069 -> /opt/sge/crc/spool/d32cepyc195/job_scripts/814069: line 41: 674954 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500
## job 814070 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_size_500_doped21-r3_S12_L001_combined_fmt/clustered_doped21-r3_S12_L001_combined_fmt
## job 814071 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_size_500_doped21-r2_S11_L001_combined_fmt/clustered_doped21-r2_S11_L001_combined_fmt
## job 814072 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_size_500_doped21-r1_S10_L001_combined_fmt/clustered_doped21-r1_S10_L001_combined_fmt

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# make a new directory for analysis
clusterOut=$outputsPath"/clustered_size_500_"$nameTag
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

# filter to keep sequences with matching up- and down-stream sequences
clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500

# status message
echo "Analysis complete!"