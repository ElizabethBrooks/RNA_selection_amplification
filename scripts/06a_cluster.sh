#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_clustalo_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 08a_cluster.sh inputFile
# usage ex: qsub 08a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r8_S8_L001_combined.fmt.fa
## job 816015 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 08a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r7_S7_L001_combined.fmt.fa
## job 816016 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 08a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r6_S6_L001_combined.fmt.fa
## job 816017
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 08a_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 08a_cluster.sh "${fileList[$i]}"; done
## job 814037 -> Couldn't allocate MPI memory
## job 814038 -> Couldn't allocate MPI memory
## job 814039 -> /opt/sge/crc/spool/d32cepyc207/job_scripts/814039: line 41: 2143111 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814040 -> Distance calculation within sub-clusters done. CPU time: 2166.49u 45.15s 00:36:51.64 Elapsed: 00:15:11
## job 814041 -> /opt/sge/crc/spool/d32cepyc189/job_scripts/814041: line 41: 3944610 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814042 -> /opt/sge/crc/spool/d32cepyc214/job_scripts/814042: line 41: 2583840 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814043 -> Distance calculation within sub-clusters done. CPU time: 2447.68u 48.10s 00:41:35.77 Elapsed: 00:17:54
## job 814044 -> /opt/sge/crc/spool/d32cepyc217/job_scripts/814044: line 41: 81853 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814045 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r3_S12_L001_combined_fmt/clustered_doped21-r3_S12_L001_combined_fmt
## job 814046 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r2_S11_L001_combined_fmt/clustered_doped21-r2_S11_L001_combined_fmt
## job 814047 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r1_S10_L001_combined_fmt/clustered_doped21-r1_S10_L001_combined_fmt
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_s4q20/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 08a_cluster.sh "${fileList[$i]}"; done
## job 868912 -> FATAL: Memory allocation for distance matrix failed
## job 868914 -> FATAL: Memory allocation for distance matrix failed
## job 868915 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r6_S6_L001_cleaned_fmt/clustered_r6_S6_L001_cleaned_fmt
## job 868916 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r5_S5_L001_cleaned_fmt/clustered_r5_S5_L001_cleaned_fmt
## job 868917 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r4_S4_L001_cleaned_fmt/clustered_r4_S4_L001_cleaned_fmt
## job 868918 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r3_S3_L001_cleaned_fmt/clustered_r3_S3_L001_cleaned_fmt
## job 868919 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r2_S2_L001_cleaned_fmt/clustered_r2_S2_L001_cleaned_fmt
## job 868920 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_r1_S1_L001_cleaned_fmt/clustered_r1_S1_L001_cleaned_fmt
## job 868921 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r3_S12_L001_cleaned_fmt/clustered_doped21-r3_S12_L001_cleaned_fmt
## job 868922 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r2_S11_L001_cleaned_fmt/clustered_doped21-r2_S11_L001_cleaned_fmt
## job 868923 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r1_S10_L001_cleaned_fmt/clustered_doped21-r1_S10_L001_cleaned_fmt

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisType=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_"$analysisType
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

# filter to keep sequences with matching up- and down-stream sequences
clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/"$nameTag

# status message
echo "Analysis complete!"
