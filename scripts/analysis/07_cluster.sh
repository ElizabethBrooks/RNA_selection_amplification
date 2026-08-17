#!/bin/bash

# created by: Elizabeth Brooks

# script to perform CD-Hit clustering of RNA sequence data
# usage: bash 07_cluster.sh

# perform clustering
#cd-hit-est -i /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/r8_S8_L001_formatted_above2.fa -o /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered/r8_S8_L001_formatted_above2  -c 0.90 -n 8

# summarize the results
#plot_len1.pl /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr 1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999 > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr_data.txt

# add header to the output results
echo "run_name,sequence_ID,read_counts,cluster_ID,sequence_counts,sequence" > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr_peak_data.txt

# initialize counters and flags
clustFlag=0
seqCount=0

# loop over each cluster
while read LINE; do
	# check if the current line is the cluster header
	if [[ $LINE == ">Cluster"* ]]; then
		# check if the current cluster is not the first
		if [[ $clustFlag -eq 1 ]]; then
			# format the peak sequnece data
			outData=$peakData","$currCluster","$seqCount","$peakSeq
			# output the previous cluster data
			echo $outData >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr_peak_data.txt
		fi
		# record current cluster
		currCluster=$(echo $LINE | sed "s/>//g" | cut -d" " -f2)
		# trigger cluster flag for outputs
		clustFlag=1
		# re-set the sequence counter
		seqCount=0
	else
		# check if the current sequence is a representative peak sequence
		if [[ $LINE == *"... *" ]]; then
			# retrieve the peak sequence data
			peakData=$(echo $LINE | cut -d">" -f2 | cut -d"." -f1)
			# retrieve the peak sequence
			peakSeq=$(cat /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/r8_S8_L001_formatted_above2.fa | grep -A1 $peakData | tail -n+2)
			# increment the sequence counter
			seqCount=$(($seqCount+1))
		else
			# increment the sequence counter
			seqCount=$(($seqCount+1))
		fi
	fi
done < /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr
