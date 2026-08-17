#!/bin/bash
#$ -r n
#$ -N RNA_summarize_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to summarize clustering information
# usage: qsub 08_summarize.sh sampleTag
# usage ex: bash 08_summarize.sh r8_S8_L001

# retrieve input sample tag
sampleTag=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")
#analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/07_clustered_1500_full"

# make a new directory for analysis
tablesOut=$outputsPath"/08_summarized_1500_full"
mkdir $tablesOut

# move to outputs directory
cd $tablesOut

# status message
echo "Processing $sampleTag ..."

# set input file
f1=$inputsPath"/"$sampleTag"_clustered.aux"

# retrieve cluster names
cat $f1 | cut -d":" -f1 | sed "s/ /_/g" > $tablesOut"/"$sampleTag"_cluster_numbers.tmp.txt"
# retrieve sequence headers
cat $f1 | cut -d":" -f2 | tr -s " " | cut -d" " -f8 > $tablesOut"/"$sampleTag"_seq_headers.tmp.txt"
# combine cluster numbers with sequence headers
paste -d "," $tablesOut"/"$sampleTag"_seq_headers.tmp.txt" $tablesOut"/"$sampleTag"_cluster_numbers.tmp.txt" > $tablesOut"/"$sampleTag"_cluster_sequences_table.tmp.csv"
# add header to the output sequences file
echo "run_name,sequence_ID,read_counts,cluster_ID,sequence" > $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv"

# loop over each sequence line and retrieve sequence data
while read line; do
	# retrieve sequence header
	seqHeader=$(echo $line | cut -d"," -f1-3)
	# status message
	echo "Processing $seqHeader ..."
	# retrieve sequence data
	seqData=$(cat $outputsPath"/06_formatted/"$sampleTag"_formatted_above2.fa" | grep -A1 $seqHeader | tail -1)
	# add sequence data to the output file
	echo $line","$seqData >> $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv"
done < $tablesOut"/"$sampleTag"_cluster_sequences_table.tmp.csv"

# filter to keep only the peak sequences
cat $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv" | awk -F',' '!seen[$4]++' > $tablesOut"/"$sampleTag"_cluster_peaks_table.tmp.csv"
# remove the header
cat $tablesOut"/"$sampleTag"_cluster_peaks_table.tmp.csv" | tail -n+2 > $tablesOut"/"$sampleTag"_cluster_peaks_table.noHeader.tmp.csv"
# add header to the output sequences file
echo "run_name,sequence_ID,read_counts,cluster_ID,sequence_counts,sequence" > $tablesOut"/"$sampleTag"_cluster_peaks_table.tmp.csv"

# loop over each peak line and add the number of sequences and reads in the cluster
while read line; do
	# retrieve cluster number
	clusterName=$(echo $line | cut -d"," -f4)
	# status message
	echo "Processing $clusterName ..."
	# count the number of sequences in the current cluster
	clusterCount=$(cat $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv" | grep -c $clusterName)
	#clusterCount=$(cat $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv" | awk -F',' -v clusterIn="$clusterName" '$4==clusterIn')
	# count the number of reads in the cluster
	clusterReads=$(cat $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv" | grep $clusterName | cut -d"," -f3 | awk '{s+=$1} END {print s}')
	# separate sequence from the other data
	seqHeader=$(echo $line | cut -d"," -f1,2)
	clusterNum=$(echo $line | cut -d"," -f4)
	seqData=$(echo $line | cut -d"," -f5)
	# add cluster count to outputs
	echo $seqHeader","$clusterReads","$clusterNum","$clusterCount","$seqData >> $tablesOut"/"$sampleTag"_cluster_peaks_table.tmp.csv"
done < $tablesOut"/"$sampleTag"_cluster_peaks_table.noHeader.tmp.csv"

# re-format cluster names
cat $tablesOut"/"$sampleTag"_cluster_sequences_table.clust.tmp.csv" | sed "s/Cluster_//g" > $tablesOut"/"$sampleTag"_cluster_sequences_table.csv"
cat $tablesOut"/"$sampleTag"_cluster_peaks_table.tmp.csv" | sed "s/Cluster_//g" > $tablesOut"/"$sampleTag"_cluster_peaks_table.csv"

# clean up
rm $tablesOut"/"$sampleTag"_cluster_numbers.tmp.txt"
rm $tablesOut"/"$sampleTag"_seq_headers.tmp.txt"
rm $tablesOut"/"$sampleTag"_cluster_sequences_table"*".tmp.csv"
rm $tablesOut"/"$sampleTag"_cluster_peaks_table"*".tmp.csv"

# combine tables from each run
#cat $tablesOut"/"*"_cluster_sequences_table.csv" | head -1 > $tablesOut"/cluster_sequences_table.csv"
#for i in $tablesOut"/"*"_cluster_sequences_table.csv"; do tail -n+2 $i >> $tablesOut"/cluster_sequences_table.csv"; done
#cat $tablesOut"/"*"_cluster_peaks_table.csv" | head -1 > $tablesOut"/cluster_peaks_table.csv"
#for i in $tablesOut"/"*"cluster_peaks_table.csv"; do tail -n+2 $i >> $tablesOut"/cluster_peaks_table.csv"; done

# status message
echo "Analysis complete!"
