#!/bin/bash

# script to calculate the percent identity of cluster sequences compared with the peaks
# usage: bash 09_identify.sh inputFile
# usage ex: bash 09_identify.sh 07a_clustered
# usage ex: bash 09_identify.sh 07b_clustered

# retrieve input file
inputFile=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/08_summarized/"$inputFile

# make a new directory for analysis
tablesOut=$outputsPath"/09_identified"
mkdir $tablesOut

# make a new directory for analysis
tablesOut=$tablesOut"/"$inputFile
mkdir $tablesOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $tablesOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd $tablesOut

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*_above9_cluster_peaks_table\.csv; do
	# trim to sample tag
	nameTag=$(basename $f1 | sed 's/_peaks_table\.csv//g')
	# setup peak sequences file name
	f2=$inputsPath"/"$nameTag"_sequences_table.csv"
	# status message
	echo "Processing $nameTag ..."
	# initialize counter
	peakCount=0
	# loop over each cluster name number
	while read peakData; do
		# increment peak counter
		peakCount=$(($peakCount+1))
		# check if header line
		if [ $peakCount -eq 1 ]; then
			# adjust and output the header
			echo "run_name,sequence_ID,read_counts,cluster_ID,sequence,percent_ID" > $tablesOut"/"$nameTag"_sequences_identity_table.csv"
		else
			# get the current cluster number
			clusterName=$(echo $peakData | cut -d"," -f4)
			# status message
			echo "Processing Cluster $clusterName ..."
			# retrieve sequences for the current cluster
			#cat $f2 | grep ",$clusterName," > $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			cat $f2 | awk -F',' -v clusterIn="$clusterName" '$4==clusterIn' > $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			# initialize counters
			numSeqs=0
			totalID=0
			totalSqrVar=0
			# loop over each sequence line and retrieve sequences
			while read seqData; do
				# retrieve the current cluster peak sequence
				peakSeq=$(echo $peakData | cut -d"," -f6)
				# retrieve the current sequence in the current cluster
				currSeq=$(echo $seqData | cut -d"," -f5)
				# initialize mismatch counter
				numMismatch=0
				# loop over each character base in the sequence
				#echo $sequence | awk '{for (i=0; ++i <= length($0);) printf "%s", substr($0, i, 1)}'
				for (( i=0; i<${#currSeq}; i++ )); do
					# compare each character with the peak
					if [ "${currSeq:$i:1}" != "${peakSeq:$i:1}" ]; then
						# increment the number of mismatches
						numMismatch=$(($numMismatch+1))
					fi
				done
				# calculate the number of matches
				seqLength=40
				# calculate the percent identity
				distID=$(($seqLength - $numMismatch))
				propID=$(echo "scale=4; $distID / $seqLength" | bc)
				percentID=$(echo "scale=4; $propID * 100" | bc)
				# add the percent identity to the running total for the cluster
				totalID=$(echo "scale=4; $totalID + $percentID" | bc)
				# increase the counter for the number of sequences in the cluster
				numSeqs=$(($numSeqs + 1))
				# add the percent identity to the current sequence information and output to a new table
				echo $seqData","$percentID >> $tablesOut"/"$nameTag"_sequences_identity_table.csv"
			done < $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			# clean up
			rm $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
		fi
	done < $f1
done

# filter input sequences to keep only those with >= 90% ID to peak
for i in $tablesOut"/"*"_sequences_identity_table.csv"; do fileOut=$(echo $i | sed "s/\.csv//g"); cat $i | awk -F',' '$6 >= 90' > $fileOut"_atLeast90.csv"; done

# combine tables from each run
#head -1 $tablesOut"/"*"_sequences_identity_table.csv" > $tablesOut"/sequences_identity_table.csv"
#for i in $tablesOut"/"*"_sequences_identity_table.csv"; do tail -n+2 $i >> $tablesOut"/sequences_identity_table.csv"; done
#head -1 $tablesOut"/"*"_sequences_identity_table_atLeast90.csv" > $tablesOut"/sequences_identity_table_atLeast90.csv"
#for i in $tablesOut"/"*"_sequences_identity_table_atLeast90.csv"; do tail -n+2 $i >> $tablesOut"/sequences_identity_table_atLeast90.csv"; done

# status message
echo "Analysis complete!"
