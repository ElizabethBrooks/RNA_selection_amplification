#!/bin/bash
#$ -r n
#$ -N RNA_quantify_d_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to count the number of sequences shared across the top 10 sequences for the runs
# usage: qsub 09b_quantify_top10.sh inputRun
# usage ex: bash 09b_quantify_top10.sh r8_S8_L001

# retrieve input run name
inputRun=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")
#analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/06_formatted"

# retrieve input sequences
#inputSeqs=$inputsPath"/"$inputRun"_formatted.fa" ## quantification of all sequences
inputSeqs=$inputsPath"/"$inputRun"_formatted_above2.fa"

# process just the top 10 most abundant sequences
tablesOut=$outputsPath"/09d_quantified_top10_above2"

# make a new directory for analysis
mkdir $tablesOut

# directory for tmp inputs
#seqsInput=$tablesOut"/tmp"

# make a new directory for tmp inputs
#mkdir $seqsInput
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $seqsInput directory already exsists... skipping creation of tmp inputs."
#else
	# setup tmp inputs data
#	cp $outputsPath"/05_combined/"*_combined\.RC\.fa > $seqsInput
#fi

# move to outputs directory
cd $tablesOut

# name formatted sequences file
fmtSeqs=$tablesOut"/"$inputRun"_formatted.tmp.fa"

# re-format input sequences for processing
# process just the top 10 most abundant sequences
#cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' | sort -k3 -n | head -10 > $fmtSeqs
cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' | head -10 > $fmtSeqs

# add final new line
#echo "" >> $fmtSeqs

# name output file
#countsOut=$tablesOut"/"$inputRun"_counts_table.csv"
countsPlotOut=$tablesOut"/"$inputRun"_counts_plot_table.csv"

# retrieve header
inputHeader="run_name,sequence_ID,read_counts,sequence"

# add the run names to the header 
#header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data outputs files
#echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis of $inputRun ..."

# loop over round sequences
while read data; do
	# clean up the run name
	seqData=$(echo $data | sed "s/>r//g")
	# remove the sequence
	dataNoSeq=$(echo $seqData | cut -d"," -f1-3)
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f4)
	# reverse compliment the sequence
	#seqRev=$(echo $seq | tr ACGTacgt TGCAtgca | rev)
	# update the count data
	#countData=$dataNoSeq","$seqRev
	countData=$dataNoSeq","$seq
	countDataOut=$countData
	# status message
	echo "Processing $seq ..."
	# loop over each round sequences file
	for f2 in $outputsPath"/05_combined/"*_combined\.RC\.fa; do
		# retrieve run name
		runName=$(basename $f2 | sed "s/_S.*_L001_combined\.RC\.fa//g" | sed "s/r//g" | sed "s/21-/_/g")
		# count the number of seq occurances in each round
		numReads=$(cat $f2 | grep -wc $seq)
		# add the number of seqs for the round
		countDataOut=$(echo $countDataOut","$numReads)
		# add the counts data to the outputs file
		echo $countData","$numReads","$runName >> $countsPlotOut
	done
	# add the counts data to the outputs file
	#echo $countDataOut >> $countsOut
done < $fmtSeqs

# clean up
rm $fmtSeqs
#rm $tablesOut"/"*_combined\.RC\.fa

# combine run data
#head -1 /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/r8_S8_L001_counts_plot_table.csv > /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/counts_plot_table.csv
#for i in /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/*_counts_plot_table.csv; do echo $i; tail -n+2 $i >> /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/counts_plot_table.csv; done
#cat /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/counts_plot_table.csv | grep -v "doped" > /scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified/counts_plot_table_noDoped.csv

# status message
echo "Analysis complete!"
