#!/bin/bash
#$ -r n
#$ -N RNA_quantify_b_jobOutput
#$ -q largemem

# created by: Elizabeth Brooks

# script to count the number of sequences shared across runs
# usage: qsub 09a_quantify.sh inputRun runName
# usage ex: bash 09a_quantify.sh r8_S8_L001 r8_S8_L001

# retrieve input run name
inputRun=$1

# retrieve input run name
runName=$2

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")
#analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# setup inputs run data
inputData=$outputsPath"/05_combined/"$runName"_combined.RC.fa"

# retrieve input sequences
#inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted.fa" ## quantification of all sequencess
inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted_above2.fa"

# name of a new directory for analysis
tablesOut=$outputsPath"/09b_quantified_above2"

# make a new directory for analysis
mkdir $tablesOut

# setup tmp inputs data
inputRunData=$tablesOut"/"$inputRun"_"$runName"_combined.RC.tmp.fa"
cat $inputData > $inputRunData

# move to outputs directory
cd $tablesOut

# name formatted sequences file
fmtSeqs=$tablesOut"/"$inputRun"_"$runName"_formatted.tmp.fa"

# re-format input sequences for processing
cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' > $fmtSeqs

# add final new line
echo "" >> $fmtSeqs

# name output file
#countsOut=$tablesOut"/"$inputRun"_in_"$runName"_counts_table.csv"
countsPlotOut=$tablesOut"/"$inputRun"_in_"$runName"_counts_plot_table.csv"

# retrieve header
inputHeader="run_name,sequence_ID,read_counts,sequence"

# add the run names to the header 
#header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data outputs files
#echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis of $inputSeqs over $inputRunData ..."

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
	# count the number of seq occurances in each round
	numReads=$(cat $inputRunData | grep -wc $seq)
	# add the number of seqs for the round
	countDataOut=$(echo $countDataOut","$numReads)
	# add the counts data to the outputs file
	echo $countData","$numReads","$runName >> $countsPlotOut
	# add the counts data to the outputs file
	#echo $countDataOut >> $countsOut
done < $fmtSeqs

# clean up
rm $fmtSeqs
rm $inputRunData

# add run tags to sequence IDs
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv; do runTag=$(basename $i | cut -d"_" -f1 | sed "s/r//g"); tail -n+2 $i | awk -v runIn=$runTag 'BEGIN{FS=OFS=","}{$2 = runIn"_"$2; print}' > $i.fmt; done
# after processing the last round of data, combine all plotting data files
#head -1 /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r8_S8_L001_counts_plot_table.csv > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv.fmt; do tail -n+2 $i | grep -v "doped" >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv; done
# clean up
#rm /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv.fmt

# status message
echo "Analysis complete!"
