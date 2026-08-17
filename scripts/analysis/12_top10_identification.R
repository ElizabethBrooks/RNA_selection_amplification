#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

# R script to identifiy the top 10 most abundant sequneces
# and their similarity to the peak sequences

# turn of scientific notation
options(scipen=10000)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/12_top10_identification"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# read in cluster family sequence data
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv")

# read in sequence count data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09c_quantified_top10_all/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
  
# remove duplicate sequence data
seqs_input <- seqs_input[!duplicated(seqs_input$sequence),]

# set the sequence length
seqLength <- 40

# setup data frame length
data_length <- nrow(seqs_input) * nrow(r8_peaks)

# data frame of sequence data
seqs_data <- data.frame(
  run_name = rep(NA, data_length),
  sequence_ID = rep(NA, data_length),
  sequence = rep(NA, data_length),
  peak_cluster_ID = rep(NA, data_length),
  peak_identity = rep(NA, data_length)
)

# loop over each sequence and compare with the peak to note if it has >= 90% similarity
for (seq_num in 0:(data_length-1)) {
  # loop over rach peak sequence
  for (peak_num in 1:nrow(r8_peaks)) {
    # compare the current sequence with the current peak
    numMatch <- mapply(
      function(x, y) {
        len <- length(x)
        sum(x[1:len] == y[1:len])
      }, 
      strsplit(seqs_input$sequence[seq_num+1], ''), 
      strsplit(r8_peaks$sequence[peak_num], '')
    )
    # determine identity percent
    identity_perc <- 100*numMatch/seqLength
    # set the index
    index <- ((seq_num)*nrow(r8_peaks)) + peak_num
    # record data
    seqs_data$run_name[index] <- seqs_input$run_name[seq_num+1]
    seqs_data$sequence_ID[index] <- seqs_input$sequence_ID[seq_num+1]
    seqs_data$sequence[index] <- seqs_input$sequence[seq_num+1]
    seqs_data$peak_cluster_ID[index] <- r8_peaks$cluster_ID[peak_num]
    seqs_data$peak_identity[index] <- identity_perc
  }
}

# export data
write.csv(seqs_data, file = paste(out_dir, "/family_identities_top10.csv", sep = ""), row.names = FALSE, quote = FALSE)

# check how many sequences have >= 90% identity to each peak
for (cluster_num in 0:(nrow(r8_peaks)-1)) {
  print(cluster_num)
  print(nrow(seqs_data[seqs_data[seqs_data$peak_cluster_ID == cluster_num,]$peak_identity >= 90,]))
}

# keep sequences that have >= 90% identity to any peak
seqs_out <- seqs_data[seqs_data$peak_identity >= 90,]

# export data
write.csv(seqs_out, file = paste(out_dir, "/family_identities_atLeast90_top10.csv", sep = ""), row.names = FALSE, quote = FALSE)
