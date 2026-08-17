#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

# R script to create top 10 sequence abundance plots
# usage: 03_F2D_F3_T2_top10_sequences.R

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(dplyr)
library(ComplexHeatmap)
library(stringr)
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

# supress cowplot package messages
suppressMessages( require(cowplot) )

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F2D_F3_T2_top10_sequences_above2_log10"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F2D_F3_T2_top10_sequences_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F2D_F3_T2_top10_sequences_all"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000", "#D55E00", "#0072B2", "#F0E442", "#E69F00", "#009E73", "#999999")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
#diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
#quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09d_quantified_top10_above2/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character", "sequence_ID"="character"))
#seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09c_quantified_top10_all/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character", "sequence_ID"="character"))

# replace counts with 0 for sequences with less than 2 reads
seqs_counts[seqs_counts$counts < 3,"counts"] <- 0

# reverse complement the sequences
#seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# add the fraction abundances
#seqs_counts$frac_abundance <- NA
#for (run_num in 1:8) {
#  seqs_counts[seqs_counts$counts_run_name == run_num, "frac_abundance"] <- seqs_counts[seqs_counts$counts_run_name == run_num, "counts"]/quality[run_num]
#}

# add log values
seqs_counts$log_counts <- log10(seqs_counts$counts)
#seqs_counts$log_counts <- log(seqs_counts$counts)
#seqs_counts$log_frac_abundance <- log(seqs_counts$frac_abundance)

# add percent abundance
#seqs_counts$perc_abundance <- 100*seqs_counts$frac_abundance

# set infinite and NA values equal to zero
is.na(seqs_counts)<-sapply(seqs_counts, is.infinite)

# re-order the data for plotting
#seqs_counts <- seqs_counts[order(seqs_counts$log_counts, decreasing=TRUE),]

# setup midpoint values for plotting
seqs_counts_noNA <- seqs_counts
seqs_counts_noNA[is.na(seqs_counts_noNA)] <- 0
mid_log_counts <- max(seqs_counts_noNA$log_counts)/2
#mid_log_frac_abundance <- log(min(seqs_counts_noNA[seqs_counts_noNA$frac_abundance != 0,"frac_abundance"]))/2
#mid_perc_abundance <- max(seqs_counts$perc_abundance)/2

# create ID factor list
ID_list <- rep("1",8)
#ID_list <- rep(1,8)
for (i in 2:10) {
  ID_list <- c(ID_list, rep(i,8))
}
ID_data <- as.factor(ID_list)
#ID_data <- factor(ID_data, levels = levels(ID_data)[c(1, 3:10, 2)])
ID_data <- factor(ID_data, levels = levels(ID_data)[rev(c(1, 3:10, 2))])

# initialize plot data
counts_heatmap_list <- vector('list', 8)

# loop over each round and create heatmaps
for (run_num in 1:8) {
  # subset seq data
  seqs_counts_subset <- seqs_counts[seqs_counts$run_name == run_num,]
  # remove run tags
  seqs_counts_subset[c("run", "ID")] <- do.call(rbind, strsplit(seqs_counts_subset$sequence_ID, "_"))
  #seqs_counts_subset[c("run", "Sequence")] <- do.call(rbind, strsplit(seqs_counts_subset$sequence, "_"))
  # add sequence top 10 ID
  seqs_counts_subset$ID <- ID_data
  # set round plot title
  run_title <- paste("Round", run_num)
  # set run indexes for plotting
  run_start <- 0.5 + (run_num - 1)
  run_end <- 1.5 + (run_num - 1)
  # To-do: consider using factors
  # create heatmap of log counts
  #counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(ID), log_counts), fill= log_counts)) + 
  counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, ID, fill= log_counts)) + 
  #counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(Sequence), log_counts), fill= log_counts)) + 
    theme_classic(base_size = 16) +
    geom_tile(colour = "black") +
    ggtitle(run_title) +
    theme(#axis.title.x = element_text(size = 14), 
          #axis.title.y = element_text(size = 14), 
          #axis.text.x = element_text(size = 14), 
          #axis.text.y = element_text(size = 14), 
          #text = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5)) +
    #ylab("Sequence ID") +
    ylab("Sequence") +
    xlab("Round") +
    scale_fill_gradient2(name = "Log Counts",
                         low = "#F0E442",
                         mid = safe_colors[7],
                         high = safe_colors[5],
                         midpoint = mid_log_counts,
                         na.value = "white") +
    coord_fixed() +
    annotate("rect", xmin = run_start, xmax = run_end, ymin = 0.5, ymax = 10.5, 
             colour = safe_colors[2], fill = "transparent", linewidth = 1.25)
    #scale_y_discrete(labels=seqs_counts_subset$ID)
  # save the plot
  exportFile <- paste(out_dir, "/r", run_num, "_top10_sequence_log_counts.png", sep = "")
  #png(exportFile, units="in", width=5, height=5, res=300)
  png(exportFile, units="in", width=10, height=5, res=300)
  print(counts_heatmap_subset)
  dev.off()
  # add the log counts heatmap to the plot list
  counts_heatmap_list[[run_num]] <- counts_heatmap_subset
}

# combine heatmaps with the log counts for each of the top 10 sequences per round
grid_log_counts_plot <- plot_grid(plotlist=counts_heatmap_list,  ncol = 4, align = 'v')
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_log_counts.png", sep = "")
png(exportFile, units="in", width=20, height=10, res=300)
print(grid_log_counts_plot + theme_classic(base_size = 14))
dev.off()

# export plotting data
write.csv(seqs_counts, file = paste(out_dir, "/top10_sequence_counts.csv", sep = ""), row.names = FALSE, quote = FALSE)

# create table of top 10 sequences per run
sequence_data <- data.frame(
  round = c(rep("1", 10), 
            rep("2", 10), 
            rep("3", 10), 
            rep("4", 10), 
            rep("5", 10), 
            rep("6", 10),
            rep("7", 10), 
            rep("8", 10)),
  ranking = c(rep(seq(from=1, to=10, by=1), 8)),
  shared = c(rep("1", 10), 
             rep("2", 10), 
             rep("3", 10), 
             rep("4", 10), 
             rep("5", 10), 
             rep("6", 10),
             rep("7", 10), 
             rep("8", 10)),
  sequence = c(seqs_counts[seqs_counts$run_name == 1 & seqs_counts$counts_run_name == 1, "sequence"],
               seqs_counts[seqs_counts$run_name == 2 & seqs_counts$counts_run_name == 2, "sequence"],
               seqs_counts[seqs_counts$run_name == 3 & seqs_counts$counts_run_name == 3, "sequence"],
               seqs_counts[seqs_counts$run_name == 4 & seqs_counts$counts_run_name == 4, "sequence"],
               seqs_counts[seqs_counts$run_name == 5 & seqs_counts$counts_run_name == 5, "sequence"],
               seqs_counts[seqs_counts$run_name == 6 & seqs_counts$counts_run_name == 6, "sequence"],
               seqs_counts[seqs_counts$run_name == 7 & seqs_counts$counts_run_name == 7, "sequence"],
               seqs_counts[seqs_counts$run_name == 8 & seqs_counts$counts_run_name == 8, "sequence"])
)

# determine overlap of top 10 sequences
# loop over each sequence
for (seq_num in 1:nrow(sequence_data)) {
  # initialize outputs variable
  match_list <- ""
  # loop over each round
  for (round_num in 1:8) {
    # loop over each top sequence
    for (top_num in 1:10) {
      # check if current sequence is in the current round top sequences
      if (sequence_data$sequence[seq_num] == seqs_counts[seqs_counts$run_name == round_num & seqs_counts$counts_run_name == round_num, "sequence"][top_num]) {
        match_list <- paste(match_list, round_num, sep="_")
      }
    }
  }
  # add matching rounds to outputs
  sequence_data$shared[seq_num] <- substr(match_list, 2, nchar(match_list))
}

# export sequence data
#write.csv(sequence_data, file = paste(out_dir, "/top10_sequences.csv", sep = ""), row.names = FALSE, quote = FALSE)

# re-organize sequence ranking data 
# sequence, Round1_rank, Round2_rank, Round3_rank, etc.
unique_seqs <- unique(sequence_data$sequence)
data_length <- length(unique_seqs)
sequence_ranks <- data.frame(
  sequence = unique_seqs,
  Round1 = rep(NA, data_length),
  Round2 = rep(NA, data_length),
  Round3 = rep(NA, data_length),
  Round4 = rep(NA, data_length),
  Round5 = rep(NA, data_length),
  Round6 = rep(NA, data_length),
  Round7 = rep(NA, data_length),
  Round8 = rep(NA, data_length)
)

# loop over each run
for (seq_num in 1:data_length) {
  # retrieve the current sequence
  curr_seq <- sequence_ranks$sequence[seq_num]
  # add the round rankings
  sequence_ranks$Round1[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 1,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 1, "ranking"])
  sequence_ranks$Round2[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 2,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 2, "ranking"])
  sequence_ranks$Round3[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 3,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 3, "ranking"])
  sequence_ranks$Round4[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 4,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 4, "ranking"])
  sequence_ranks$Round5[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 5,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 5, "ranking"])
  sequence_ranks$Round6[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 6,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 6, "ranking"])
  sequence_ranks$Round7[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 7,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 7, "ranking"])
  sequence_ranks$Round8[seq_num] <- ifelse(nrow(sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 8,]) == 0, NA, sequence_data[sequence_data$sequence ==  curr_seq & sequence_data$round == 8, "ranking"])
}

# add column of ranked seq names
#seq_names <- paste("Rank", seq(1, 19), sep ="_")
#seq_names <- seq(1, 19)
seq_names <- c("Seq1_1", "Seq1_2", "Seq1_3", "Seq1_4", "Seq1_5", "Seq1_6", "Seq1_7", "Seq1_8", "Seq1_9", "Seq1_10", "Seq2_9", "Seq2_10", "Seq3_6", "Seq3_9", "Seq4_10", "Seq5_8", "Seq6_10", "Seq8_7", "Seq8_8")
sequence_ranks$rank_ID <- seq_names

# export sequence data
write.csv(sequence_ranks, file = paste(out_dir, "/top10_sequences_rankings.csv", sep = ""), row.names = FALSE, quote = FALSE)

# retrieve data for ranked sequences
seqs_counts_ranked <- seqs_counts[seqs_counts$sequence %in% sequence_ranks$sequence,]
seqs_counts_ranked <- seqs_counts_ranked[ , !(names(seqs_counts_ranked) %in% c("run_name","sequence_ID","read_counts"))]
seqs_counts_ranked <- seqs_counts_ranked %>% distinct()
seqs_counts_ranked$rank_ID <- rep("NA", nrow(seqs_counts_ranked))

# additional color pallete
ranked_color_pallete <- c(safe_colors[1:10], palette.colors(palette = "Okabe-Ito"))

# loop over each ranked sequence
for (rank_num in 1:nrow(sequence_ranks)) {
  # add column of ranked seq names
  seqs_counts_ranked[seqs_counts_ranked$sequence == sequence_ranks$sequence[rank_num], "rank_ID"] <- sequence_ranks$rank_ID[rank_num]
  # add plotting colorss
  seqs_counts_ranked[seqs_counts_ranked$sequence == sequence_ranks$sequence[rank_num], "rank_color"] <- ranked_color_pallete[rank_num]
}

# make a factor variable
seqs_counts_ranked$rank_ID <- as.factor(seqs_counts_ranked$rank_ID)
#levels(seqs_counts_ranked$rank_ID)
seqs_counts_ranked$rank_ID <- factor(seqs_counts_ranked$rank_ID, levels = levels(seqs_counts_ranked$rank_ID)[rev(c(1, 12:19, 2:11))])
#levels(seqs_counts_ranked$rank_ID)

# calculate percent abundance
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
seqs_counts_ranked$read_abun <- rep("NA", nrow(seqs_counts_ranked))
# loop over each run
for (run_num in 1:8) {
  # add read abundances
  seqs_counts_ranked[seqs_counts_ranked$counts_run_name == run_num, "read_abun"] <- 100*seqs_counts_ranked[seqs_counts_ranked$counts_run_name == run_num, "counts"]/quality[run_num]
}
seqs_counts_ranked$read_abun <- as.numeric(seqs_counts_ranked$read_abun)

# create heatmap of read abundances
#counts_heatmap_ranked <- ggplot(data = seqs_counts_ranked, aes(counts_run_name, rank_ID, fill= read_abun)) + 
  #theme_classic(base_size = 16) +
  #geom_tile(colour = "black") +
  #geom_line() +
  #ggtitle(run_title) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ylab("Sequence ID") +
  #ylab("Sequence") +
  #xlab("Round") +
  #scale_fill_gradient2(name = "Log Counts",
  #                     low = safe_colors[3],
  #                     mid = safe_colors[4],
  #                     high = safe_colors[5],
                       #midpoint = mid_log_counts,
  #                     na.value = "white") +
  #coord_fixed()

# create line plot of non-fam ranked seqs read abundances
seqs_counts_ranked_subset <- seqs_counts_ranked[seqs_counts_ranked$rank_ID %in% c(seq_names[5], seq_names[6], seq_names[8], seq_names[10], seq_names[12], seq_names[13], seq_names[14], seq_names[16], seq_names[17], seq_names[19]),]
#sequence_ranks_subset <- sequence_ranks[sequence_ranks$rank_ID %in% c(5, 6, 8, 10, 12, 13, 14, 16, 17, 19),]
# loop over each ranked sequence
#for (rank_num in 1:nrow(sequence_ranks_subset)) {
#  seqs_counts_ranked_subset[seqs_counts_ranked_subset$sequence == sequence_ranks_subset$sequence[rank_num], "rank_color"] <- safe_colors[rank_num]
#}
counts_plot_ranked_subset <- ggplot(data=seqs_counts_ranked_subset, aes(x=as.character(counts_run_name), y=read_abun, group=rank_ID, color=rank_color))+
  geom_line(size = 1.25) +
  geom_point(size = 2.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Seq ID", labels = seqs_counts_ranked_subset$rank_ID, breaks = seqs_counts_ranked_subset$rank_color, guide = "legend") +
  scale_y_continuous(limits=c(0, 4), breaks=seq(0, 4, 1)) +#, labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Percent Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/ranked_sequences_non_family_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(counts_plot_ranked_subset)
dev.off()

# add ranks
#seqs_counts_ranked$rank <- str_split_fixed(seqs_counts_ranked$rank_ID, "_", 2)[,2]
#seqs_counts_ranked[seqs_counts_ranked$read_abun == 0,]$read_abun <- NA

#counts_plot_ranked <- ggplot(data = seqs_counts_ranked, aes(x = counts_run_name, y = rank, group = rank_ID)) +
#  geom_line(aes(color = rank_ID, alpha = 1), size = 2) +
#  geom_point(aes(color = rank_ID, alpha = 1), size = 4) #+
  #scale_y_reverse(breaks = 1:nrow(seqs_counts_ranked))
# save the plot
#exportFile <- paste(out_dir, "/ranked_sequences_non_family_percent_abundances.png", sep = "")
#png(exportFile, units="in", width=5, height=4, res=300)
#print(counts_plot_ranked)
#dev.off()

# export plotting data
write.csv(seqs_counts_ranked, file = paste(out_dir, "/ranked_sequences_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
