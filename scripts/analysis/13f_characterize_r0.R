#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

# R script to characterize the round 0 sequences

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(rcartocolor)
library(stringr)
library(dplyr)
library(ggrepel)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set the number of rounds
num_rounds <- 1

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_t0_run1"
#out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE)

# read in complement data
seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13e_conservation_r0_run1/overhang_data_wobble.csv"
#seqsFile <- args[2]
complement_data <- read.csv(seqsFile)

# replace NAs with zeros
complement_data[is.na(complement_data)] <- 0

# numbers of high quality reads
quality <- nrow(complement_data)

# keep sequences with at least a 4bp consecutive match
complement_data_dissimilar <- complement_data[complement_data$identity < 25, c("sequence")]

# output dissimilar sequence data
write.csv(complement_data_dissimilar, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_t0.csv", sep = ""), row.names = FALSE, quote = FALSE)

# replace 0's and 3's with <4
complement_data[complement_data$tag == 0, "tag"] <- "<5"
complement_data[complement_data$tag == 3, "tag"] <- "<5"
complement_data[complement_data$tag == 4, "tag"] <- "<5"
complement_data[complement_data$tag == 2, "tag"] <- "<5"

# get lists of identities and types
#tag_list <- unique(complement_data$tag)
tag_list <- c("8", "7", "4_3", "3_4", "6", "3_3", "5", "<5")
#tag_list <- c("8", "7", "4_3", "3_4", "6", "3_3", "5", "4", "3", "2", "0")

# add mapping table
identity_mappings <- data.frame(
  tag = tag_list,
  colors = safe_colors[1:8]
  #colors = safe_colors[1:11]
)

# add placeholder columns
complement_data$colors <- NA

# add run num column
complement_data$counts_run_name <- "r0_S0_L001"

# loop over each possible tag label
for (label_num in 1:nrow(identity_mappings)) {
  # set tag colors for plotting
  ifelse(
    complement_data$tag == identity_mappings$tag[label_num], 
    complement_data[complement_data$tag == identity_mappings$tag[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# determine the frequency of the numbers of bases (tags)
max_bases <- max(complement_data$tag, na.rm = TRUE)
min_bases <- min(complement_data$tag, na.rm = TRUE)
base_tbl <- table(complement_data$tag)
#base_names <- c("<4", "4", "5", "3_3", "6", "3_4", "4_3", "7", "8")
base_names <- names(base_tbl)
base_data <- data.frame(
  num_bases = as.character(base_names),
  freq_bases = c(unname(base_tbl)),
  #base_color = rev(safe_colors[1:length(base_tbl)])
  base_color = rev(c("#000", safe_colors[1:5], safe_colors[11:12]))
  #base_color = rev(c("#000", safe_colors[1:5], safe_colors[9:12]))
)
total_freq <- sum(base_data$freq_bases)
base_data$perc_freq <- (base_data$freq_bases/total_freq)*100
base_data_cleaned <- base_data %>% mutate(across(where(is.numeric), round, 2))

# bar chart of region frequencies
base_counts_plot <- ggplot(base_data_cleaned, aes(x=num_bases, y=perc_freq)) +
  #geom_bar(stat="identity", fill = base_data_cleaned$base_color) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(fill = 'black')) +
  #scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10), labels = function(x) paste0(x, "%")) +
  #guides(y = guide_axis(cap = "upper")) +
  geom_text(aes(label=paste0(sprintf("%1.1f", perc_freq),"%")),
            position=position_stack(vjust=0.5), size = 3, color = "white") +
  ylab("Proportion") +
  xlab("Complementarity")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total_t0_perc.png", sep = "")
#exportFile <- paste(out_dir, "/overhang_percent_abundance_total_t0_perc_all.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# output plotting data
write.csv(base_data_cleaned, file = paste(out_dir, "/overhang_conservation_t0.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(base_data_cleaned, file = paste(out_dir, "/overhang_conservation_t0_all.csv", sep = ""), row.names = FALSE, quote = FALSE)

# vectors of bins (total, consecutive, gaped)
#identity_bins <- unique(complement_data$identity)
tag_bins <- c("0", "2", "3", "4", "5", "3_3", "6", "3_4", "4_3", "7", "8")
plot_bins <- c(paste(tag_bins, "T", sep = "-"), paste(tag_bins, "C", sep = "-"), paste(tag_bins, "G", sep = "-"))

# set data lengths
data_length <- length(plot_bins)
seq_data_length <- nrow(complement_data)

# set the run num
run_num <- "0"

# set the counts run name
counts_run_num <- paste("r", run_num, "_S", run_num, "_L001", sep = "")

# initialize data frame for identity bin counts
complement_counts <- data.frame(
  run_name = rep(NA, data_length),
  tag = rep(NA, data_length),
  type = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length)
)

# loop over tag bins
for (bin_index in 1:data_length) {
  # retrieve current tag
  cur_tag <- strsplit(plot_bins[bin_index], split = "-")[[1]][1]
  # retrieve current type
  cur_type <- strsplit(plot_bins[bin_index], split = "-")[[1]][2]
  # add run ID
  complement_counts$run_name[bin_index] <- run_num
  # add tag
  complement_counts$tag[bin_index] <- cur_tag
  # add type
  complement_counts$type[bin_index] <- cur_type
  # check the type of complement
  if (cur_type == "T") { # total (including wobble)
    # add overhang complement unique counts
    complement_counts$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag,])
  }else if (cur_type == "C") { # consecutive (including wobble)
    # add overhang complement unique counts
    complement_counts$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag & complement_data$gap == "no",])
  }else if (cur_type == "G") { # gaped (including wobble)
    # add overhang complement unique counts
    complement_counts$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag & complement_data$gap == "yes",])
  }
  # add overhang complement counts
  complement_counts$counts[bin_index] <- complement_counts$counts_unique[bin_index]
  # add fraction abundance
  complement_counts$frac_abundance_unique[bin_index] <- complement_counts$counts_unique[bin_index]/seq_data_length
  complement_counts$frac_abundance[bin_index] <- complement_counts$counts[bin_index]/quality
}

# add percent counts
complement_counts$perc_abundance_unique <- 100*complement_counts$frac_abundance_unique
complement_counts$perc_abundance <- 100*complement_counts$frac_abundance

# update columns for plotting
complement_data_out <- data.frame(
  run_name = complement_data$counts_run_name,
  sequence_ID = paste("0", seq(1:nrow(complement_data)), sep = "_"),
  sequence = complement_data$sequence,
  counts = complement_data$counts,
  counts_run_name = complement_data$counts_run_name,
  complement = complement_data$complement,
  identity = complement_data$identity,
  identity_subset = complement_data$identity_subset,
  tag = complement_data$tag,
  tag_subset = complement_data$tag_subset,
  gap = complement_data$gap,
  wobble = complement_data$wobble,
  location = complement_data$location,
  all_locations = complement_data$all_locations,
  all_identities = complement_data$all_identities,
  all_tags = complement_data$all_tags
)

# export data
write.csv(complement_data_out, file = paste(out_dir, "/", counts_run_num, "_overhang_data_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts, file = paste(out_dir, "/", counts_run_num, "_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)

# create bar plot of total overhang identity percent
#base_counts_plot <- ggplot(complement_counts, aes(y=perc_abundance_unique, x=tag)) + 
#  geom_bar(stat="identity") +
#  theme_classic(base_size = 16) +
  #scale_fill_manual(breaks = unique(complement_counts$tag), values = unique(complement_counts$colors), labels = unique(complement_counts$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  #guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
#  labs(fill = "Complementary\nBases") +
#  ylab("Proportion") +
#  xlab("Complementarity")
# save the plot
#exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_r0_chart.png", sep = "")
#png(exportFile, units="in", width=5, height=4, res=300)
#print(base_counts_plot)
#dev.off()
