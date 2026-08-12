#!/usr/bin/env Rscript

# R script to create overhang conservation plots
# usage: 07_family_overhang_conservation.R

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(rcartocolor)
library(stringr)
library(dplyr)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set the number of rounds
num_rounds <- 8

# numbers of high quality reads
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#above_two_reads <- c(18, 19, 26, 27, 1585, 10626, 7230, 6315)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2"
dir.create(out_dir, showWarnings = FALSE)

# set outputs sub-directory
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family1"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family3"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family4"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family5"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family6"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family7"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family8"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family9"
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_above2/family10"

# create outputs sub-directory
dir.create(out_dir, showWarnings = FALSE)

# read in identity data
identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_data_wobble.csv"
complement_data_rounds <- read.csv(identityFile)

# read in family data
familyFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2_r8_unique/family_sequences.csv"
family_data <- read.csv(familyFile)

# retrieve family data
#family_data <- family_data[family_data$family_ID == 1,]
#family_data <- family_data[family_data$family_ID == 2,]
#family_data <- family_data[family_data$family_ID == 3,]
#family_data <- family_data[family_data$family_ID == 4,]
#family_data <- family_data[family_data$family_ID == 5,]
#family_data <- family_data[family_data$family_ID == 6,]
#family_data <- family_data[family_data$family_ID == 7,]
#family_data <- family_data[family_data$family_ID == 8,]
#family_data <- family_data[family_data$family_ID == 9,]
family_data <- family_data[family_data$family_ID == 10,]

# subset to family data
complement_data_family <- merge(complement_data_rounds, family_data, by = "sequence")

# remove unnecessary columns
complement_data_family <- subset(complement_data_family, select = -c(run_name.y, sequence_ID.y, counts.y, counts_run_name.y))

# re-name columns
names(complement_data_family)[names(complement_data_family) %in% c("run_name.x", "sequence_ID.x", "counts.x", "counts_run_name.x")] <- c("run_name", "sequence_ID", "counts", "counts_run_name")

# replace NAs with zeros
complement_data_family[is.na(complement_data_family)] <- 0
#complement_data_rounds[is.na(complement_data_rounds)] <- 0

# copy data frame for later plotting
complement_data <- complement_data_family
#complement_data <- complement_data_rounds

#keep sequences with at least a 4bp consecutive match (100*2/8 = 25)
#complement_data_similar <- complement_data[complement_data$identity >= 25, c("sequence", "counts", "counts_run_name")]
complement_data_dissimilar <- complement_data[complement_data$identity <= 25, c("run_name", "sequence_ID", "sequence", "counts", "counts_run_name", "complement", "identity")]
complement_data_dissimilar_round8 <- complement_data_dissimilar[complement_data_dissimilar$counts_run_name == "r8_S8_L001", c("run_name", "sequence_ID", "sequence", "counts", "complement", "identity")]

# read in sequences that have at least 90% identity to any peak
#seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max_atLeast90.csv")
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max.csv")
seqs_identities <- seqs_identities[c("run_name", "sequence_ID", "sequence", "peak_cluster_ID", "peak_identity")]

# read in family data
fam_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/ST2_family_table/r8_family_count_data.csv")

# setup cluster list
cluster_list <- fam_data$cluster_ID

# loop over each cluster
seqs_identities$family_ID <- NA
for (cluster_num in min(cluster_list):max(cluster_list)) {
  # add family IDs
  seqs_identities[seqs_identities$peak_cluster_ID == cluster_num, "family_ID"] <- fam_data[fam_data$cluster_ID == cluster_num, "family_ID"]
}

# add families
complement_data_dissimilar_round8_identities <- right_join(seqs_identities, complement_data_dissimilar_round8, by=c("run_name", "sequence_ID"), keep = FALSE)

# clean up data
#complement_data_dissimilar_round8_out <- complement_data_dissimilar_round8_identities[c("counts", "sequence.y", "family_ID", "peak_identity", "complement", "identity")]
#colnames(complement_data_dissimilar_round8_out) <- c("counts", "sequence", "family_ID", "peak_identity", "complement", "identity")
complement_data_dissimilar_round8_out <- complement_data_dissimilar_round8_identities[c("counts", "sequence.y", "family_ID", "peak_identity")]
colnames(complement_data_dissimilar_round8_out) <- c("counts", "sequence", "family_ID", "peak_identity")

# output dissimilar sequence data
write.csv(complement_data_dissimilar_round8_identities, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_round8_identities.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_dissimilar_round8_out, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_round8.csv", sep = ""), row.names = FALSE, quote = FALSE)

# replace 0's through 4's with <5
complement_data[complement_data$tag == 0, "tag"] <- "<5"
complement_data[complement_data$tag == 2, "tag"] <- "<5"
complement_data[complement_data$tag == 3, "tag"] <- "<5"
complement_data[complement_data$tag == 4, "tag"] <- "<5"

# vectors of bins (total, consecutive, gaped)
#tag_bins <- unique(complement_data$tag)
tag_bins <- c("<5", "5", "3_3", "6", "4_3", "3_4", "7", "8")

# set data lengths
data_length <- length(tag_bins)

# list of round IDs
round_list <- unique(complement_data$run_name)

# initialize data frame for identity bin counts
complement_counts_sorted <- data.frame(
  run_name = rep(NA, data_length),
  tag = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length)
)
complement_counts_out <- data.frame()

# loop over each run
for (run_num in min(round_list):max(round_list)) {
  # loop over tag bins
  for (bin_index in 1:data_length) {
    # retrieve current tag
    cur_tag <- tag_bins[bin_index]
    # add run ID
    complement_counts_sorted$run_name[bin_index] <- run_num
    # add tag
    complement_counts_sorted$tag[bin_index] <- cur_tag
    # set the counts run name
    counts_run_num <- paste("r", run_num, "_S", run_num, "_L001", sep = "")
    # set the number of unique sequences
    seq_data_length <- nrow(complement_data[complement_data$run_name == run_num,])
    # add overhang complement counts
    complement_counts_sorted$counts[bin_index] <- sum(complement_data[complement_data$tag == cur_tag & complement_data$run_name == run_num & complement_data$counts_run_name == counts_run_num, "counts"])
    # add overhang complement unique counts
    complement_counts_sorted$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag & complement_data$run_name == run_num & complement_data$counts_run_name == counts_run_num,])
    # add fraction abundance
    complement_counts_sorted$frac_abundance_unique[bin_index] <- complement_counts_sorted$counts_unique[bin_index]/seq_data_length
    complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/quality[run_num]
    #complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/above_two_reads[run_num]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts_sorted)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique
complement_counts_out$perc_abundance <- 100*complement_counts_out$frac_abundance

# export data
write.csv(complement_counts_out, file = paste(out_dir, "/", counts_run_num, "_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)

# loop over each run
#for (run_num in 1:8) {
# calculate fraction abundance
#  complement_counts_sorted[complement_counts_sorted$run_name == run_num, "frac_abundance_unique"] <- complement_counts_sorted[complement_counts_sorted$run_name == run_num, "counts_unique"]/unique_reads[run_num]
#}
# calculate percent abundance
#complement_counts_sorted$perc_abundance_unique <- complement_counts_sorted$frac_abundance_unique*100

# get lists of tags
#tag_list <- unique(complement_counts_out$tag)
tag_list <- c("8", "7", "4_3", "3_4", "6", "3_3", "5", "<5")

# add mapping table
identity_mappings <- data.frame(
  tag = tag_list,
  colors = safe_colors[1:8]
)

# add placeholder columns
complement_counts_out$colors <- NA

# loop over each possible tag label
for (label_num in 1:nrow(identity_mappings)) {
  # set tag colors for plotting
  ifelse(
    complement_counts_out$tag == identity_mappings$tag[label_num], 
    complement_counts_out[complement_counts_out$tag == identity_mappings$tag[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# setup data for plotting
complement_counts_out$tag <- as.factor(complement_counts_out$tag)
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 3,]
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 4,]

# sort the data for plotting
complement_counts_total <-  complement_counts_out %>% arrange(factor(tag, levels = tag_list))
complement_counts_total$tag <- factor(complement_counts_total$tag, levels = c("8", "7", "4_3", "3_4", "6", "3_3", "5", "<5"))

## plots using sequencing read counts

# create line plot of total overhang identity percent
base_abun_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance, group=tag, color=colors))+
  geom_line(size = 1.25) +
  geom_point(size = 2.75) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Complementarity", labels = complement_counts_total$tag, breaks = complement_counts_total$colors, guide = "legend") +
  scale_y_continuous(limits=c(0, 40), breaks=seq(0, 40, 10)) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_abun_plot)
dev.off()

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=tag, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$tag), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  labs(fill = "Complementarity") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=tag, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(fill = 'black')) +
  scale_fill_manual(breaks = unique(complement_counts_total$tag), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  geom_text(aes(label=paste0(sprintf("%1.1f", perc_abundance_unique),"%")),
            position=position_stack(vjust=0.5), size = 3, color = "white") +
  labs(fill = "Complementarity") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart_perc.png", sep = "")
png(exportFile, units="in", width=6, height=4, res=300)
print(base_counts_plot)
dev.off()

# export data
write.csv(complement_counts_sorted, file = paste(out_dir, "/overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_total, file = paste(out_dir, "/overhang_conservation_wobble_total.csv", sep = ""), row.names = FALSE, quote = FALSE)
