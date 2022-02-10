# Looking at DifCover output (.DNAcopyout) using male and female samples to determine X and Y linked scaffolds. Doing this for 4 male v female runs (M1F1, M1F2, M2F1, M2F2) and 2 controls (M1M2, F1F2). Script for the 2 controls towards the bottom.
# Guidance and code help from Phil Grayson (https://github.com/phil-grayson/SexFindR), with some modifications

library(tidyverse)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/MFcoverage/series")

# enter in file infos
DNAcopyout ="M1F1sample1_sample2.ratio_per_w_CC0_a3_A30_b3_B27_v1000_l500.log2adj_0.9.DNAcopyout"
type="M1F1"
scaffold_info="Northern_bottlenose_whale_051018_shortLabel.fasta.fai"

# load in the raw data (including non-significant regions), renaming columns, and adding bp spanned
difcover <- read_tsv(file = DNAcopyout ,col_names = F) %>% rename(scaf=X1,start=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases spanned" = stop-start)

# parse samtools faidx output for scaffold name and length
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","length"))

# join the two 
proportion <- full_join(difcover,scaffold_lengths) %>% mutate(proportion = `bases spanned`/length)

# plot proportion of scaffold versus log2(male/female) coverage
proportion %>% ggplot(aes(x=`log2(Male coverage/Female coverage)`,y=proportion)) + geom_point()+xlab("log2(Male coverage/Female coverage)")#+geom_vline(xintercept=-0.5, col="white")
#ggsave(paste("coverage.", type, ".png",sep=""), width = 11, height = 8, dpi = 300)

# interested in difCover regions that are likely to be Y (enrichement greater than 2) or X (enrichment less than -0.7369656)
filtered_proportion <- proportion %>% filter(`log2(Male coverage/Female coverage)` >= 2 | `log2(Male coverage/Female coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

# plot proportions to see how many are substantial proportions of their scaffold; based on plot, may want to only consider scaffolds that fit a cutoff like 0.25 (set as red line here)
filtered_proportion %>% ggplot(aes(x=scaf,y=`total chromosome proportion with significantly different coverage`)) + geom_point(alpha=0.2) + geom_hline(yintercept=0.25, colour="red") + ggtitle(paste("total scaf proportion with sig different coverage -", type, sep=" ")) + ylab("total proportion") + xlab("scaffold")
ggsave(paste("scaffoldproportion_sig_", type, ".png",sep=""), width = 7, height = 5, dpi = 300)
# +geom_text(aes(label=scaf), alpha=0.5) 

# keeping only scaffolds that have proportion of at least 0.25
filtered_proportion %>% filter(`total chromosome proportion with significantly different coverage` >= 0.25) %>% ggplot(aes(x=scaf,y=`total chromosome proportion with significantly different coverage`)) + geom_point(alpha=0.2)

# number of scaffolds remaining using 0.25 cutoff
likely_sex_linked <- filtered_proportion %>% filter(`total chromosome proportion with significantly different coverage` >= 0.25)
nrow(likely_sex_linked)

# number of scaffolds in X and Y
likely_X_linked <- likely_sex_linked %>% filter(`log2(Male coverage/Female coverage)` <= -0.7369656)
nrow(likely_X_linked)

likely_Y_linked <- likely_sex_linked %>% filter(`log2(Male coverage/Female coverage)` >= 2)
nrow(likely_Y_linked)

# number of bases in X and Y
X_length <- sum(likely_X_linked$`bases spanned`)
X_length

Y_length <- sum(likely_Y_linked$`bases spanned`)
Y_length

# proportion of X and Y in whole genome
genome_length <- sum(scaffold_lengths$length)
X_length/genome_length #something around 4-6% should be good for this species
Y_length/genome_length #Y is expected to be very small

# give the X and Y linked scaffolds appropriate labels
likely_X_linked_list <- likely_X_linked %>% mutate(chr_type="X") %>% select(scaf,chr_type)
likely_Y_linked_list <- likely_Y_linked %>% mutate(chr_type="Y") %>% select(scaf,chr_type)
linked_with_names <- bind_rows(likely_X_linked_list,likely_Y_linked_list)

# combine that with the original "proportion" tibble and if there isn't a chr_type then label as autosome
proportion_with_names <- left_join(proportion,linked_with_names) %>% replace_na(list(chr_type = "Autosome"))

# plot proportion of scaffold versus log2(male/female) coverage with color codes
proportion_with_names %>% ggplot(aes(x=`log2(Male coverage/Female coverage)`,y=proportion)) + geom_point(aes(color=chr_type), alpha=0.3) +xlab("log2(Male coverage/Female coverage)") + theme_classic() + ggtitle(paste("coverage of scaffolds between male and female sample -", type, sep=" "))# + geom_vline(xintercept =0)
ggsave(paste("Coverage plot_", type, ".png",sep=""), width = 10, height = 6, dpi = 300)

# prepare data for making bed files from the likely_X/Y_linked data (to use later for bedtools annotate)
likely_X_linked$sexchr <- "X"
forbed_X <- likely_X_linked %>% select(scaf, start, stop, sexchr)
write.table(forbed_X, paste("X_linked_scafwindows.", type, ".bed",sep=""), sep="\t", row.names=FALSE,col.names=FALSE)

likely_Y_linked$sexchr <- "Y"
forbed_Y <- likely_Y_linked %>% select(scaf, start, stop, sexchr)
write.table(forbed_Y, paste("Y_linked_scafwindows.", type, ".bed",sep=""), sep="\t", row.names=FALSE,col.names=FALSE)



######## looking at same-sex control data (M1M2, F1F2):

# enter file infos
DNAcopyout ="M1M2_sample1_sample2.ratio_per_w_CC0_a3_A30_b3_B27_v1000_l500.log2adj_0.9.DNAcopyout"
type="M1M2"
scaffold_info="Northern_bottlenose_whale_051018_shortLabel.fasta.fai"

# load in the raw data (including non-significant regions), renaming columns, and adding bp spanned
difcover <- read_tsv(file = DNAcopyout ,col_names = F) %>% rename(scaf=X1,start=X2,stop=X3,windows=X4,"log2(samsex coverage/samsex coverage)"=X5) %>% mutate("bases spanned" = stop-start)

# parse samtools faidx output for scaffold name and length & join with difcover
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","length"))
proportion <- full_join(difcover,scaffold_lengths) %>% mutate(proportion = `bases spanned`/length)

# plot proportion of scaffold versus log2(male/female) coverage
proportion %>% ggplot(aes(x=`log2(samsex coverage/samsex coverage)`,y=proportion)) + geom_point()+xlab("log2(Female coverage/Female coverage)")#+geom_vline(xintercept=0, col="white")
ggsave(paste("coverage.", type, ".png",sep=""), width = 11, height = 8, dpi = 300)

# to keep not-enriched sites:
filtered_proportion_not_enrich <- proportion %>% filter(`log2(samsex coverage/samsex coverage)` < 2 & `log2(samsex coverage/samsex coverage)` > -2) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

# plotting not-enriched scafs
filtered_proportion_not_enrich %>% ggplot(aes(x=scaf,y=`total chromosome proportion with significantly different coverage`)) + geom_point(alpha=0.2) + geom_hline(yintercept=0.25, colour="red") + ggtitle(paste("total scaf proportion without sig different coverage -", type, sep=" ")) + ylab("total proportion") + xlab("scaffold")
ggsave(paste("scaffoldproportion_notenrich_", type, ".png",sep=""), width = 7, height = 5, dpi = 300)

# filtering scaffolds with 0.25 cutoff
likely_not_enriched <- filtered_proportion_not_enrich %>% filter(`total chromosome proportion with significantly different coverage` >= 0.25)
nrow(likely_not_enriched)
notenrich_length <- sum(likely_not_enriched$`bases spanned`)
notenrich_length
genome_length <- sum(scaffold_lengths$length)
notenrich_length/genome_length #should be high proportion, here was 95%

likely_not_enriched_list <- likely_not_enriched %>% mutate(chr_type="Control") %>% select(scaf,chr_type)
proportion_with_names <- left_join(proportion,likely_not_enriched_list) %>% replace_na(list(chr_type = "Outlier"))
proportion_with_names %>% ggplot(aes(x=`log2(samsex coverage/samsex coverage)`,y=proportion)) + geom_point(aes(color=chr_type), alpha=0.3) +xlab("log2(Female coverage/Female coverage)") + theme_classic() + ggtitle(paste("coverage of scaffolds between male and male sample -", type, sep=" ")) #+ geom_vline(xintercept =0)
ggsave(paste("Coverage plot_control_", type, ".png",sep=""), width = 9, height = 5, dpi = 300)

# prepare data for making bed file (to use later for bedtools annotate)
likely_not_enriched$chr <- "Control"
likely_not_enriched$rownum  <- 1:nrow(likely_not_enriched)
forbed_control <- likely_not_enriched %>% select(scaf, start, stop, chr, rownum)

# before saving bed file for the control runs, check that no window "stop"s are at 0, otherwise bedtools will give error. Check the "likely_not_enriched" dataframe and sort the "stop" column in order by low-> high and then if there's a scaffold window with start and stop both at 0, then remove the row.

# F1F2: scaffold 132977 has a window 0 to 0 but also has 3 more windows in next line. Row 3515 in "forbed_control" should be removed (0 to 0). Not sure why this row exists, but I suspect this extra empty row happened in the difcover run.
# M1M2: scaffold 204716 same thing. Row 21806

#forbed_control <- forbed_control[-3515,-5] 
forbed_control <- forbed_control[-21806,-5]
write.table(forbed_control, paste("Control_scafwindows.", type, ".bed",sep=""), sep="\t", row.names=FALSE,col.names=FALSE)
