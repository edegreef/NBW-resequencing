# look at satsuma outputs and create a list of the query genome scaffolds that make up each target/model chromosome

library(tidyverse)
library(stringr)
library(dplyr)

# set directory that contains all the satsuma_summary.chained.out files (assuming split by chr)
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/satsuma_outputs")

# make list of filenames to import
filenames <- list.files(pattern="satsuma_summary.chained*")

# import them all together
files <- lapply(filenames,function(i){
  read.delim(i, header=FALSE)
})

# merge files from a large list format into a dataframe
satsuma_raw <- do.call(rbind, files)

# rename columns
satsuma_raw <- rename(satsuma_raw,target_chr=V1,start_target_base=V2,end_target_base=V3,query_scaf=V4,start_query_base=V5, end_query_base=V6, identity=V7, orientation=V8)

# also add columns for bp in each row/chunk
satsuma_raw <- mutate(satsuma_raw,"target_bp_spanned"=end_target_base-start_target_base)
satsuma_raw <- mutate(satsuma_raw,"query_bp_spanned"=end_query_base-start_query_base)

# change target sequence name (from "NC..." to chr#). I believe this will go in order of chr1, chr10, chr11, ... then chr2, chr20, chr21, chr3..), check this before running things below.

# extract "NC..." name from each chrom file
target_chr_list <- as.data.frame(unique(satsuma_raw$target_chr))
c1 <- target_chr_list[1,1]
c10 <- target_chr_list[2,1]
c11 <- target_chr_list[3,1]
c12 <- target_chr_list[4,1]
c13 <- target_chr_list[5,1]
c14 <- target_chr_list[6,1]
c15 <- target_chr_list[7,1]
c16 <- target_chr_list[8,1]
c17 <- target_chr_list[9,1]
c18 <- target_chr_list[10,1]
c19 <- target_chr_list[11,1]
c2 <- target_chr_list[12,1]
c20 <- target_chr_list[13,1]
c21 <- target_chr_list[14,1]
c3 <- target_chr_list[15,1]
c4 <- target_chr_list[16,1]
c5 <- target_chr_list[17,1]
c6 <- target_chr_list[18,1]
c7 <- target_chr_list[19,1]
c8 <- target_chr_list[20,1]
c9 <- target_chr_list[21,1]
cx <- target_chr_list[22,1]
cy <- target_chr_list[23,1]

# replace "NC..." with "chr#"
satsuma2 <- satsuma_raw %>% mutate(target_chr=ifelse(as.character(target_chr) == c1, "chr1", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c2, "chr2", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c3, "chr3", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c4, "chr4", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c5, "chr5", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c6, "chr6", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c7, "chr7", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c8, "chr8", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c9, "chr9", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c10, "chr10", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c11, "chr11", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c12, "chr12", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c13, "chr13", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c14, "chr14", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c15, "chr15", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c16, "chr16", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c17, "chr17", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c18, "chr18", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c19, "chr19", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c20, "chr20", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == c21, "chr21", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == cx, "chrx", as.character(target_chr)))
satsuma2 <- satsuma2 %>% mutate(target_chr=ifelse(as.character(target_chr) == cy, "chry", as.character(target_chr)))

# quick look at identity
hist(satsuma2$identity, breaks=50)

# filter out ones with identity less than 0.6 (keeping ones above 0.6)
satsuma3 <- satsuma2 %>% filter(identity > 0.6) 
hist(satsuma3$identity, breaks=50)

# next, see if there are scaffolds double-dipping in multiple chrs. Referred to here as 'duplicates'. Want to keep the ones with higher % match and discard the duplicate scaffold hit with very little match

# create unique ID for chr$scaf info
satsuma3$id <- paste(satsuma3$target_chr, satsuma3$query_scaf, sep="_")

# look at total query_bp_spanned per unique 'id'
scaf_sum <- satsuma3 %>% group_by(id) %>%  summarise(query_bp_spanned = sum(query_bp_spanned))

# list of unique chr+scaf ID
scaf_uniqueid <- satsuma3[!duplicated(satsuma3$id), ]

# pull out unique id and scaf number
scaf_number <- scaf_uniqueid[,c(4,11)]
scaf_sum <- scaf_sum %>% left_join(scaf_number, by=c("id"="id")) %>% na.omit()

# add actual scaffold length to see percent
scaffold_lengths <- read_tsv("C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/Northern_bottlenose_whale_051018_shortLabel.fasta.fai",col_names = c("scaf","scaffold_length"))
scaf_sum <- scaf_sum %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()
scaf_sum$percent <- scaf_sum$query_bp_spanned/scaf_sum$scaffold_length

# then see which are 'duplicate' query_scafs
scaf_dups <- scaf_sum[duplicated(scaf_sum$query_scaf)|duplicated(scaf_sum$query_scaf, fromLast=TRUE), ]

# 206 obs -> looks like there are about 100 duplicates then
# want to keep the one with higher percent (covering more bp in matched chr)
# will be filtering for % cover later on full set, but should first ID which duplicate scaffold entries to remove

# order duplicates by query_scaf then percent in each query_scaf (higher percent listed first)
scaf_dups <- scaf_dups[order(scaf_dups$query_scaf, -scaf_dups$percent),]

# make list of the chr_scaf id to remove
scaf_dups_remove <- scaf_dups[duplicated(scaf_dups$query_scaf),]              
id_remove <- as.data.frame(scaf_dups_remove$id)
colnames(id_remove) <- "id"

# from satsuma3 dataframe, remove rows with the duplicate scaffolds listed in id_remove
satsuma4 <- anti_join(satsuma3, id_remove, by=c("id"))
 
# check that it worked, no dups should be there now. 
scaf_sum_check <- satsuma4 %>% group_by(id) %>%  summarise(query_bp_spanned = sum(query_bp_spanned))
scaf_uniqueid_check <- satsuma4[!duplicated(satsuma4$id), ]
scaf_number_check <- scaf_uniqueid_check[,c(4,11)]
scaf_sum_check <- scaf_sum_check %>% left_join(scaf_number_check, by=c("id"="id")) %>% na.omit()
scaf_sum_check <- scaf_sum_check %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()
scaf_sum_check <- scaf_sum_check[duplicated(scaf_sum_check$query_scaf)|duplicated(scaf_sum_check$query_scaf, fromLast=TRUE), ]

# scaf_sum_check should have 0 obs if 'duplicate' removal worked. 

# moving on to making usable scaf/chr list

# see how many bp spanned total for each query scaffold and see how many from full scaffold length. probably will need to filter some of these scaffolds out that have little portion actually matched
group_scaf <- satsuma4 %>% group_by(query_scaf) %>%  summarise(query_bp_spanned = sum(query_bp_spanned))

# add in scaffold length info for looking at genome proportion
group_scaf <- group_scaf %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()

# check % covered
group_scaf$percent <- group_scaf$query_bp_spanned/group_scaf$scaffold_length
hist(group_scaf$percent)

# some are above 100%, let's see what is going on
# looking at just chr1 first
chrom1 <- subset(satsuma4, target_chr=="chr1")
chrom1_scaf <- chrom1 %>% group_by(query_scaf) %>%  summarise(query_bp_spanned = sum(query_bp_spanned))
chrom1_scaf <- chrom1_scaf %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()
chrom1_scaf$percent <- chrom1_scaf$query_bp_spanned/chrom1_scaf$scaffold_length
hist(chrom1_scaf$percent)

# pull out a scaffold that has over 100%
test<- subset(satsuma3, satsuma3$query_scaf==232488)
sum(test$query_bp_spanned)

# scaf 232488 should have 39,936 bp but has 59,913 covered here. Looking at start & end query base, looks like some overlap within chr. So that's why some of the summaries are over 100%. This is probably ok, since we just want to filter out the low % matches.

# filter out scaffolds that didn't hit at least at least half synteny since these may not be very certain
group_scaf_filtered <- group_scaf %>% filter(percent > 0.5)

# see how much of NBW genome is retained
NBW_ref=2353367117
kept <- sum(group_scaf_filtered$scaffold_length)
kept/NBW_ref

# looks like 93% retained, meaning 7% of NBW genome scaffolds from satsuma are unlocalized

# make list of each scaffold and it's respective chromosome, first need to subset satsuma results with the scaffolds we want to keep.
satsuma_filtered <- group_scaf_filtered %>% full_join(satsuma4, by=c("query_scaf"="query_scaf")) %>% na.omit()

# pulling out chr/scaf matches, then selecting specific columns
scaf_list_prep <- satsuma_filtered[!duplicated(satsuma_filtered$query_scaf),] %>% select("query_scaf", "scaffold_length", "target_chr", "start_target_base", "orientation")

# order rows by target_chr, then by target_position
scaf_list_prep <- scaf_list_prep[order(scaf_list_prep$target_chr, scaf_list_prep$start_target_base),]

# need to order with manual set since otherwise it'll go by chr1, chr10, chr11:
good_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chrx", "chry")

# proper chr order
scaf_list <- scaf_list_prep %>% arrange(factor(target_chr, levels=good_order))

# noticing that Y is gone (probably since filter percent > 60%). just check this to make sure that is indeed what happened:
check_y <- subset(satsuma3, target_chr=="chry")
check_y <- check_y %>% group_by(query_scaf) %>% summarise(query_bp_spanned = sum(query_bp_spanned)) %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()
check_y$percent <- check_y$query_bp_spanned/check_y$scaffold_length
hist(check_y$percent)

# highest percent in Y was 41%. Since we determined X and Y separately before with coverage things, we already know which scaffolds hit these chrs.

# could use the scaf_list as final list, but I want to see if there are any scaffolds double dipping with X and Y results from difcover, so planning to use satsuma results here for autosomal chrs

# check any double dipping scaffolds on X and Y 

# load lists of X and Y
X_scafs <- read.table("C:/Users/Evelien de Greef/Dropbox/NBW-me/MFcoverage/series/final_scaffolds_Xlinked.txt", header=FALSE, col.names="query_scaf")
Y_scafs <- read.table("C:/Users/Evelien de Greef/Dropbox/NBW-me/MFcoverage/series/final_scaffolds_Ylinked.txt", header=FALSE, col.names="query_scaf")

# subsetting scaf_list for no x chrs
scaf_list_no_x <- subset(scaf_list, target_chr != "chrx")

X_test <- X_scafs %>% left_join(scaf_list_no_x, by="query_scaf") %>% na.omit()
# only 2 scaffolds double dipping here

Y_test <- Y_scafs %>% left_join(scaf_list_no_x, by="query_scaf") %>% na.omit()
# no double dips here

# remove the 2 X-chr double dips from the scaf_list. Since these two are small scaffolds, I wouldn't consider these two an issue. I've already determined X and Y with coverage-based method so going to keep these 2 scafs as "X" and then remove the duplicate from scaf_list
# prep ID for removal
X_dip <- as.data.frame(X_test$query_scaf)
colnames(X_dip) <- "query_scaf"
scaf_list2 <- anti_join(scaf_list, X_dip, by=c("query_scaf"))

# removing satsuma-determined X chr because will be using the difcover results for X in the end. The downside here is that the ordering of X and Y will not be included
scaf_list3 <- subset(scaf_list2, target_chr != "chrx")

# can save a file here for just satsuma-autosomes
write.csv(scaf_list3, "nbw_scaf_chr_info.csv", row.names = FALSE)

# now for making 'master' list

# add an ordering id
scaf_list3$order_id <- 1:nrow(scaf_list3)

# merge sex-scaf info
X_scafs$target_chr <- "chrx"
Y_scafs$target_chr <- "chry"
sex_scafs <- rbind(X_scafs, Y_scafs)
sex_scafs <- sex_scafs %>% left_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()

# add empty columns in sex_scaf list, just as a placeholder so we can rbind with scaf_list
sex_scafs$orientation <- NA
sex_scafs$start_target_base <- NA
sex_scafs <- sex_scafs %>% select("query_scaf", "scaffold_length", "target_chr", "start_target_base", "orientation")
sex_scafs$order_id <- NA

# make a final list that has info for scaffolds belonging to each chrs with X and Y at the end
scaf_list_final <- rbind(scaf_list3, sex_scafs)

# save list for nbw things!
write.csv(scaf_list_final, "nbw_scaf_chr_info_allchrs.csv", row.names = FALSE)


### comparing total chromosome lengths in bottlenose whale with blue whale model
# prep file
chr_lengths <- scaf_list_final %>% group_by(target_chr) %>%  summarise(scaffold_length = sum(scaffold_length)) %>% arrange(factor(target_chr, levels=good_order))
colnames(chr_lengths) <- c("chr", "nbw_length")

# getting lengths for unplaced scaffolds
unplaced <- scaffold_lengths %>% 
  anti_join(scaf_list_final, by=c("scaf"="query_scaf")) %>% na.omit()

unplaced_length <- sum(unplaced$scaffold_length)
unplaced_length / NBW_ref
# 6% unplaced, not too bad.

# add these to the chr_length
chr_temp <- data.frame(first_column = "unplaced",
                  second_column = unplaced_length)
colnames(chr_temp) <- c("chr", "nbw_length")
chr_lengths <- rbind(chr_lengths, chr_temp)

# total sum should be same as total ref genome length
sum(chr_lengths$nbw_length) == NBW_ref

# pulling blue whale chr length from ncbi (https://www.ncbi.nlm.nih.gov/assembly/GCF_009873245.2/#/st)
chr_lengths$blw_length <- c(184938300, #chr1
                            175897734, #chr2
                            171266408, #chr3
                            144968589, #chr4
                            140689829, #chr5
                            116510015, #chr6
                            113414938, #chr7
                            110314666, #chr8
                            107421550, #chr9
                            104744437, #chr10
                            104068540, #chr11
                            91445419, #chr12
                            90635089, #chr13
                            90457838, #chr14
                            88470553, #chr15
                            86152963, #chr16
                            81207215, #chr17
                            79663398, #chr18
                            60735208, #chr19
                            60304989, #chr20
                            36241783, #chr21
                            128877148, #chrx,
                            2349494, #chry
                            4076438 #un
                            )

# total blw sum should be same as total blw genome length
BLW_ref=2374852541
sum(chr_lengths$blw_length) == BLW_ref

# add genome proportions too? why not
chr_lengths$nbw_proportion <- chr_lengths$nbw_length/NBW_ref
chr_lengths$blw_proportion <- chr_lengths$blw_length/BLW_ref

# save chr list. Just be aware the the X and Y are specifically from difcover analyses, so this list is reflecting overall comparison with nbw and blw, but not solely satsuma resuts.
write.csv(chr_lengths, "NBW_BLW_compare_chrlengths.csv", row.names=FALSE)

# plotting to compare for fun
library(ggplot2)
library(ggExtra)
library(reshape2)
library(ggrepel)

chr_lengths <- read.csv("NBW_BLW_compare_chrlengths.csv", header=TRUE)

plot <- ggplot(chr_lengths, aes(nbw_proportion, blw_proportion)) +
  geom_abline(slope=1, intercept=0, linetype="dotted", col="gray50")+
  geom_point(size=2.5, col="#00008b", alpha=0.8, pch=16) + 
  theme_bw() +
  geom_label_repel(aes(label=chr), max.overlaps = 20)+
  xlab("NBW genome proportion")+
  ylab("BLW genome proportion")

plot_m <- ggExtra::ggMarginal(
  plot,
  type = 'histogram',
  margins = 'both',
  xparams = list(bins=30),
  yparams = list(bins=30),
  size = 10,
  colour = 'gray20',
  fill = '#6699cc'
)

plot_m
ggsave(filename="NBW_BLW_compare_plot.png", plot=plot_m,width=9,height=7,dpi=2000)
