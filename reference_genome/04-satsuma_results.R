# look at satsuma outputs

library(tidyverse)
library(stringr)

# set directory that contains all the satsuma_summary.chained.out files (assuming split by chr)
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/satsuma_outputs")

# make list of filenames to import
filenames <- list.files(pattern="satsuma_summary.chained*")

# import them all together
files <- lapply(filenames,function(i){
  read.delim(i, header=FALSE)
})

# merge files from a large list format into a dataframe
satsuma <- do.call(rbind, files)

# rename columns
satsuma <- rename(satsuma,target_chr=V1,start_target_base=V2,end_target_base=V3,query_scaf=V4,start_query_base=V5, end_query_base=V6, identity=V7, orientation=V8)

# also add columns for bp in each row/chunk
satsuma <- mutate(satsuma,"target_bp_spanned"=end_target_base-start_target_base)
satsuma <- mutate(satsuma,"query_bp_spanned"=end_query_base-start_query_base)

# change target sequence name (from "NC..." to chr#). I believe this will go in order of chr1, chr10, chr11, ... then chr2, chr20, chr21, chr3..), check this before running things below.

# extract "NC..." name from each chrom file
target_chr_list <- as.data.frame(unique(satsuma$target_chr))
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
satsuma2 <- satsuma %>% mutate(target_chr=ifelse(as.character(target_chr) == c1, "chr1", as.character(target_chr)))
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
hist(satsuma$identity)

# see how many unique scaffolds on query genome
uniq_scaf <- as.data.frame(unique(satsuma$query_scaf))

# will need to see if these scaffolds hit other chromosomes