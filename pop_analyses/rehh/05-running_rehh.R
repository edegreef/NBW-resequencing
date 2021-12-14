# Using rehh for haplotype analyses to identify regions under selection
# Coding help from Matt Thorstensen

library(rehh)
library(vcfR)
library(tidyverse)
library(ggplot2)
library(DescTools)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

# create map file
vcf <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.imputed.vcf.gz")
chroms <- unique(vcf@fix[,1])
colnames(vcf@fix)
# Create a map file from the vcf for use with rehh
map_file <- tibble(marker_ID = vcf@fix[,3], chromosome = vcf@fix[,1], position = vcf@fix[,2], ref = vcf@fix[,4], alt = vcf@fix[,5])
#write_delim(map_file, "rehh_mapfile_chr10.inp", delim = "\t", col_names = FALSE)
#rm(vcf)

map_file$marker_ID <- paste(map_file$chromosome, map_file$position, sep="-")


# #### did this loop step for each chr separately and for each group, did on biology-01 computer cuz it took a while, then reloaded inputs further below.
#for(i in 1:length(chroms)) {
# create internal representation
#  hh <- data2haplohh(hap_file = "NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss01.n36.notinyscaf.imputed.mac2.min100kb.northern.recode.chr1.vcf",
#                     map_file = "rehh_mapfile_chr1_only.inp", 
#                     chr.name = chroms[i],
#                     polarize_vcf = FALSE,
#                     vcf_reader = "vcfR")
# perform scan on a single chromosome (calculate iHH values)
#  scan <- scan_hh(hh, polarized = FALSE)
# concatenate chromosome-wise data frames to
# a data frame for the whole genome
# (more efficient ways certainly exist...)
#  if (i == 1) {
#    NOR_chr1_wgscan <- scan
#  } else {
#    NOR_chr1_wgscan <- rbind(NOR_chr1_wgscan, scan)
#  }
#}

#write.csv(NOR_chr1_wgscan, "NOR_chr1_wgscan.csv")


# load loop output scan file, chr1-21 and chru (unlocalized)
# make list of filenames to import
setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh/wgscan")

# Northern data
NOR_filenames <- list.files(pattern="NOR*") # make filename pattern
NOR_files <- lapply(NOR_filenames,function(i){
  read.csv(i, header=TRUE, row.names=1)
}) # load all files
NOR_scan_merged <- do.call(rbind, NOR_files) # merge files

# Scotian Shelf data
SS_filenames <- list.files(pattern="SS*") # make filename pattern
SS_files <- lapply(SS_filenames,function(i){
  read.csv(i, header=TRUE, row.names=1)
}) # load all files
SS_scan_merged <- do.call(rbind, SS_files) # merge files



# all NBW as one pop
#ALL_filenames <- list.files(pattern="ALL*") # make filename pattern
#ALL_files <- lapply(ALL_filenames,function(i){
#  read.csv(i, header=TRUE, row.names=1)
#}) # load all files
#ALL_scan_merged <- do.call(rbind, ALL_files) # merge files

# swtich back to previous directory
setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

####################
# Calculate the XP-EHH statistic between each pair of populations (Northern vs Scotian Shelf)
xpehh.NOR_SS <- ies2xpehh(scan_pop1 = NOR_scan_merged,
                          scan_pop2 = SS_scan_merged,
                          popname1 = "Northern",
                          popname2 = "ScotianShelf",
                          p.adjust.method = "fdr")

# add snp ID
xpehh.NOR_SS$SNP <- paste(xpehh.NOR_SS$CHR, xpehh.NOR_SS$POSITION, sep="-")

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.01. -log10(q) > 2 because -log10(0.01) = 2
# q < 0.05  (-log10(q) > 1.30103
NORSS_candidates <- calc_candidate_regions(xpehh.NOR_SS, 
                                           threshold = 1.30103, 
                                           pval = TRUE, 
                                           window_size = 1E5, 
                                           overlap = 1E4, 
                                           min_n_extr_mrk = 2)
NORSS_candidates

# add chromosome info to xpehh.NOR_SS  (for plotting)
scafs  <- read.csv("C:/Users/eveli/Dropbox/NBW-me/reference/satsuma_outputs/nbw_scaf_chr_info_autosomes.csv")
scafs$contig <- paste("contig", scafs$query_scaf, sep="")

# add chr info to xpehh
xpehh.NOR_SS.chr <- xpehh.NOR_SS %>% left_join(scafs, by=c("CHR"="contig"))
xpehh.NOR_SS.chr$target_chr[is.na(xpehh.NOR_SS.chr$target_chr)] <- "chru"

# adding new position number for within chromosome
xpehh.NOR_SS.chr <- xpehh.NOR_SS.chr %>% group_by(target_chr) %>% mutate(new_pos = row_number())

# order chr1, chr2... instead of chr1, chr10..
good_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chru")

# proper chr order
xpehh.NOR_SS.chr <- xpehh.NOR_SS.chr %>% arrange(factor(target_chr, levels=good_order))
xpehh.NOR_SS.chr$plotorder <- 1:nrow(xpehh.NOR_SS.chr)

# Filter for significant hits (q < 0.01) among the NOR_SS XP-EHH SNPs
# q < 0.05, 1.30103
sig_xpehh.NOR_SS <- filter(xpehh.NOR_SS.chr, LOGPVALUE > 1.30103)

# Check which SNPs appear within a candidate region
sig_xpehh.NOR_SS$cand_region <- ifelse(sapply(sig_xpehh.NOR_SS$POSITION, function(p) 
  any(NORSS_candidates$START <= p & NORSS_candidates$END >= p)),"cand", NA)

# remove single marks
sig_xpehh.NOR_SS2<- sig_xpehh.NOR_SS[sig_xpehh.NOR_SS$CHR %in% sig_xpehh.NOR_SS$CHR[duplicated(sig_xpehh.NOR_SS$CHR)],]

sig_xpehh.NOR_SS$cand_region <- ifelse(sapply(seq_along(sig_xpehh.NOR_SS$POSITION), function(i) {
  inds <- NORSS_candidates$START <= sig_xpehh.NOR_SS$POSITION[i] & NORSS_candidates$END >= sig_xpehh.NOR_SS$POSITION[i]
  any(inds) & (sig_xpehh.NOR_SS$CHR[i] == NORSS_candidates$CHR[which.max(inds)])
}), "cand", NA)

# Pull positive (selection in Northern) and negative (selection in Scotian Shelf) SNPs for highlighting on the manhattan plot
xpehh_NORSS_pos <- filter(sig_xpehh.NOR_SS2, LOGPVALUE > 1.30103 & `XPEHH_Northern_ScotianShelf` > 0 & cand_region == "cand")

xpehh_NORSS_neg <- filter(sig_xpehh.NOR_SS2, LOGPVALUE > 1.30103 & `XPEHH_Northern_ScotianShelf` < 0 & cand_region == "cand")

# prepare plot
chr_colors <- c(chr1="gray60",
                chr2="gray70", 
                chr3="gray60", 
                chr4="gray70", 
                chr5="gray60",
                chr6="gray70",
                chr7="gray60",
                chr8="gray70",
                chr9="gray60",
                chr10="gray70",
                chr11="gray60",
                chr12="gray70",
                chr13="gray60",
                chr14="gray70",
                chr15="gray60",
                chr16="gray70",
                chr17="gray60",
                chr18="gray70",
                chr19="gray60",
                chr20="gray70",
                chr21="gray60",
                chru="gray70")

# work on split later
# regular gray dot size 1.5, and bold dot size 2.5
# changing for squish
# intercept +- 4.56 for log 2, +- 4.1302 for log 1.30103
xpehh_NORSS_plot <- ggplot()+
  geom_point(data=xpehh.NOR_SS.chr, aes(x=plotorder,y=XPEHH_Northern_ScotianShelf,color=target_chr), alpha=0.5, size=0.7)+
  scale_color_manual(values=chr_colors)+
 geom_point(data=xpehh_NORSS_pos, aes(x=plotorder, y=XPEHH_Northern_ScotianShelf), size=1.5,colour="#313695") +
 geom_point(data=xpehh_NORSS_neg, aes(x=plotorder, y=XPEHH_Northern_ScotianShelf), size=1.5,colour="#A50026") +
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
        #panel.grid.minor.y=element_blank(),
        #panel.grid.major.y=element_blank(),
        text=element_text(family="serif", size=12),axis.text.x=element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  ylab("XP-EHH")+
  ylim(c(-9.5,9.5))

xpehh_NORSS_plot
ggsave("rehh_NBW2_2M_ggplot_NORSS_log1.3_darkblue.png", width=10,height=2.5,dpi=1000)

write.csv(xpehh_NORSS_pos, "sigsnps_xpehh_northern_logp1.3.csv")
write.csv(xpehh_NORSS_neg, "sigsnps_xpehh_scotianshelf_logp1.3.csv")

##################################################################                                  
###################### Did same thing for other group comparisons:
##################################################################
                                              
# create map file
vcf <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.imputed.vcf.gz")

chroms <- unique(vcf@fix[,1])
colnames(vcf@fix)
# Create a map file from the vcf for use with rehh
map_file <- tibble(marker_ID = vcf@fix[,3], chromosome = vcf@fix[,1], position = vcf@fix[,2], ref = vcf@fix[,4], alt = vcf@fix[,5])
#write_delim(map_file, "rehh_mapfile_chr10.inp", delim = "\t", col_names = FALSE)
#rm(vcf)

map_file$marker_ID <- paste(map_file$chromosome, map_file$position, sep="-")

# load loop output scan file, chr1-21 and chru (unlocalized)
# make list of filenames to import
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh/wgscan")

# IC data
IC_filenames <- list.files(pattern="IC*") # make filename pattern
IC_files <- lapply(IC_filenames,function(i){
  read.csv(i, header=TRUE, row.names=1)
}) # load all files
IC_scan_merged <- do.call(rbind, IC_files) # merge files

# ARLBNF
ARLBNF_filenames <- list.files(pattern="ARLBNF*") # make filename pattern
ARLBNF_files <- lapply(ARLBNF_filenames,function(i){
  read.csv(i, header=TRUE, row.names=1)
}) # load all files
ARLBNF_scan_merged <- do.call(rbind, ARLBNF_files) # merge files

# Scotian Shelf data
SS_filenames <- list.files(pattern="SS*") # make filename pattern
SS_files <- lapply(SS_filenames,function(i){
  read.csv(i, header=TRUE, row.names=1)
}) # load all files
SS_scan_merged <- do.call(rbind, SS_files) # merge files

# swtich back to previous directory
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

###########
# Calculate the XP-EHH statistic between each pair of populations (Northern vs Scotian Shelf)
xpehh.IC_SS <- ies2xpehh(scan_pop1 = IC_scan_merged,
                          scan_pop2 = SS_scan_merged,
                          popname1 = "Iceland",
                          popname2 = "ScotianShelf",
                          p.adjust.method = "fdr")


# add snp ID
xpehh.IC_SS$SNP <- paste(xpehh.IC_SS$CHR, xpehh.IC_SS$POSITION, sep="-")

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.01. -log10(q) > 2 because -log10(0.01) = 2
# q < 0.05  (-log10(q) > 1.30103
ICSS_candidates <- calc_candidate_regions(xpehh.IC_SS, 
                                           threshold = 1.30103, 
                                           pval = TRUE, 
                                           window_size = 1E5, 
                                           overlap = 1E4, 
                                           min_n_extr_mrk = 2)
ICSS_candidates

# add chromosome info to xpehh.NOR_SS  (for plotting)
scafs  <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/satsuma_outputs/nbw_scaf_chr_info_autosomes.csv")
scafs$contig <- paste("contig", scafs$query_scaf, sep="")

# add chr info to xpehh
xpehh.IC_SS.chr <- xpehh.IC_SS %>% left_join(scafs, by=c("CHR"="contig"))
xpehh.IC_SS.chr$target_chr[is.na(xpehh.IC_SS.chr$target_chr)] <- "chru"

# adding new position number for within chromosome
xpehh.IC_SS.chr <- xpehh.IC_SS.chr %>% group_by(target_chr) %>% mutate(new_pos = row_number())

# order chr1, chr2... instead of chr1, chr10..
good_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chru")

# proper chr order
xpehh.IC_SS.chr <- xpehh.IC_SS.chr %>% arrange(factor(target_chr, levels=good_order))
xpehh.IC_SS.chr$plotorder <- 1:nrow(xpehh.IC_SS.chr)

# Filter for significant hits (q < 0.01) among the NOR_SS XP-EHH SNPs
# q < 0.05, 1.30103
sig_xpehh.IC_SS <- filter(xpehh.IC_SS.chr, LOGPVALUE > 1.30103)

# Check which SNPs appear within a candidate region
sig_xpehh.IC_SS$cand_region <- ifelse(sapply(sig_xpehh.IC_SS$POSITION, function(p) 
  any(ICSS_candidates$START <= p & ICSS_candidates$END >= p)),"cand", NA)

# remove single marks
sig_xpehh.IC_SS2<- sig_xpehh.IC_SS[sig_xpehh.IC_SS$CHR %in% sig_xpehh.IC_SS$CHR[duplicated(sig_xpehh.IC_SS$CHR)],]

sig_xpehh.IC_SS$cand_region <- ifelse(sapply(seq_along(sig_xpehh.IC_SS$POSITION), function(i) {
  inds <- ICSS_candidates$START <= sig_xpehh.IC_SS$POSITION[i] & ICSS_candidates$END >= sig_xpehh.IC_SS$POSITION[i]
  any(inds) & (sig_xpehh.IC_SS$CHR[i] == ICSS_candidates$CHR[which.max(inds)])
}), "cand", NA)

# Pull positive (selection in Northern) and negative (selection in Scotian Shelf) SNPs for highlighting on the manhattan plot
xpehh_ICSS_pos <- filter(sig_xpehh.IC_SS2, LOGPVALUE > 1.30103 & `XPEHH_Iceland_ScotianShelf` > 0 & cand_region == "cand")

xpehh_ICSS_neg <- filter(sig_xpehh.IC_SS2, LOGPVALUE > 1.30103 & `XPEHH_Iceland_ScotianShelf` < 0 & cand_region == "cand")

xpehh_ICSS_plot <- ggplot()+
  geom_point(data=xpehh.IC_SS.chr, aes(x=plotorder,y=XPEHH_Iceland_ScotianShelf,color=target_chr), alpha=0.5, size=0.7)+
  scale_color_manual(values=chr_colors)+
  geom_point(data=xpehh_ICSS_pos, aes(x=plotorder, y=XPEHH_Iceland_ScotianShelf), size=1.5,color="#4575B4")+
  geom_point(data=xpehh_ICSS_neg, aes(x=plotorder, y=XPEHH_Iceland_ScotianShelf), size=1.5,colour="#A50026") +
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
        #panel.grid.minor.y=element_blank(),
        #panel.grid.major.y=element_blank(),
        text=element_text(family="serif", size=12),axis.text.x=element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  ylab("XP-EHH")+
  ylim(c(-9.5,9.5))

xpehh_ICSS_plot
ggsave("rehh_NBW2_2M_ggplot_ICSS_log1.3_test6.png", width=10,height=2.5,dpi=1000)


#------------------------------------ARLBNF and SS
####################
# Calculate the XP-EHH statistic between each pair of populations (Northern vs Scotian Shelf)
xpehh.ARLBNF_SS <- ies2xpehh(scan_pop1 = ARLBNF_scan_merged,
                         scan_pop2 = SS_scan_merged,
                         popname1 = "ARLBNF",
                         popname2 = "ScotianShelf",
                         p.adjust.method = "fdr")


# add snp ID
xpehh.ARLBNF_SS$SNP <- paste(xpehh.ARLBNF_SS$CHR, xpehh.ARLBNF_SS$POSITION, sep="-")

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.01. -log10(q) > 2 because -log10(0.01) = 2
# q < 0.05  (-log10(q) > 1.30103
ARLBNFSS_candidates <- calc_candidate_regions(xpehh.ARLBNF_SS, 
                                          threshold = 1.30103, 
                                          pval = TRUE, 
                                          window_size = 1E5, 
                                          overlap = 1E4, 
                                          min_n_extr_mrk = 2)
ARLBNFSS_candidates

# add chromosome info to xpehh.NOR_SS  (for plotting)
scafs  <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/satsuma_outputs/nbw_scaf_chr_info_autosomes.csv")
scafs$contig <- paste("contig", scafs$query_scaf, sep="")

# add chr info to xpehh
xpehh.ARLBNF_SS.chr <- xpehh.ARLBNF_SS %>% left_join(scafs, by=c("CHR"="contig"))
xpehh.ARLBNF_SS.chr$target_chr[is.na(xpehh.ARLBNF_SS.chr$target_chr)] <- "chru"

# adding new position number for within chromosome
xpehh.ARLBNF_SS.chr <- xpehh.ARLBNF_SS.chr %>% group_by(target_chr) %>% mutate(new_pos = row_number())

# order chr1, chr2... instead of chr1, chr10..
good_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chru")

# proper chr order
xpehh.ARLBNF_SS.chr <- xpehh.ARLBNF_SS.chr %>% arrange(factor(target_chr, levels=good_order))
xpehh.ARLBNF_SS.chr$plotorder <- 1:nrow(xpehh.ARLBNF_SS.chr)

# Filter for significant hits (q < 0.01) among the NOR_SS XP-EHH SNPs
# q < 0.05, 1.30103
sig_xpehh.ARLBNF_SS <- filter(xpehh.ARLBNF_SS.chr, LOGPVALUE > 1.30103)

# Check which SNPs appear within a candidate region
sig_xpehh.ARLBNF_SS$cand_region <- ifelse(sapply(sig_xpehh.ARLBNF_SS$POSITION, function(p) 
  any(ARLBNFSS_candidates$START <= p & ARLBNFSS_candidates$END >= p)),"cand", NA)

# remove single marks
sig_xpehh.ARLBNF_SS2<- sig_xpehh.ARLBNF_SS[sig_xpehh.ARLBNF_SS$CHR %in% sig_xpehh.ARLBNF_SS$CHR[duplicated(sig_xpehh.ARLBNF_SS$CHR)],]

sig_xpehh.ARLBNF_SS$cand_region <- ifelse(sapply(seq_along(sig_xpehh.ARLBNF_SS$POSITION), function(i) {
  inds <- ARLBNFSS_candidates$START <= sig_xpehh.ARLBNF_SS$POSITION[i] & ARLBNFSS_candidates$END >= sig_xpehh.ARLBNF_SS$POSITION[i]
  any(inds) & (sig_xpehh.ARLBNF_SS$CHR[i] == ARLBNFSS_candidates$CHR[which.max(inds)])
}), "cand", NA)

# Pull positive (selection in Northern) and negative (selection in Scotian Shelf) SNPs for highlighting on the manhattan plot
xpehh_ARLBNFSS_pos <- filter(sig_xpehh.ARLBNF_SS2, LOGPVALUE > 1.30103 & `XPEHH_ARLBNF_ScotianShelf` > 0 & cand_region == "cand")

xpehh_ARLBNFSS_neg <- filter(sig_xpehh.ARLBNF_SS2, LOGPVALUE > 1.30103 & `XPEHH_ARLBNF_ScotianShelf` < 0 & cand_region == "cand")

xpehh_ARLBNFSS_plot <- ggplot()+
  geom_point(data=xpehh.ARLBNF_SS.chr, aes(x=plotorder,y=XPEHH_ARLBNF_ScotianShelf,color=target_chr), alpha=0.5, size=0.7)+
  scale_color_manual(values=chr_colors)+
  geom_point(data=xpehh_ARLBNFSS_pos, aes(x=plotorder, y=XPEHH_ARLBNF_ScotianShelf), size=1.5,color="#ABD9E9")+
  geom_point(data=xpehh_ARLBNFSS_neg, aes(x=plotorder, y=XPEHH_ARLBNF_ScotianShelf), size=1.5,colour="#A50026") +
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
        #panel.grid.minor.y=element_blank(),
        #panel.grid.major.y=element_blank(),
        text=element_text(family="serif", size=12),axis.text.x=element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  ylab("XP-EHH")+
  ylim(c(-9.5,9.5))

xpehh_ARLBNFSS_plot
ggsave("rehh_NBW2_2M_ggplot_ARLBNFSS_log1.3_test6.png", width=10,height=2.5,dpi=1000)
                         
                                              
