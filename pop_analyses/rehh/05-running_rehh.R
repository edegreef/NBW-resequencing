# Using rehh for haplotype analyses to identify regions under selection

library(rehh)
library(vcfR)
library(tidyverse)
library(ggplot2)
library(DescTools)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

# Create map file for whole genome
vcf <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.imputed.vcf.gz")
chroms <- unique(vcf@fix[,1])
colnames(vcf@fix)
map_file <- tibble(marker_ID = vcf@fix[,3], chromosome = vcf@fix[,1], position = vcf@fix[,2], ref = vcf@fix[,4], alt = vcf@fix[,5])
map_file$marker_ID <- paste(map_file$chromosome, map_file$position, sep="-")

# Load loop output scan file, chr1-21 and chru (unlocalized)
# Make list of filenames to import
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

# Swtich back to previous directory
setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

# Calculate the XP-EHH statistic between each pair of populations (Northern vs Scotian Shelf)
xpehh.NOR_SS <- ies2xpehh(scan_pop1 = NOR_scan_merged,
                          scan_pop2 = SS_scan_merged,
                          popname1 = "Northern",
                          popname2 = "ScotianShelf",
                          p.adjust.method = "fdr")

# Add snp ID
xpehh.NOR_SS$SNP <- paste(xpehh.NOR_SS$CHR, xpehh.NOR_SS$POSITION, sep="-")

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.05
#q < 0.05  (-log10(q) > 1.30103)
#if using q < 0.01, then -log10(q) > 2 

NORSS_candidates <- calc_candidate_regions(xpehh.NOR_SS, 
                                           threshold = 1.30103, 
                                           pval = TRUE, 
                                           window_size = 1E5, 
                                           overlap = 1E4, 
                                           min_n_extr_mrk = 2)
NORSS_candidates

# Add chromosome info to xpehh.NOR_SS  (for plotting)
scafs  <- read.csv("C:/Users/eveli/Dropbox/NBW-me/reference/satsuma_outputs/nbw_scaf_chr_info_autosomes.csv")
scafs$contig <- paste("contig", scafs$query_scaf, sep="")

# Add chr info to xpehh
xpehh.NOR_SS.chr <- xpehh.NOR_SS %>% left_join(scafs, by=c("CHR"="contig"))
xpehh.NOR_SS.chr$target_chr[is.na(xpehh.NOR_SS.chr$target_chr)] <- "chru"

# Adding new position number for within chromosome
xpehh.NOR_SS.chr <- xpehh.NOR_SS.chr %>% group_by(target_chr) %>% mutate(new_pos = row_number())

# Order chr1, chr2... instead of chr1, chr10..
good_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chru")

# Proper chr order
xpehh.NOR_SS.chr <- xpehh.NOR_SS.chr %>% arrange(factor(target_chr, levels=good_order))
xpehh.NOR_SS.chr$plotorder <- 1:nrow(xpehh.NOR_SS.chr)

# remove chru?
#xpehh.NOR_SS.chr <- subset(xpehh.NOR_SS.chr, target_chr != "chru")

# Filter for significant hits among the NOR_SS XP-EHH SNPs
sig_xpehh.NOR_SS <- filter(xpehh.NOR_SS.chr, LOGPVALUE > 1.30103)

# Check which SNPs appear within a candidate region
sig_xpehh.NOR_SS$cand_region <- ifelse(sapply(sig_xpehh.NOR_SS$POSITION, function(p) 
  any(NORSS_candidates$START <= p & NORSS_candidates$END >= p)),"cand", NA)

# Remove single marks
sig_xpehh.NOR_SS2<- sig_xpehh.NOR_SS[sig_xpehh.NOR_SS$CHR %in% sig_xpehh.NOR_SS$CHR[duplicated(sig_xpehh.NOR_SS$CHR)],]

sig_xpehh.NOR_SS$cand_region <- ifelse(sapply(seq_along(sig_xpehh.NOR_SS$POSITION), function(i) {
  inds <- NORSS_candidates$START <= sig_xpehh.NOR_SS$POSITION[i] & NORSS_candidates$END >= sig_xpehh.NOR_SS$POSITION[i]
  any(inds) & (sig_xpehh.NOR_SS$CHR[i] == NORSS_candidates$CHR[which.max(inds)])
}), "cand", NA)

# Pull positive (selection in Northern) and negative (selection in Scotian Shelf) SNPs for highlighting on the manhattan plot
xpehh_NORSS_pos <- filter(sig_xpehh.NOR_SS2, LOGPVALUE > 1.30103 & `XPEHH_Northern_ScotianShelf` > 0 & cand_region == "cand")

xpehh_NORSS_neg <- filter(sig_xpehh.NOR_SS2, LOGPVALUE > 1.30103 & `XPEHH_Northern_ScotianShelf` < 0 & cand_region == "cand")


# Prepare plot colours
chr_colors <- c(chr1="gray60", chr2="gray70", chr3="gray60", chr4="gray70", chr5="gray60", chr6="gray70", chr7="gray60", chr8="gray70", chr9="gray60", chr10="gray70", chr11="gray60", chr12="gray70", chr13="gray60", chr14="gray70", chr15="gray60", chr16="gray70", chr17="gray60", chr18="gray70", chr19="gray60", chr20="gray70", chr21="gray60", chru="gray70")

# Work on split later
# Regular gray dot size 1.5, and bold dot size 2.5

# Intercept +- 4.56 for log 2, +- 4.1302 for log 1.30103
xpehh_NORSS_plot <- ggplot()+
  geom_point(data=xpehh.NOR_SS.chr, aes(x=plotorder,y=XPEHH_Northern_ScotianShelf,color=target_chr), alpha=0.6, size=0.7)+
  scale_color_manual(values=chr_colors)+
  geom_hline(yintercept=4.1302, linetype="dashed", color="black")+
  geom_hline(yintercept=-4.1302, linetype="dashed", color="black")+
 geom_point(data=xpehh_NORSS_pos, aes(x=plotorder, y=XPEHH_Northern_ScotianShelf), size=1.5,color="#4575B4")+
geom_point(data=xpehh_NORSS_neg, aes(x=plotorder, y=XPEHH_Northern_ScotianShelf), size=1.5,colour="#D73027") +
 # geom_point(data=SV_pos, aes(x=bin_order, y=XPEHH_Northern_ScotianShelf.x), size=1.5,colour="black") +
#  geom_point(data=SV_neg, aes(x=bin_order, y=XPEHH_Northern_ScotianShelf.x), size=1.5,colour="black") +
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(),text=element_text(family="serif", size=12),axis.text.x=element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  ylab("XP-EHH")

xpehh_NORSS_plot
ggsave("rehh_NBW2_2M_ggplot_NORSS_log1.3.png", width=13,height=5,dpi=1000)

## add snp names ?
xpehh_NORSS_plot + geom_text(data=sig_xpehh.NOR_SS,aes(x=plotorder, y=XPEHH_Northern_ScotianShelf,label=SNP,alpha=0.5), hjust=0,vjust=1)
#geom_text(aes(label=sample, alpha=0.5), hjust=0,vjust=1)
#ggsave("test21ggplot_redo_check_SNPlabel.png", width=40,height=20,dpi=600)


## histogram
xpehh_hist <-ggplot(xpehh.NOR_SS.chr, aes(x=XPEHH_Northern_ScotianShelf)) + 
  geom_histogram(color="black", fill="steelblue", bins=50)+
  theme_classic()
xpehh_hist

# Save lists for significant snps
write.csv(xpehh_NORSS_pos, "sigsnps_xpehh_northern_logp1.3.csv")
write.csv(xpehh_NORSS_neg, "sigsnps_xpehh_scotianshelf_logp1.3.csv")
