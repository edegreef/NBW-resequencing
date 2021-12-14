# Looking at breakdancer output file to see if which/any regions in the xp-ehh analyses are inversions

library(tidyverse)

# From breakdancer's readme file (https://github.com/kenchen/breakdancer#readme):
#BreakDancer's output file consists of the following columns:

#1. Chromosome 1
#2. Position 1
#3. Orientation 1
#4. Chromosome 2
#5. Position 2
#6. Orientation 2
#7. Type of a SV
#8. Size of a SV
#9. Confidence Score
#10. Total number of supporting read pairs
#11. Total number of supporting read pairs from each map file
#12. Estimated allele frequency
#13. Software version
#14. The run parameters

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/breakdancer")

# Load breakdancer output
ctx_combined10 <- read.delim("combine_10_samples.redoheader.deDupRG.copy.ctx", header=T)

# Filter for minimum score of 50
ctx_combined10_2 <- subset(ctx_combined10, Score > 50)

# Pull out INV (inversions), CTX (inter-chromosomal translocation), and ITX (intra-chromosomal translocation)
ctx_combined10_2 <- subset(ctx_combined10_2, Type=="INV" | Type=="CTX" | Type=="ITX")

# Need to add in scaffold info to get right start & stops
scafs  <- read.csv("C:/Users/eveli/Dropbox/NBW-me/reference/scaffoldlengths.csv")

# If chr matches then add length. Note that Chr1 and Chr2 is referring to first and second scaffold, not chromosome#1 and chromosome#2.
ctx_combined10_3 <- ctx_combined10_2 %>% left_join(scafs, by=c("X.Chr1"="scaffold"))
colnames(ctx_combined10_3)[12] <- "bp_length_chr1"

ctx_combined10_3 <- ctx_combined10_3 %>% left_join(scafs, by=c("Chr2"="scaffold"))
colnames(ctx_combined10_3)[13] <- "bp_length_chr2"

# Make list from inversions that are on one/same scaffold
same_scaf <- subset(ctx_combined10_3, X.Chr1==Chr2)
samescaf_pos <- same_scaf[,c("X.Chr1", "Pos1", "Pos2")]

# Make list from inversions that ate on different scaffolds
diff_scaf <- subset(ctx_combined10_3, X.Chr1!=Chr2)

# This one maybe have to use chr1 pos up till bp_length. and reverse for chr2
diff_scaf$chr2start <- 1

# Make list of scaffolds and position interval for filtering
# Diff scaf
diffscaf_pos <- diff_scaf[,c("X.Chr1", "Pos1", "bp_length_chr1", "Chr2", "chr2start", "Pos2")]
diffscaf_pos <- subset(diffscaf_pos, Pos2 > 1)

# Load snps to make a CHROM POS list for filtering
vcf <- read.table("C:/Users/Evelien de Greef/Dropbox/NBW-me/demography/smcpp2/NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss04.n36.min100kb.vcf.gz")

snplist <- vcf[,c("V1", "V2")]
colnames(snplist)[1] <- "CHR"
colnames(snplist)[2] <- "POSITION"

snplist_same <- snplist
snplist_diffchr1 <- snplist
snplist_diffchr2 <- snplist

#V1=scaf
#V2=pos

# For list of same scaf svs
samescaf_pos$CHR1 <- paste("contig",samescaf_pos$X.Chr1,sep="")

snplist_same$SV <- ifelse(sapply(seq_along(snplist$POSITION), function(i) {
  inds <- samescaf_pos$Pos1 <= snplist_same$POSITION[i] & samescaf_pos$Pos2 >= snplist_same$POSITION[i]
  any(inds) & (snplist_same$CHR[i] == samescaf_pos$CHR1[which.max(inds)])
}), "SV_Schr1", "not_sv")


# For list of diff scaf svs on chr1
diffscaf_pos$CHR1 <- paste("contig",diffscaf_pos$X.Chr1,sep="")

snplist_diffchr1$SV <- ifelse(sapply(seq_along(snplist_diffchr1$POSITION), function(i) {
  inds <- diffscaf_pos$Pos1 <= snplist_diffchr1$POSITION[i] & diffscaf_pos$bp_length_chr1 >= snplist_diffchr1$POSITION[i]
  any(inds) & (snplist_diffchr1$CHR[i] == diffscaf_pos$CHR1[which.max(inds)])
}), "SV_Dchr1", "not_sv")

# For list of diff scaf svs on chr2
diffscaf_pos$CHR2 <- paste("contig",diffscaf_pos$Chr2,sep="")

snplist_diffchr2$SV <- ifelse(sapply(seq_along(snplist_diffchr2$POSITION), function(i) {
  inds <- diffscaf_pos$chr2start <= snplist_diffchr2$POSITION[i] & diffscaf_pos$Pos2 >= snplist_diffchr2$POSITION[i]
  any(inds) & (snplist_diffchr2$CHR[i] == diffscaf_pos$CHR1[which.max(inds)])
}), "SV_Dchr2", "not_sv")

# Prep list for all of them together
same_SV <- subset(snplist_same, SV=="SV_Schr1")
diffchr1_SV <- subset(snplist_diffchr1, SV=="SV_Dchr1")
diffchr2_SV <- subset(snplist_diffchr2, SV=="SV_Dchr2")

SV_list <- rbind(same_SV, diffchr1_SV, diffchr2_SV)
SV_list$chr_pos <- paste(SV_list$CHR, SV_list$POSITION, sep=" ")

# Save list
#write.csv(SV_list, "SV_list_NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss01.n36.imputed.mac2.min10kbCHROM_POS.csv")

# Also remove contig208, 230415, 230424, 230611, 231412, 231787, 232034, 232430, 233046 - these were identified as inversion as well

#x <- read.csv("SV_list_NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss04.n36.min100kb_CHROM_POS.csv", header=T)

# Pull snps from these contigs that have inversions as well
pull <- subset(snplist, CHR=="contig208" | CHR=="contig230415" | CHR=="contig230424" | CHR=="contig230611" | CHR=="contig231412" | CHR=="contig231787" | CHR=="contig232034" | CHR=="contig232430" | CHR=="contig233046")

pull$SV <- "part2"
pull$chr_pos <- paste(pull$CHR, pull$POSITION, sep=" ")

#x <- x[,-1]
#colnames(x)[4] <- "chr_pos"

new_x <- rbind(x, pull)
write.csv(new_x, "SV_list_NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss04.n36.min100kb_withPart2_CHROM_POS.csv")
