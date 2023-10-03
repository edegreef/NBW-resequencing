#!/usr/bin/Rscript

# Run data2haplohh in rehh to create the wgscan file
# Rehh vignette: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
# Help from Matt Thorstensen to set this up

library(rehh)
library(vcfR)
library(tidyverse)

# Create map file, do this for each chromosome
vcf <- read.vcfR("NBW2_SNPS_2M.abbr.37.min50kb.imputed.chr1.vcf")
chroms <- unique(vcf@fix[,1])
colnames(vcf@fix)

map_file <- tibble(marker_ID = vcf@fix[,3], chromosome = vcf@fix[,1], position = vcf@fix[,2], ref = vcf@fix[,4], alt = vcf@fix[,5])
write_delim(map_file, "sep_chrs/rehh_mapfile_chr1.inp", delim = "\t", col_names = FALSE)

map_file$marker_ID <- paste(map_file$chromosome, map_file$position, sep="-")

# I ran this loop step for each chr separately and for each group, and did this on biology-01 cluster cuz it took a while, then reloaded the scan output file later for rehh analyses.

#R --vanilla < prep_rehh_files.R

# -------------IC
for(i in 1:length(chroms)) {
 #create internal representation
  hh <- data2haplohh(hap_file = "NBW2_SNPS_2M.abbr.37.min50kb.imputed.chr1.IC.vcf",
                     map_file = "rehh_mapfile_chr1.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
# perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
# concatenate chromosome-wise data frames to a data frame for the whole genome (more efficient ways certainly exist...)
  if (i == 1) {
    IC_chr1_wgscan <- scan
  } else {
    IC_chr1_wgscan <- rbind(IC_chr1_wgscan, scan)
  }
}

write.csv(IC_chr1_wgscan, "IC_chr1_wgscan.csv")


# ------------NOR
for(i in 1:length(chroms)) {
 #create internal representation
  hh <- data2haplohh(hap_file = "NBW2_SNPS_2M.abbr.37.min50kb.imputed.chr1.NOR.vcf",
                     map_file = "rehh_mapfile_chr1.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
# perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
# concatenate chromosome-wise data frames to a data frame for the whole genome (more efficient ways certainly exist...)
  if (i == 1) {
    NOR_chr1_wgscan <- scan
  } else {
    NOR_chr1_wgscan <- rbind(NOR_chr1_wgscan, scan)
  }
}

write.csv(NOR_chr1_wgscan, "NOR_chr1_wgscan.csv")


# ------------SS
for(i in 1:length(chroms)) {
 #create internal representation
  hh <- data2haplohh(hap_file = "NBW2_SNPS_2M.abbr.37.min50kb.imputed.chr1.SS.vcf",
                     map_file = "rehh_mapfile_chr1.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
# perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
# concatenate chromosome-wise data frames to a data frame for the whole genome (more efficient ways certainly exist...)
  if (i == 1) {
    SS_chr1_wgscan <- scan
  } else {
    SS_chr1_wgscan <- rbind(SS_chr1_wgscan, scan)
  }
}

write.csv(SS_chr1_wgscan, "SS_chr1_wgscan.csv")

