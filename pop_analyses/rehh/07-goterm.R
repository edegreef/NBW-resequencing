# Prepping gene list and then using enrichR for go-term analysis
# Help from Matt Thorstensen for initial set up

#library(devtools)
#install_github("wjawaid/enrichR")
library(enrichR)
library(tidyverse)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/rehh")

# Read in magma outputs 
# Northern
xpehh_NOR_pos <- read_delim("NOR_log1.3_magma_annotation.genes.annot", delim = "\t", skip = 3, col_names = c("geneID", "gene_loc", "snp", "snp2", "snp3", "snp4", "snp5", "snp6", "snp7", "snp8", "snp9", "snp10", "snp11", "snp12", "snp13", "snp14"))

# Scotian Shelf
xpehh_SS_pos <- read_delim("SS_log1.3_magma_annotation.genes.annot", delim = "\t", skip = 3, col_names = c("geneID", "gene_loc", "snp", "snp2", "snp3", "snp4", "snp5", "snp6", "snp7", "snp8", "snp9", "snp10", "snp11", "snp12", "snp13", "snp14"))

# Read in gff file with just genes
genes_raw <- read_delim("C:/Users/eveli/Dropbox/NBW-me/reference/annotation/NBW.genes.full.gff", delim = "\t", col_names=c("scaffold","prog","type","bp_start","bp_end","X", "X1", "X2", "attributes"))

# Re-format gene list
# Extract GMOD id's
split=str_split_fixed(genes_raw$attributes, ";", 2)
split[,1]=gsub("ID=","",split[,1])

# Extract actual gene names & descriptions
split2=str_split_fixed(genes_raw$attributes,"Note=",2)
split2[,2]=gsub("Similar to ","",split2[,2])
split3=str_split_fixed(split2[,2],":",2)

# Merge to new file
genes <- as.data.frame(cbind(split[,1], split3))
colnames(genes) <- c("GMOD", "gene", "description")

rm(split, split2, split3)

xpehh_NOR_genes <- xpehh_NOR_pos %>% left_join(genes, by=c("geneID"="GMOD"))
xpehh_SS_genes <- xpehh_SS_pos %>% left_join(genes, by=c("geneID"="GMOD"))

xpehh_NOR_genes$group <- "Northern"
xpehh_SS_genes$group <- "Scotian"

# Re-order
col_order <- c("group","geneID", "gene_loc", "gene", "description","snp", "snp2")#, "snp3", "snp4", "snp5", "snp6")
xpehh_NOR_genes2 <- xpehh_NOR_genes[, col_order]
xpehh_SS_genes2 <- xpehh_SS_genes[, col_order]

# Extract first 5 columns
xpehh_NOR_genes2 <- xpehh_NOR_genes2[,c(1:5)]
xpehh_SS_genes2 <- xpehh_SS_genes2[,c(1:5)]

merged <- rbind(xpehh_NOR_genes2, xpehh_SS_genes2)

#write.csv(merged, "merged_list_candidategenes_full_20kb_withdups.csv")

# Extract contig value
merged$contig <- sub("\\:.*", "", merged$gene_loc)

# Add chr info
scafs  <- read.csv("C:/Users/eveli/Dropbox/NBW-me/reference/satsuma_outputs/nbw_scaf_chr_info_autosomes.csv")
scafs$query_scaf <- as.character(scafs$query_scaf)

# Add chr info to xpehh
merged.chr <- merged %>% left_join(scafs, by=c("contig"="query_scaf"))
merged.chr$target_chr[is.na(merged.chr$target_chr)] <- "chru"

# Re-order
col_order <- c("group","geneID","target_chr", "contig", "scaffold_length","start_target_base","orientation","gene_loc", "gene", "description")
#"snp", "snp2", "snp3", "snp4", "snp5", "snp6", "snp7", "snp8", "snp9", "snp10", "snp11", "snp12", "snp13", "snp14")

merged.chr <- merged.chr[, col_order]
#merged.chr <- merged.chr[order(-merged.chr$snp_count),]

#write.csv(merged.chr, "merged_list_candidategenes_full_20kb_aed0.25_withdups_info.csv")

# Remove "duplicates"
xpehh_NOR_genes3 <- xpehh_NOR_genes2[!duplicated(xpehh_NOR_genes2$gene), ]
xpehh_SS_genes3 <- xpehh_SS_genes2[!duplicated(xpehh_SS_genes2$gene), ]

# Merge and save
xpehh_NOR_genes3$group <- "Northern"
xpehh_SS_genes3$group <- "Scotian"

merge <- rbind(xpehh_NOR_genes3, xpehh_SS_genes3)
#write.csv(merge, "merged_list_20kb_aed0.25_candidategenes.csv")

# Take distinct/unique gene Ids
NOR_genes <- distinct(xpehh_NOR_genes, gene)
NOR_genes <- as.vector(NOR_genes)
NOR_genes <- subset(NOR_genes, gene != "PSMD12" & gene != "PMM1")

SS_genes <- distinct(xpehh_SS_genes, gene)
SS_genes <- as.vector(SS_genes)
SS_genes <- subset(SS_genes, gene != "PSMD12" & gene != "PMM1")

write.csv(NOR_genes, "NOR_genes_list.csv")
write.csv(SS_genes, "SS_genes_list.csv")

# Taking a look at what databases are available in EnrichR
listEnrichrDbs()

# Run enrichR, searching the biological process and molecular function databases
NOR_enriched <- enrichr(NOR_genes$gene, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
NOR_enriched_BP <- as_tibble(NOR_enriched[["GO_Biological_Process_2018"]])
NOR_enriched_MF <- as_tibble(NOR_enriched[["GO_Molecular_Function_2018"]])
NOR_enriched_CC <- as_tibble(NOR_enriched[["GO_Cellular_Component_2018"]])

# Run enrichR, searching the biological process and molecular function databases
SS_enriched <- enrichr(SS_genes$gene, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
SS_enriched_BP <- as_tibble(SS_enriched[["GO_Biological_Process_2018"]])
SS_enriched_MF <- as_tibble(SS_enriched[["GO_Molecular_Function_2018"]])
SS_enriched_CC <- as_tibble(SS_enriched[["GO_Cellular_Component_2018"]])

# Selecting results with adjusted.p.value < 0.05

NOR_enriched_BP_filter <- dplyr::filter(NOR_enriched_BP, Adjusted.P.value < 0.05)
NOR_enriched_BP_filter$type <- "NOR_BP"
NOR_enriched_MF_filter <- dplyr::filter(NOR_enriched_MF, Adjusted.P.value < 0.05)
NOR_enriched_MF_filter$type <- "NOR_MF"
NOR_enriched_CC_filter <- dplyr::filter(NOR_enriched_CC, Adjusted.P.value < 0.05)
NOR_enriched_CC_filter$type <- "NOR_CC"

SS_enriched_BP_filter <- dplyr::filter(SS_enriched_BP, Adjusted.P.value < 0.05)
SS_enriched_BP_filter$type <- "SS_BP"
SS_enriched_MF_filter <- dplyr::filter(SS_enriched_MF, Adjusted.P.value < 0.05)
SS_enriched_MF_filter$type <- "SS_MF"
SS_enriched_CC_filter <- dplyr::filter(SS_enriched_CC, Adjusted.P.value < 0.05)
SS_enriched_CC_filter$type <- "SS_CC"

all <- rbind(NOR_enriched_BP_filter,NOR_enriched_MF_filter,NOR_enriched_CC_filter,SS_enriched_BP_filter,SS_enriched_MF_filter,SS_enriched_CC_filter)

write.csv(all, "GO_results_20kb.csv")
