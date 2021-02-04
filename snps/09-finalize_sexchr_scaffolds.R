# 1) look at annotated bed file to create final list of sex-linked scaffolds, 2) make upset plots to show overlap of scaffolds for the combos, and 3) check final X and Y with satsuma alignment 

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/MFcoverage/series")

# 1) create final list of sex-linked scaffolds (doing X and Y separately)
library(tidyverse)

bed="NBW_genome.10KBwindows.Xlinked.bed"
type="Xlinked"
scaffold_info="Northern_bottlenose_whale_051018_shortLabel.fasta.fai"

# bedtools command for ref to remember order of columns
#bedtools annotate -i NBW_genome.10KBwindows.bed -files X_linked_scafwindows.M1F1.bed X_linked_scafwindows.M1F2.bed X_linked_scafwindows.M2F1.bed X_linked_scafwindows.M2F2.bed Control_scafwindows_remove1scaf.M1M2.bedControl_scafwindows.F1F2.bed > NBW_genome.10KBwindows.Xlinked.bed

# load data and rename scaffolds
sexlinked_windows <- read_tsv(file = bed ,col_names = F) %>% rename(scaf=X1,start=X2,stop=X3,M1F1=X4,M1F2=X5, M2F1=X6, M2F2=X7, M1M2=X8, F1F2=X9)

# want to keep scaffolds that are enriched in the 4 MF combos and not enriched in controls (looking for non-zero values across)
# for the Y chr, do not include the FF control
filtered_windows <- sexlinked_windows %>% filter(M1F1 !=0 & M1F2 !=0 & M2F1 !=0 & M2F2 !=0 & M1M2 !=0)# & F1F2 !=0)

# add scaffold lengths
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","scaffold_length"))
filtered_windows <- full_join(filtered_windows,scaffold_lengths) %>% na.omit()

# keep only unique scaf names for final list
sex_linked_scaffolds <- filtered_windows[!duplicated(filtered_windows$scaf), ]
sex_linked_scaffolds_list <- sex_linked_scaffolds[,1]

# check total scaffold lengths to see that it is reasonable (should be close to the checks when looking at difcover results)
chrlength <- sum(sex_linked_scaffolds$scaffold_length)
genome_length <- sum(scaffold_lengths$scaffold_length)
genome_length
chrlength/genome_length

# save list of scaffolds to use as reference for the sex chrom in the genome!
write.table(sex_linked_scaffolds_list,paste("final_scaffolds_", type, ".txt",sep=""), col.names = FALSE, row.names = FALSE)


# 2) make upset plot to see how many sex-linked scaffolds were identified in combos of M/F runs

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(RColorBrewer)

# load X or Y bed
sexlinked_windows <- read_tsv(file = bed ,col_names = F) %>% rename(scaf=X1,start=X2,stop=X3,M1F1=X4,M1F2=X5, M2F1=X6, M2F2=X7, M1M2=X8, F1F2=X9)

# make dataframe into matrix with 0's and 1's for each combo
windows <- sexlinked_windows[,1:3]
sexcombo <- sexlinked_windows[,4:9]
convert01 <- sexcombo %>% mutate_if(is.numeric, ~1 * (. != 0))
table <- cbind(windows, convert01)

# make dataframe to include by scaffold only (removing window info so no "duplicate" scaffold names). Do for each run then merge later
run1 <- table %>% filter(M1F1 !=0)
unique_run1 <- run1[!duplicated(run1$scaf), ]
unique_run1 <- select(unique_run1, scaf, M1F1)

run2 <- table %>% filter(M1F2 !=0)
unique_run2 <- run2[!duplicated(run2$scaf), ]
unique_run2 <- select(unique_run2, scaf, M1F2)

run3 <- table %>% filter(M2F1 !=0)
unique_run3 <- run3[!duplicated(run3$scaf), ]
unique_run3 <- select(unique_run3, scaf, M2F1)

run4 <- table %>% filter(M2F2 !=0)
unique_run4 <- run4[!duplicated(run4$scaf), ]
unique_run4 <- select(unique_run4, scaf, M2F2)

run5 <- table %>% filter(M1M2 !=0)
unique_run5 <- run5[!duplicated(run5$scaf), ]
unique_run5 <- select(unique_run5, scaf, M1M2)

run6 <- table %>% filter(F1F2 !=0)
unique_run6 <- run6[!duplicated(run6$scaf), ]
unique_run6 <- select(unique_run6, scaf, F1F2)

scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","scaffold_length"))
scaf_list <- scaffold_lengths[,1]

# merge columns for each run
merge_scaf_list <- scaf_list %>% full_join(unique_run1, by=c("scaf"="scaf")) %>% full_join(unique_run2, by=c("scaf"="scaf")) %>% full_join(unique_run3, by=c("scaf"="scaf")) %>% full_join(unique_run4, by=c("scaf"="scaf")) %>% full_join(unique_run5, by=c("scaf"="scaf"))# %>% full_join(unique_run6, by=c("scaf"="scaf"))

# convert NAs to 0's
merge_scaf_list[is.na(merge_scaf_list)] <- 0
df <- as.data.frame(merge_scaf_list)

# check that the values are correct
sum(df$M1F1)

# convert to makecombmat format
df <- df[,-1]
df <- make_comb_mat(df, mode="intersect")

# see which combos want to extract
df[1:4]

# all 6 runs: code 111111
#df_extracted <- extract_comb(df, "111111", "111110")
#UpSet(df_extracted)

# x-linked
brewer.pal(n=6, name="Blues")
UpSet((df[comb_size(df) <= 10000]), 
      set_order=c("M1F1", "M1F2", "M2F1", "M2F2", "M1M2", "F1F2"),
      lwd=4, pt_size=unit(5, "mm"),
      comb_col=c("#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "black")[comb_degree((df[comb_size(df) <= 10000]))],
      border=TRUE,
      top_annotation=upset_top_annotation((df[comb_size(df) <= 10000]), annotation_name_rot=90, annotation_name_side="left", axis_param=list(side="left")),
      right_annotation=upset_right_annotation((df[comb_size(df) <= 10000]), gp=gpar(fill="gray", border=FALSE), annotation_name_side="top", axis_param=list(side="top")),
      column_title="X-linked scaffolds")


# y-linked
brewer.pal(n=5, name="Blues")
UpSet((df[comb_size(df) <= 10000]), 
      set_order=c("M1F1", "M1F2", "M2F1", "M2F2", "M1M2"),
      lwd=4, pt_size=unit(5, "mm"),
      comb_col=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "black")[comb_degree((df[comb_size(df) <= 10000]))],
      border=TRUE,
      top_annotation=upset_top_annotation((df[comb_size(df) <= 10000]), annotation_name_rot=90, annotation_name_side="left", axis_param=list(side="left")),
      right_annotation=upset_right_annotation((df[comb_size(df) <= 10000]), gp=gpar(fill="gray", border=FALSE), annotation_name_side="top", axis_param=list(side="top")),
      column_title="Y-linked scaffolds")


# 3) see how well satsuma synteny with blue whale matches

# load satsuma_summary.chained.out file and rename columns
satsuma <- "C:/Users/Evelien de Greef/Dropbox/NBW-me/reference/satsuma_outputs/satsuma_summary.chained.chrX.out"
satsuma <- read.delim(satsuma, header=F) %>% rename(target_chr=V1,start_target_base=V2,end_target_base=V3,query_scaf=V4,start_query_base=V5, end_query_base=V6, identity=V7, orientation=V8)

satsuma$bases_spanned <- satsuma$end_query-satsuma$start_query

# add scaffold length info
satsuma2 <- satsuma %>% full_join(scaffold_lengths, by=c("query_scaf"="scaf")) %>% na.omit()
hist(satsuma2$identity)

# look at proportion of scaf with identity > 0.6
satsuma3 <- satsuma2 %>% filter(identity > 0.6) 

list_satsum <- satsuma3[!duplicated(satsuma3$query_scaf), ]
satsum_length <- sum(list_satsum$scaffold_length)
satsum_length/genome_length

list_satsum$satsuma <- "YES"
merge_check <- sex_linked_scaffolds %>% full_join(list_satsum, by=c("scaf"="query_scaf")) %>% na.omit()

total_crosscheck <- sum(merge_check$scaffold_length.y)
total_crosscheck/genome_length

percent_match <- (total_crosscheck/genome_length)/(satsum_length/genome_length)
percent_match
