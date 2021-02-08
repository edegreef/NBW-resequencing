# 1) Looking at FST between populations. Help from Matt Thorstensen to for code to run Hierfstat.
# 2) Making fst heatmap
# 3) Examining isolation-by-distance pattern & running mantel test


# 1) estimate FST 
library(hierfstat)
library(vcfR)
library(tidyverse)
library(adegenet)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/fst")

# Read in vcf
vcf <- read.vcfR("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.vcf") 

# Reformat vcfR object into a genlight object
vcf <- vcfR2genlight(vcf)

# Reformat to a matrix for input into hierfstat
mat <- as.matrix(vcf)

# Prepare 0, 1, 2 matrix in hierfstat format
snps <- mat
snps[snps==0] <- 11
snps[snps==1] <- 12
snps[snps==2] <- 22

# using partial data as test run
#subset <- snps[,1:1000]

df <- rownames_to_column(as.data.frame(snps), var="Ind")

# bring in sample info
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info.csv", header=T)
sample_info_subset <- subset(sample_info, remove_indiv!="Y") #removing individuals with high miss and kin pair

# extract pop info
pops <- tibble(ID=sample_info_subset$sample, region=sample_info_subset$region)
colnames(pops)[1] <- "Ind"

# make pops info numeric (automatically goes in alphabetical order)
# 1=Arctic, 2=Iceland, 3=Labrador, 4=Newfoundland, 5=Scotian_shelf
pops$region<- as.numeric(as.factor(pops$region))

# Left_join info to snps
df <- left_join(pops, df) %>%
  column_to_rownames(var="Ind")

save(df, file="NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.df.leftjoined.RData")

# Look at pairwise Fst with hierfstat
wc <- pairwise.WCfst(df, diploid = TRUE)
wc
write.csv(wc,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.csv")

# Estimate confidence intervals
wc_bootstrp_confidence <- boot.ppfst(dat = df, nboot = 1000)
wc_bootstrp_confidence$ll
wc_bootstrp_confidence$ul

# save outputs
write.csv(wc_bootstrp_confidence$ll,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ll.csv")
write.csv(wc_bootstrp_confidence$ul,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ul.csv")
save(wc_bootstrp_confidence,file="NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.wc_bootstrp_CI.RData")


# 2) fst heatmap
library(reshape2)
library(ggplot2)

# load in fst matrix - used output from hierfstat, then re-ordered pops by latitude and adjusted fst to slatkin's linearized fst (fst-(1-fst))
fst <- read.csv("autosomes.LDpruned.n36.FST.slat.matrix.csv", header=T)

# edit column and row names
colnames(fst) <- c("X", "IC", "AR", "LB", "NF", "SS")
fst$X <- c("IC", "AR", "LB", "NF", "SS")

# convert matrix to melted dataframe to use in geom_tile
melt_fst <- melt(fst)

# define order for plot axis
x_pop_order <- c("SS", "NF", "LB", "AR", "IC")
y_pop_order <-  c("SS", "NF", "LB", "AR", "IC")

fst <- ggplot(melt_fst, aes(x=X, y=variable, fill=value)) + 
  geom_tile()+
  scale_x_discrete(position="top", limits=x_pop_order)+
  xlab("")+
  ylab("")+
  # xlim(level=x_pop_order)+
  ylim(level=y_pop_order)+
  scale_fill_distiller(palette="Blues", trans="reverse", na.value = "white")+
  theme_minimal()+
  ggtitle("FST between regions")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(hjust = 0))+
  geom_text(aes(label = round(value, 4)))+
  labs(fill="FST")

fst
ggsave("FST_autosomes_LDpruned_n36.png", width=6, height=5, dpi=300)


# 3) Isolation-by-distance
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(RColorBrewer)

# load in data (in data, each pop pair is a separate row, with FST values and distance as columns)
data <- read.csv("site_mean_coords_fordist_newsnps_updated.csv", header=T)
data$pair_abbrev

# make a dataset excluding labrador
data_noLB <- data[-c(1,3,7,9),]
data_noLB$pair_abbrev

# looking at color palettes
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(9, "Blues")
brewer.pal(9, "Blues")

# manually setting color palette (using 'Blues' + one black)
pal <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "black")

# plot IBD
IBD <- ggplot(data, aes(x=shelf_distance_km,y=slatkins_fst))+
  geom_smooth(method = "lm", se=TRUE, fill="gray85", formula=y~x, linetype=0, color="gray")+
  geom_pointrange(aes(ymin=slatkins_ll, ymax=slatkins_ul), size=1, color="black")+
  geom_point(aes(fill=pair_abbrev), colour="black",pch=21,size=4) +
  theme_classic()+
  #ggtitle("Isolation by distance")+
  #annotate("text", x=5200, y=0.002, label= "Mantel R2 = 0.54", size=4, fontface="italic", color="black")+
  scale_fill_manual(values=rev(pal), breaks=c("IC-LB", "IC-SS", "LB-SS", "IC-NF", "NF-SS", "AR-SS", "AR-LB", "AR-NF", "LB-NF", "IC-AR"))+
  xlab("Distance (km)")+
  ylab("FST")+
  xlim(900,5700)+
  geom_text(aes(label=pair_abbrev), color="black",size=3, hjust=-0.3,vjust=0.5)+
  labs(fill="Region pair")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
IBD
ggsave("IBD_plot_blues_pointlabels.png", width=7, height=6, dpi=300)

# plot IBD without LB site
IBD_noLB <- ggplot(data_noLB, aes(x=shelf_distance_km,y=slatkins_fst))+
  geom_smooth(method = "lm", se=TRUE, fill="gray85", formula=y~x, color="gray", linetype=0)+
  geom_pointrange(aes(ymin=slatkins_ll, ymax=slatkins_ul), size=1, color="black")+
  geom_point(aes(fill=pair_abbrev), colour="black",pch=21,size=4) +
  theme_classic()+
  #ggtitle("Isolation by distance")+
  #annotate("text", x=5200, y=0.002, label= "Mantel R2 = 0.54", size=4, fontface="italic", color="black")+
  scale_fill_manual(values=rev(pal), breaks=c("IC-LB", "IC-SS", "LB-SS", "IC-NF", "NF-SS", "AR-SS", "AR-LB", "AR-NF", "LB-NF", "IC-AR"))+
  xlab("Distance (km)")+
  ylab("FST")+
  xlim(2000,5700)+
  geom_text(aes(label=pair_abbrev), color="black",size=3, hjust=-0.3,vjust=0.5)+
  labs(fill="Region pair")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

IBD_noLB
ggsave("IBD_plot_blues_pointlabels_noLB.png", width=7, height=6, dpi=300)


# if want to combine plots:
library(patchwork)
IBD / IBD_noLB + plot_annotation(tag_levels = 'A')

ggsave("IBD_with_and_withoutLB_vertical.png", width=7, height=10, dpi=600)

# run Mantel's test to test correlation of FST & distance
library(ade4)

# load data (distance matrix, and fst matrix)
distance_data <- read.csv("shelf_dist_matrix.csv")
fst_data <- read.csv("autosomes.LDpruned.n36.FST.slat.matrix.csv")

distance_dist <- dist(distance_data)
fst_dist <- dist(fst_data)

distance_mantel <- mantel.rtest(distance_dist, fst_dist, nrepet=9999)
distance_mantel

# run again without LB
distance_data_noLB <- distance_data[-3,]
distance_data_noLB <- distance_data_noLB[,-4]
fst_data_noLB <- fst_data[-3,]
fst_data_noLB <- fst_data_noLB[,-4]

distance_dist_noLB <- dist(distance_data_noLB)
fst_dist_noLB <- dist(fst_data_noLB)

distance_mantel_noLB <- mantel.rtest(distance_dist_noLB, fst_dist_noLB, nrepet=9999)
distance_mantel_noLB
