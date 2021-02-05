# running t-SNE 

library(adegenet)
library(Rtsne)
library(tidyverse)
library(vcfR)
library(patchwork)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/tsne")

# load in snps
all_snps <- read.vcfR("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.miss01.imputed.vcf")

# load sample info
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info.csv", header=T)
sample_info_36 <- subset(sample_info, remove_indiv!="Y")
poplist.names <- sample_info_36$region

# add pop info to genlight object
all_snps <- vcfR2genlight(all_snps)
all_snps@pop <- as.factor(poplist.names)

# running rtsne, first trying low perplexity value 
all_snps_tsne_plex2 <- Rtsne(t(na.omit(t((as.matrix(all_snps))))), dims = 2, initial_dims = 100, perplexity = 2, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 5000, verbose = TRUE)

# trying out perplexity 5
all_snps_tsne_plex5 <- Rtsne(t(na.omit(t((as.matrix(all_snps))))), dims = 2, initial_dims = 100, perplexity = 5, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 5000, verbose = TRUE)

# trying out perplexity 10 (highest I can go is 11 with this sample size)
all_snps_tsne_plex10 <- Rtsne(t(na.omit(t((as.matrix(all_snps))))), dims = 2, initial_dims = 100, perplexity = 10, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 5000, verbose = TRUE)

# quick plot
#plot(all_snps_tsne$Y)

# create plots of t-SNE results

# first run
perplex_2 <- ggplot(as.data.frame(all_snps_tsne_plex2$Y), aes(x=all_snps_tsne_plex2$Y[,1], y=all_snps_tsne_plex2$Y[,2], colour = factor(poplist.names))) + 
  geom_point(size = 4) + 
  labs(x = "Dimension 1", y = "Dimension 2") +
  labs(color = "Region")+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  ggtitle("perplexity: 2")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))
perplex_2

# second run
perplex_5 <- ggplot(as.data.frame(all_snps_tsne_plex5$Y), aes(x=all_snps_tsne_plex5$Y[,1], y=all_snps_tsne_plex5$Y[,2], colour = factor(poplist.names))) + 
  geom_point(size = 4) + 
  labs(x = "Dimension 1", y = "Dimension 2") +
  labs(color = "Region")+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  ggtitle("perplexity: 5")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))
perplex_5

# third run
perplex_10 <- ggplot(as.data.frame(all_snps_tsne_plex10$Y), aes(x=all_snps_tsne_plex10$Y[,1], y=all_snps_tsne_plex10$Y[,2], colour = factor(poplist.names))) + 
  geom_point(size = 4) + 
  labs(x = "Dimension 1", y = "Dimension 2") +
  labs(color = "Region")+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  ggtitle("perplexity: 10")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))

perplex_10

# combine into one figure
perplex_2 / perplex_5 / perplex_10 + plot_layout(guides = 'collect') + plot_annotation(title="t-SNE using autosomes, LD pruned, max miss 0.1, imputed for the rest, n=36", theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("tsne_perplex2_5_10_latorder.png", width=7, height=12, dpi=300)
