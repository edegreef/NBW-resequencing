# First part - running PCA using help from the PCAdapt vignette: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
# Second half - making covariance matrix heatmap

# Running PCA using help from the PCAdapt vignette: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

library(pcadapt)
library(ggplot2)
library(patchwork)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/pca")

# Load data
snp_data <- read.pcadapt("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.37.bed", type = "bed")

# Load sample info
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/Hyperoodon_popinfo4_reseq_Oct2021.csv", header=T)

# Filter out individuals not not in snp data
sample_info_subset <- subset(sample_info, remove!="Y")

# Run pcadapt. K value will be how many eigenvectors to be produced
x <- pcadapt(input = snp_data, K = 20)

# Screeplot
plot(x, option = "screeplot")

# Plot quick PCA
plot(x, option = "scores", pop = sample_info_subset$region)

# Plot other eigenvectors
plot(x, option = "scores", i = 5, j = 6, pop = sample_info_subset$region)

# Look at pca scores
scores <- as.data.frame(x$scores)

# Look at loadings
loadings <- as.data.frame(x$loadings)

# Z scores
z_scores <- as.data.frame(x$zscores)

# Look at proportion variance
proportion <- as.data.frame(x$singular.values)
PC1_proportion <- (round(proportion[1,], digits=4))*100
PC2_proportion <- (round(proportion[2,], digits=4))*100
PC3_proportion <- (round(proportion[3,], digits=4))*100
PC4_proportion <- (round(proportion[4,], digits=4))*100
PC5_proportion <- (round(proportion[5,], digits=4))*100
PC6_proportion <- (round(proportion[6,], digits=4))*100

# Save the scores as separate data file to adjust the pca plot for colors, etc
evec <- cbind(sample_info_subset$Genome_ID, scores)
colnames(evec)[1] <- "sample"
#write.csv(evec,"NBW2_SNPS_2M_pops_37.evec.csv", row.names=FALSE)

library(ggplot2)
library(ggrepel)

# plot in ggplot
pca <- ggplot(data=evec, aes(x=V1,y=V2))+
  geom_point(aes(color=sample_info_subset$region),size=5, alpha=0.85)+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(family="serif", size=14))+
  theme(panel.border=element_blank(), 
        axis.line=element_line(),
        text=element_text(family="serif", size=14),
        legend.position=c(0.88,0.88), #(0,0)=bottom left, (1,1)=top right
        legend.background=element_rect(fill=NA, color="gray30"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.size=unit(0.8,'lines'))+
  #guides(color=guide_legend(override.aes=list(size=2)))+
  xlab(paste("PC1 (", PC1_proportion, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion, "%)", sep=""))+
  scale_color_manual(values=c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026"), 
                     breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"),
                     labels=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian Shelf"))+
  labs(color= "Region")
  #geom_label_repel(aes(label=sample_info_subset$Genome_ID), max.overlaps = 24)
  
pca

ggsave("PCA_NBW2_SNPS_2M_37_noprune.png", width=6, height=5, dpi=1000)

##################
### Trying something out (PCA heatmap)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# PCA output in order by pop
# evec_poporder <- read.csv("evec_n36_poporder.csv", header=T)
# add pop info to evec data file
#id_pop <- sample_info_subset[,c("Genome_ID","region")]
#evec_pop <- cbind(id_pop, evec)
#evec_pop <-evec_pop[,-1]
#evec_pop <- evec_pop[order(evec_pop$region),]
  
#write.csv(evec_pop,"evec_pop_prep.csv", row.names=FALSE)
# Fix order for AR/IC
evec_poporder <- read.csv("NBW2_SNPS_2M_pops_37_orderprep.evec.csv", header=T)

# Remove pop
evec_prep <- evec_poporder[,-2]

# Edit dataframe so it has sample ID in column 1 and the eigenvalues in other columns.
#evec_poporder <- evec_poporder[order(evec_poporder[,2]),]
#evec_poporder <- evec_poporder[,-2]
#evec_poporder <- evec_poporder[,-2] 

# dst function help from: https://stackoverflow.com/questions/44731442/create-heat-map-from-pca-coordinates-in-r

# Define a distance function based on euclidean norm
# calculated between PCA values of the i-th and j-th items
dst <- Vectorize(function(i,j,dtset) sqrt(sum((dtset[i,2:3]-dtset[j,2:3])^2)), vectorize.args=c("i","j"))

# Here is the distance between echocardiography and infarction
dst(1,2,evec_prep)

# Calculate the distance matrix
nr <- nrow(evec_prep)
mtx <- outer(1:nr, 1:nr, "dst", dtset=evec_prep)
colnames(mtx) <- rownames(mtx) <- evec_prep[,1]

# Change dataframe format to melted matrix
mtx.long <- melt(mtx)

# Plot the heatmap using ggplot2
ggplot(mtx.long, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile()+
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette="RdYlBu")+
  theme(axis.text.x=element_text(angle=90),
        text=element_text(family="serif"))+
  ggtitle("distance matrix using PC1 and PC2")+
  scale_y_discrete(position = "right")
ggsave("PC_distance_mat_NBW2_SNPS_2M_37-IC-AR.png", width=9, height=8, dpi=1000)

