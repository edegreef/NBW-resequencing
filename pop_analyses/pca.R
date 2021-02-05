# First part - running PCA using help from the PCAdapt vignette: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
# Second half - making covariance matrix heatmap

library(pcadapt)
library(ggplot2)
library(patchwork)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/pca")

# load data
snp_data <- read.pcadapt("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.bed", type = "bed")

# load sample info
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info.csv", header=T)

# filter out individuals not in snp data
sample_info_36 <- subset(sample_info, remove_indiv!="Y")

# run pcadapt. K value will be how many eigenvectors to be produced
x <- pcadapt(input = snp_data, K = 20)

# screeplot
plot(x, option = "screeplot")

# plot quick PCA
plot(x, option = "scores", pop = sample_info_36$region)

# plot other eigenvectors
plot(x, option = "scores", i = 5, j = 6, pop = sample_info_36$region)

# look at pca scores
scores <- as.data.frame(x$scores)

# look at loadings
loadings <- as.data.frame(x$loadings)

# z scores
z_scores <- as.data.frame(x$zscores)

# look at proportion variance
proportion <- as.data.frame(x$singular.values)
PC1_proportion <- (round(proportion[1,], digits=4))*100
PC2_proportion <- (round(proportion[2,], digits=4))*100
PC3_proportion <- (round(proportion[3,], digits=4))*100
PC4_proportion <- (round(proportion[4,], digits=4))*100
PC5_proportion <- (round(proportion[5,], digits=4))*100
PC6_proportion <- (round(proportion[6,], digits=4))*100

# save the scores as separate data file to adjust the pca plot for colors, etc
evec <- cbind(sample_info_36$sample, scores)
colnames(evec)[1] <- "sample"
#write.csv(evec,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.evec.csv", row.names=FALSE)

# plot in ggplot
pca <- ggplot(data=evec, aes(x=V1,y=V2))+
  geom_point(aes(fill=sample_info_36$region), colour="black",pch=21,size=4)+#, shape=sample_info$updated_sex)) +
  theme_classic()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  ggtitle("NBW PCA, autosomes, LD pruned, n=36")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC1 (", PC1_proportion, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), 
                    breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))+
  labs(fill= "Region")
  #geom_text(aes(label=sample, alpha=0.5), hjust=0,vjust=1)
  #stat_ellipse(aes(group = sample_info$region, color=sample_info$region))
pca

# save plot
ggsave("PCA_platypus_autosomes_LDpruned_n36_classic_legendorder.png", width=9, height=6, dpi=300)

# to plot other PCs
pca_p3p4 <- ggplot(data=evec, aes(x=V3,y=V4))+
  geom_point(aes(fill=sample_info_36$region), colour="black",pch=21,size=4)+
  theme_bw()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC3 (", PC3_proportion, "%)", sep=""))+
  ylab(paste("PC4 (", PC4_proportion, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"))+
  labs(fill= "Region")
pca_p3p4

pca_p5p6 <- ggplot(data=evec, aes(x=V5,y=V6))+
  geom_point(aes(fill=sample_info_36$region), colour="black",pch=21,size=4)+
  theme_bw()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC5 (", PC5_proportion, "%)", sep=""))+
  ylab(paste("PC6 (", PC6_proportion, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"))+
  labs(fill= "Region")
pca_p5p6

pca_p3p4 + pca_p5p6 + plot_layout(guides='collect') + plot_annotation(title = "NBW PCAs, autosomes, LD pruned, n=36")#,theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("PCA_PC3-6_platypus_autosomes_LDpruned_n36.png", width=12, height=5, dpi=300)


# PCAs on X an Y chromosomes
# load data
Xchr <- read.pcadapt("NBW_platypus_SNPs.filter1.filter2.ID.Xchr.LDpruned.femalesX.bed", type = "bed")
Ychr <- read.pcadapt("NBW_platypus_SNPs.filter1.filter2.ID.Ychr.LDpruned.malesY.bed", type="bed")
Xchr_36 <- read.pcadapt("NBW_platypus_SNPs.filter1.filter2.ID.Xchr.LDpruned.n36.bed", type = "bed")

sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info.csv", header=T)
sample_info_X <- subset(sample_info, remove_indiv!="Y" & updated_sex =="F")
sample_info_Y <- subset(sample_info, remove_indiv!="Y" & updated_sex =="M")
sample_info_36 <- subset(sample_info, remove_indiv!="Y")

# run pcadapt. K value will be how many eigenvectors to be produced
x <- pcadapt(input = Xchr, K = 10)
y <- pcadapt(input = Ychr, K = 10)
x36 <- pcadapt(input = Xchr_36, K=10)

# look at  pca scores
x_scores <- as.data.frame(x$scores)
y_scores <- as.data.frame(y$scores)
x36_scores <- as.data.frame(x36$scores)

# look at proportion variance
proportion_x <- as.data.frame(x$singular.values)
PC1_proportion_x <- (round(proportion_x[1,], digits=4))*100
PC2_proportion_x <- (round(proportion_x[2,], digits=4))*100

proportion_y <- as.data.frame(y$singular.values)
PC1_proportion_y <- (round(proportion_y[1,], digits=4))*100
PC2_proportion_y <- (round(proportion_y[2,], digits=4))*100

proportion_x36 <- as.data.frame(x36$singular.values)
PC1_proportion_x36 <- (round(proportion_x36[1,], digits=4))*100
PC2_proportion_x36 <- (round(proportion_x36[2,], digits=4))*100

# save the scores as separate data file so I can adjust the pca plot for colors, etc
evec_x <- cbind(sample_info_X$sample, x_scores)
colnames(evec_x)[1] <- "sample"
evec_y <- cbind(sample_info_Y$sample, y_scores)
colnames(evec_y)[1] <- "sample"
evec_x36 <- cbind(sample_info_36$sample, x36_scores)
colnames(evec_x36)[1] <- "sample"

# plot in ggplot
x_pca <- ggplot(data=evec_x, aes(x=V1,y=V2))+
  geom_point(aes(fill=sample_info_X$region), colour="black",pch=21,size=4)+
  theme_classic()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  ggtitle("X chr, females, n=16")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC1 (", PC1_proportion_x, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion_x, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), 
                    breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))+
  labs(fill= "Region")
x_pca

y_pca <- ggplot(data=evec_y, aes(x=V1,y=V2))+
  geom_point(aes(fill=sample_info_Y$region), colour="black",pch=21,size=4)+
  theme_classic()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  ggtitle("Y chr, males, n=20")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC1 (", PC1_proportion_y, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion_y, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), 
                    breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))+
  labs(fill= "Region")
y_pca

x36_pca <- ggplot(data=evec_x36, aes(x=V1,y=V2))+
  geom_point(aes(fill=sample_info_36$region), colour="black",pch=21,size=4)+
  theme_classic()+
  theme(panel.border=element_blank(), 
        axis.line=element_line())+
  ggtitle("X chr, females & males, n=36")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(paste("PC1 (", PC1_proportion_x36, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion_x36, "%)", sep=""))+
  scale_fill_manual(values=c("#C51B7D", "#481567FF", "#287D8EFF", "#95D840FF", "#FDE725FF"), 
                    breaks=c("Iceland", "Arctic", "Labrador", "Newfoundland", "Scotian_shelf"))+
  labs(fill= "Region")
x36_pca

x_pca + y_pca + x36_pca + guide_area() + plot_layout(guides='collect')
ggsave("PCA_platypus_XYchrs_LDpruned_multiplot_classic_latorder.png", width=10, height=8, dpi=300)


### trying something out (PCA heatmap)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# PCA output in order by pop
evec_poporder <- read.csv("evec_n36_poporder.csv", header=T)

# edit dataframe so it has sample ID in column 1 and the eigenvalues in other columns.
evec_poporder <- evec_poporder[order(evec_poporder[,2]),]
evec_poporder <- evec_poporder[,-2]
evec_poporder <- evec_poporder[,-2] 

# dst function help from: https://stackoverflow.com/questions/44731442/create-heat-map-from-pca-coordinates-in-r

# Define a distance function based on euclidean norm
# calculated between PCA values of the i-th and j-th items
dst <- Vectorize(function(i,j,dtset) sqrt(sum((dtset[i,2:3]-dtset[j,2:3])^2)), vectorize.args=c("i","j"))

# Here is the distance between echocardiography and infarction
dst(1,2,evec_poporder)

# Calculate the distance matrix
nr <- nrow(evec_poporder)
mtx <- outer(1:nr, 1:nr, "dst", dtset=evec_poporder)
colnames(mtx) <- rownames(mtx) <- evec_poporder[,1]

# change dataframe format to melted matrix
mtx.long <- melt(mtx)

# Plot the heatmap using ggplot2
ggplot(mtx.long, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile()+
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette="RdYlBu")+
  theme(axis.text.x=element_text(angle=90))+
  ggtitle("distance matrix using PC1 and PC2")
ggsave("PC_distance_mat.png", width=9, height=8, dpi=300)

#display.brewer.all(colorblindFriendly=TRUE)


## trying another thing out
# covariance matrix from snp data directly (like relatedness matrix)       
plink <- read.csv("matrix.sq.poporder.rel.csv", row.names = 1, header=T)
plink <- as.matrix(plink)
                     
melt_plink <- melt(plink)

ggplot(melt_plink, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile()+
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette="RdYlBu")+
  theme(axis.text.x=element_text(angle=90))
                           
# same thing but using SNPRelate (ends up with same output)  
library(SNPRelate)
bed.in <- "C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/pca/NBW_snps_auto_LDprune_n36_poporder.bed"

snpgdsBED2GDS(bed.in, out.gdsfn="nbw.gds")
snpgdsSummary("nbw.gds")
                           
open <- snpgdsOpen("nbw.gds")
test <- snpgdsGRM(open, autosome.only=FALSE)
                           
grm <- as.data.frame(test$grm)
grm <- as.matrix(grm)
grm[grm >=0.75] <- NA
melt_grm <- melt(grm)
                           
ggplot(melt_grm, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile()+
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette="RdYlBu")+
  theme(axis.text.x=element_text(angle=90))
                         
