# Making fst heatmap and examining isolation-by-distance pattern with mantel test

# 1) fst heatmap
library(reshape2)
library(ggplot2)
library(extrafont)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/fst")

# load in fst matrix
fst <- read.csv("Reichs_fst_output_NBW2_SNPS_2M_plus_up.csv", header=T)

# edit column and row names
colnames(fst) <- c("X", "JM", "AR", "LB", "NF", "SS")
fst$X <- c("JM", "AR", "LB", "NF", "SS")

# convert matrix to melted dataframe to use in geom_tile
melt_fst <- melt(fst)

# define order for plot axis
x_pop_order <- c("SS", "NF", "LB", "AR", "JM")
y_pop_order <-  c("SS", "NF", "LB", "AR", "JM")

fst_heatmap <- ggplot(melt_fst, aes(x=X, y=variable, fill=value)) + 
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
  geom_text(aes(label = round(value, 3)))+
  labs(fill="FST")

fst_heatmap
#ggsave("reich_fst_heatmap_round3.png", width=6, height=5, dpi=500)


# 2) Isolation-by-distance
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(RColorBrewer)
library(ggrepel)

# load in data (in data, each pop pair is a separate row, with FST values and distance as columns)
data <- read.csv("Reichs_fst_output_NBW2_SNPS_2M_plus_up.csv", header=T)
data$pair

# looking at color palettes
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(9, "Blues")
brewer.pal(9, "Blues")

# manually setting color palette (using 'Blues' + one black)
pal <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "black")

# plot IBD with distance
IBD <- ggplot(data, aes(x=shelf_distance_km,y=slatkin_fst))+
  geom_smooth(method = "lm", se=TRUE, fill="gray85", formula=y~x, linetype=0, color="gray")+
  geom_pointrange(aes(ymin=slatkin_ll, ymax=slatkin_ul), size=1, color="black")+
  geom_point(aes(fill=pair), colour="black",pch=21,size=4) +
  theme_classic()+
  ggtitle("Isolation by distance")+
  annotate("text", x=4500, y=0.002, label= "Mantel R2 = 0.80", size=4, fontface="italic", color="black", family="serif")+
  annotate("text", x=4700, y=0.001, label= "p =0.033", size=4, fontface="italic", color="black", family="serif")+
  scale_fill_manual(values=rev(pal), breaks=data$pair)+
  xlab("Distance (km)")+
  ylab(expression(paste(italic("F")~""[ST])))+
  xlim(1000,5000)+
  scale_y_continuous(breaks=c(0.001, 0.003, 0.005, 0.007, 0.009, 0.011))+
 # geom_text(aes(label=pair), color="black",size=4, hjust=-0.3,vjust=0.5, family="serif")+
  #geom_label_repel(aes(label=pair),family="serif", max.overlaps=24)+
  labs(fill="Region pair")+
  theme(plot.title = element_text(hjust = 0), legend.position = "none",
        text=element_text(family="serif", size=12))
IBD
#ggsave("IBD_plot_reich.png", width=8, height=6, dpi=500)
ggsave("IBD_plot_reich_adjust_up_nolabel.png", width=6, height=5, dpi=1000)

# plot IBD with latitud diff (IBL?)
IBL <- ggplot(data, aes(x=latitudinal_diff,y=slatkin_fst))+
  geom_smooth(method = "lm", se=TRUE, fill="gray85", formula=y~x, linetype=0, color="gray")+
  geom_pointrange(aes(ymin=slatkin_ll, ymax=slatkin_ul), size=1, color="black")+
  geom_point(aes(fill=pair), colour="black",pch=21,size=4) +
  theme_classic()+
  ggtitle("Isolation by latitude")+
  annotate("text", x=24, y=0.002, label= "Mantel R2 = 0.62", size=4, fontface="italic", color="black", family="serif")+
  annotate("text", x=25, y=0.001, label= "p =0.026", size=4, fontface="italic", color="black", family="serif")+
  scale_fill_manual(values=rev(pal), breaks=data$pair)+
  xlab("Latitude between sites")+
  ylab(expression(paste(italic("F")~""[ST])))+
  #xlim(3,28)+
#  ylim(0,0.012)+
  scale_x_continuous(breaks=c(5,10,15,20,25), limits=c(3,28))+
#   geom_text(aes(label=pair), color="black",size=4, hjust=-0.3,vjust=0.5, family="serif")+
 # geom_label_repel(aes(label=pair),family="serif", max.overlaps=24)+
  
  labs(fill="Region pair")+
  theme(plot.title = element_text(hjust = 0), legend.position = "none",
        text=element_text(family="serif", size=12))
IBL
#ggsave("IBL_plot_reich.png", width=8, height=6, dpi=300)
ggsave("IBL_plot_reich_up_nolabel.png", width=6, height=5, dpi=1000)


# run Mantel's test to test correlation of FST & distance and FST & latitude
library(ade4)

# load data (distance matrix, and fst matrix)
distance_data <- read.csv("shelf_dist_matrix_500m_marmap_up.csv")
latitude_data <- read.csv("lat_diff_matrix_up.csv")
fst_data <- read.csv("fst_reich_slat_matrix.csv")

distance_dist <- dist(distance_data)
latitude_dist <- dist(latitude_data)
fst_dist <- dist(fst_data)

# fst and distance
distance_mantel <- mantel.rtest(distance_dist, fst_dist, nrepet=9999)
distance_mantel

# fst and latitude
latitude_mantel <- mantel.rtest(latitude_dist, fst_dist, nrepet=9999)
latitude_mantel
