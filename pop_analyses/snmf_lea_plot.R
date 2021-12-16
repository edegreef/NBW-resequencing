# Results from SNMF

library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(patchwork)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/admixture/min50kb_miss01/LDprunedr05")

# load sample info 
sample_info <- read.csv("sample_info_pop_order_n37_standardID.csv", header=T)

# order by pop - using this for label in coverage plot
sample_info2 <- sample_info[order(sample_info[,c("pop_order")]),]

# load in snmf results (K2)
snmf_q <- read.csv("K2_qlongmatrix_min50kb_LDprunedr05_miss01.csv", header=T)
snmf_q <- snmf_q[,2:3]

# add individual IDs
snmf_q$Ind = sample_info$pop_order

# convert dataframe to long format
qlong_snmf = melt(snmf_q, id.vars=c("Ind"))
head(qlong_snmf)

# plot snmf barplot 
snmf_admix_K2 = ggplot(data=qlong_snmf, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=2")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.2" = "#FEE090", "Cluster.1" = "#A50026", "Cluster.3" = "#59A14F", "Cluster.4" = "#EDC948", "Cluster.5" = "#BAB0AC"))+  
  # ggtitle("K=2")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), text=element_text(family="serif"))+
  theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())
snmf_admix_K2

# load in snmf results (K3)
snmf_q <- read.csv("K3_qlongmatrix_min50kb_LDprunedr05_miss01.csv", header=T)
snmf_q <- snmf_q[,2:4]
snmf_q$Ind = sample_info$pop_order
qlong_snmf = melt(snmf_q, id.vars=c("Ind"))
head(qlong_snmf)


######### add this bit for sample axis if want to add sample ID to k3 plot, otherwise skip and save for K5.
sample_axis <- as.vector(sample_info2$sample)
##########

# plot snmf barplot 
snmf_admix_K3 = ggplot(data=qlong_snmf, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = FALSE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=3")+
  xlab("Individual")+
  theme_classic()+
  #  ggtitle("K=3")+
  scale_fill_manual(values = c("Cluster.1" = "#FEE090", "Cluster.3" = "#A50026", "Cluster.2" = "#4575B4", "Cluster.4" = "#EDC948", "Cluster.5" = "#BAB0AC"))+    
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), text=element_text(family="serif"))+
  theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())+
 scale_x_continuous(breaks=1:37, labels=sample_axis)+
 theme(axis.text.x=element_text(angle=90))
snmf_admix_K3

snmf_admix_K2/snmf_admix_K3
#ggsave("leaK2-3_min50kb_miss01_LDprunedr05_formega_withyellow.png", width=7, height=2.7, dpi=700)

# load in snmf results (K4)
snmf_q <- read.csv("K4_qlongmatrix_min50kb_LDprunedr05_miss01.csv", header=T)
snmf_q <- snmf_q[,2:5]
snmf_q$Ind = sample_info$pop_order
qlong_snmf = melt(snmf_q, id.vars=c("Ind"))
head(qlong_snmf)

# plot snmf barplot 
snmf_admix_K4 = ggplot(data=qlong_snmf, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = FALSE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=4")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.1" = "#ABD9E9", "Cluster.4" = "#A50026", "Cluster.2" = "#4575B4", "Cluster.3" = "#EDC948", "Cluster.5" = "#BAB0AC"))+
  # ggtitle("K=4")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), text=element_text(family="serif"))+
  theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())
snmf_admix_K4

# load in snmf results (K5)
snmf_q <- read.csv("K5_qlongmatrix_min50kb_LDprunedr05_miss01.csv", header=T)
snmf_q <- snmf_q[,2:6]

# add individual IDs
snmf_q$Ind = sample_info$pop_order

# convert dataframe to long format
qlong_snmf = melt(snmf_q, id.vars=c("Ind"))
head(qlong_snmf)

sample_axis <- as.vector(sample_info2$sample)

# plot snmf barplot 
snmf_admix_K5 = ggplot(data=qlong_snmf, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = FALSE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=5")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.2" = "#ABD9E9", "Cluster.5" = "#A50026", "Cluster.3" = "#4575B4", "Cluster.4" = "#EDC948", "Cluster.1" = "#BAB0AC"))+   
  
 # ("Cluster.1" = "#ABD9E9", "Cluster.4" = "#A50026", "Cluster.2" = "#4575B4", "Cluster.3" = "#EDC948", "Cluster.5" = "#BAB0AC"))+
  #ggtitle("K=5")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), text=element_text(family="serif"))+
  theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())+
  scale_x_continuous(breaks=1:37, labels=sample_axis)+
  theme(axis.text.x=element_text(angle=90))
snmf_admix_K5


snmf_admix_K2 / snmf_admix_K3 / snmf_admix_K4 / snmf_admix_K5
#ggsave("leaK2-5_min50kb_miss01_LDprunedr05_colours.png", width=7, height=5, dpi=600)


# add error plot
lea_er <- read.csv("snmf_lea_CVerror.csv", header=T)

lea <- ggplot(data=lea_er, aes(x=K, y=error))+
  geom_point(colour="darkblue")+
  geom_line(colour="darkblue")+
  theme_classic()  +
  ggtitle("SNMF")+
  ylab("Cross-entropy")+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(family="serif"))
lea

#ggsave("error_K1-5_miss01_LDprunedr05.png", width = 3.5, height = 3, dpi = 400)

