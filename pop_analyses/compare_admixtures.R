# Comparing results from admixture analyses (NGSadmix, SNMF, ADMIXTURE) also adding in sample coverage info

library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(patchwork)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/compare_plots/newsnps")

# load sample info 
sample_info <- read.csv("sample_info_n36.csv", header=T)

# order by pop - using this for label in coverage plot
sample_info2 <- sample_info[order(sample_info[,3]),]

# modal coverage
bar_mod <- ggplot(data=sample_info2, aes(x=factor(pop_order, labels=vcf_order), y=updated_modal_coverage))+
  geom_bar(stat="identity", fill="gray")+
  theme_classic()+
  xlab("individual")+
  ylab("Coverage")+
  ggtitle("individual modal coverage")+
  theme(plot.title=element_text(hjust=0.5),axis.text.x=element_blank(), axis.title.x=element_blank())+
  coord_cartesian(expand=FALSE)
  #scale_fill_manual( values = c( "no"="gray", "yes"="orange" ), guide = FALSE )
bar_mod


# load in NGSadmix results
qmatrix <- read.table("NBW_GL_K2_MAC2_n36.qopt")

# label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# add individual IDs
qmatrix$Ind = sample_info$pop_order

# convert dataframe to long format
qlong_ngs = melt(qmatrix, id.vars=c("Ind"))
head(qlong_ngs)

# plot admixture barplot 
NGS_admix <- ggplot(data=qlong_ngs, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Admixture")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster 1" = "lightblue", "Cluster 2" = "tomato"))+
  ggtitle("NGSadmix")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())
NGS_admix


# load in snmf results
snmf_q <- read.csv("SNMF_K2_qmatrix_platypus_auto_LDpruned_n36_miss01.csv", header=T)
snmf_q <- snmf_q[,1:2]

# add individual IDs
snmf_q$Ind = sample_info$pop_order

# convert dataframe to long format
qlong_snmf = melt(snmf_q, id.vars=c("Ind"))
head(qlong_snmf)

# plot snmf barplot 
snmf_admix = ggplot(data=qlong_snmf, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Admixture")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.1" = "tomato", "Cluster.2" = "lightblue"))+
  ggtitle("SNMF")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())
snmf_admix


# load admixture results
admixture_Q <- read.table("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.2.Q", header=F)

# label column names of qmatrix
ncol(admixture_Q)
cluster_names = c()
for (i in 1:ncol(admixture_Q)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(admixture_Q) = cluster_names
head(admixture_Q)

# add individual IDs
admixture_Q$Ind = sample_info$pop_order

# convert dataframe to long format
qlong_admix = melt(admixture_Q, id.vars=c("Ind"))
head(qlong_admix)

#prep sample name for x-axis
sample_axis <- as.vector(sample_info$sample)

# plot admixture barplot 
admix_admix = ggplot(data=qlong_admix, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Admixture")+
  xlab("Individual")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster 1" = "tomato", "Cluster 2" = "lightblue"))+
  ggtitle("ADMIXTURE")+
  theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  scale_x_continuous(breaks=1:36, labels=sample_axis)+
  theme(axis.text.x=element_text(angle=90))
admix_admix

# merge plots into one
bar_mod / NGS_admix/ snmf_admix / admix_admix
ggsave("compare_admix_plots_n36_labeled.png", width = 9, height = 9, dpi = 300)

