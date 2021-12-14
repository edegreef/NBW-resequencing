# Plot SMC++ models using the .csv outputs from 'smc++ plot' function

library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/smcpp")

# Doing 5-group first
AR <- read.csv("smc_AR_plot.csv", header=T)
IC <- read.csv("smc_IC_plot.csv", header=T)
LB <- read.csv("smc_LB_plot.csv", header=T)
NF <- read.csv("smc_NF_plot.csv", header=T)
SS <- read.csv("smc_SS_plot.csv", header=T)

# Merge
runs <- rbind(AR, IC, LB,NF, SS)

# Add unique id for label+run
runs$plot_id <- paste(runs$label, runs$plot_num, sep="_")

# Set generation time (and mutation rate for axis label)
gen=17.8
mu=paste("1.53x",expression(10^{-8}), sep="")

# Set colors
display.brewer.all(colorblindFriendly = T)
display.brewer.pal(n = 10, name = 'RdYlBu')
brewer.pal(n = 10, name = 'RdYlBu')

group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")

# Adjusting order to make legend by latitude
runs$label <- factor(runs$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# Quick test plot in ggplot
ggplot()+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label))+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  annotation_logticks(sides = "b")+
  theme_classic()+ylab("Ne")+
  xlab(paste("Years ago (g=",gen,", mu=",mu,")", sep=""))

# Extract median run per group
AR <- AR %>% group_by(plot_num) %>% mutate(row=row_number())
AR_med <- aggregate(AR[,c("x","y")], list(AR$row), FUN=median)
colnames(AR_med)[1] <- "row"
AR_med$label <- "AR"

IC <- IC %>% group_by(plot_num) %>% mutate(row=row_number())
IC_med <- aggregate(IC[,c("x","y")], list(IC$row), FUN=median)
colnames(IC_med)[1] <- "row"
IC_med$label <- "IC"

LB <- LB %>% group_by(plot_num) %>% mutate(row=row_number())
LB_med <- aggregate(LB[,c("x","y")], list(LB$row), FUN=median)
colnames(LB_med)[1] <- "row"
LB_med$label <- "LB"

NF <- NF %>% group_by(plot_num) %>% mutate(row=row_number())
NF_med <- aggregate(NF[,c("x","y")], list(NF$row), FUN=median)
colnames(NF_med)[1] <- "row"
NF_med$label <- "NF"

SS <- SS %>% group_by(plot_num) %>% mutate(row=row_number())
SS_med <- aggregate(SS[,c("x","y")], list(SS$row), FUN=median)
colnames(SS_med)[1] <- "row"
SS_med$label <- "SS"

# Merge
MED <- rbind(AR_med, IC_med, LB_med, NF_med, SS_med)

# Adjust label order by latitude (for legend)
MED$label <- factor(MED$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# Plot all runs + median in ggplot, also log scaling the y axis
smc<-ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
 # annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  geom_vline(xintercept = 500, linetype="dashed", color="gray")+
  geom_vline(xintercept = 200, linetype="dashed", color="gray")+
 geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label), alpha=0.3)+
geom_step(data=MED, aes(x=(x*gen),y=y, group=label, color=label), lwd=1.2)+
  scale_color_manual(values=group_colors)+
  annotation_logticks(sides = "lbt")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA, colour = "black", size=1))+ 
  labs(color="Region")+
  scale_y_continuous(trans='log10', breaks=c(100,300,1000,3000,10000, 30000), labels=c("100"="0.1","300"="0.3", "1000"="1 ", "3000"="3 ","10000"="10","30000"="30"), limits=c(100,20000), expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  
  scale_x_continuous(breaks=c(1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(100,300000), expand=c(0,0))+
  xlab(expression(paste("Years ago")))# (g=17.8, mu=1.53x")~10^{-8}~")"))+
 # ggtitle("Estimated by subregions")
smc

ggsave("SMC_plot_5group.png", width=7,height=4,dpi=1000)


############### Next doing 1_group_ALL files (when SMC was estimated with all samples, but I tried 5 runs to use each pop as distinguished lineage thing)

AR <- read.csv("smc_ALLAR_plot.csv", header=T)
IC <- read.csv("smc_ALLIC_plot.csv", header=T)
LB <- read.csv("smc_ALLLB_plot.csv", header=T)
NF <- read.csv("smc_ALLNF_plot.csv", header=T)
SS <- read.csv("smc_ALLSS_plot.csv", header=T)

AR$label <- "AR"
IC$label <- "IC"
LB$label <- "LB"
NF$label <- "NF"
SS$label <- "SS"

# Merge
runs2 <- rbind(AR, IC, LB, NF, SS)

# Add unique id for label+run
runs2$plot_id <- paste(runs2$label, runs2$plot_num, sep="_")

# Set generation time (and mutation rate for axis label)
gen=17.8
mu=1.53e-8

group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")

# Adjust order to make legend by latitude
runs2$label <- factor(runs2$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# Quick test plot in ggplot
ggplot()+
  geom_step(data=runs2, aes(x=(x*gen),y=y, group=plot_id, color=label))+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  annotation_logticks(sides = "b")+
  theme_classic()+ylab("Ne")+
  xlab(paste("Years ago (g=",gen,", mu=",mu,")", sep=""))

# extract median run per group
AR <- AR %>% group_by(plot_num) %>% mutate(row=row_number())
AR_med <- aggregate(AR[,c("x","y")], list(AR$row), FUN=median)
colnames(AR_med)[1] <- "row"
AR_med$label <- "AR"

IC <- IC %>% group_by(plot_num) %>% mutate(row=row_number())
IC_med <- aggregate(IC[,c("x","y")], list(IC$row), FUN=median)
colnames(IC_med)[1] <- "row"
IC_med$label <- "IC"

LB <- LB %>% group_by(plot_num) %>% mutate(row=row_number())
LB_med <- aggregate(LB[,c("x","y")], list(LB$row), FUN=median)
colnames(LB_med)[1] <- "row"
LB_med$label <- "LB"

NF <- NF %>% group_by(plot_num) %>% mutate(row=row_number())
NF_med <- aggregate(NF[,c("x","y")], list(NF$row), FUN=median)
colnames(NF_med)[1] <- "row"
NF_med$label <- "NF"

SS <- SS %>% group_by(plot_num) %>% mutate(row=row_number())
SS_med <- aggregate(SS[,c("x","y")], list(SS$row), FUN=median)
colnames(SS_med)[1] <- "row"
SS_med$label <- "SS"


# Extract median run per for all?
#MED <- runs %>% group_by(plot_num) %>% mutate(row2=row_number())
AR <- AR %>% group_by(plot_num) %>% mutate(row=row_number())
IC <- IC %>% group_by(plot_num) %>% mutate(row=row_number())
LB <- LB %>% group_by(plot_num) %>% mutate(row=row_number())
NF <- NF %>% group_by(plot_num) %>% mutate(row=row_number())
SS <- SS %>% group_by(plot_num) %>% mutate(row=row_number())
runs <- rbind(AR, IC, LB, NF, SS)

MED <- runs2 %>% group_by(plot_num) %>% mutate(row2=row_number())
MED_med <- aggregate(MED[,c("x","y")], list(MED$row2), FUN=median)
colnames(MED_med)[1] <- "row"
MED_med$label <- "ALL"

MED <- MED_med

MED$group <- rep(1:101, times=5)
MED <- aggregate(MED[,c("x","y")], list(MED$group), FUN=median)
colnames(MED)[1] <- "row"
MED$label <- "ALL"

# Color for ALL
col <- (ALL="#88419D")

# Plot all runs + median in ggplot, also log scaling the y axis
smc_all<-ggplot()+
    geom_step(data=runs2, aes(x=(x*gen),y=y, group=plot_id),color="white", alpha=0)+
annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  annotate("rect", xmin = 100, xmax =500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 100, xmax =200, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  #geom_vline(xintercept = 500, linetype="dashed", color="gray")+
  #geom_vline(xintercept = 200, linetype="dashed", color="gray")+
  geom_step(data=runs2, aes(x=(x*gen),y=y, group=plot_id),color="#8C6BB1", alpha=0.2)+
 geom_step(data=MED, aes(x=(x*gen),y=y, group=label, color=label), lwd=2)+
  scale_color_manual(values=col)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA, colour = "black", size=1),text=element_text(family="serif", size=12))+ 
  annotation_logticks(sides = "lbt")+
  scale_y_continuous(trans='log10', breaks=c(100,300,1000,3000,10000, 30000), labels=c("100"="0.1","300"="0.3", "1000"="1 ", "3000"="3 ","10000"="10","30000"="30"), limits=c(300,20000), expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(100,300000), expand=c(0,0))+
  xlab(expression(paste("Years ago")))+# (g=17.8, mu=1.53x")~10^{-8}~")"))+
  labs(color="Region")
  #ggtitle("Estimated as one population")

smc_all
#ggsave("SMC_plot_1group_ALL2_noLGM.png", width=7,height=4,dpi=1000)
#ggsave("SMC_plot_1_group_ALL_purp_properlabel.png", width=7,height=4,dpi=1000)

#smc_all / smc
#ggsave("SMC_plots_ALL_5group_2_noLGM.png", width=7,height=7.5,dpi=1000)

############## 3/4-group (AR, IC, (optionally)LB+NF, SS)
AR <- read.csv("smc_ARLBNF_group_plot.csv", header=T)
IC <- read.csv("smc_IC_plot.csv", header=T)
#LBNF <- read.csv("smc_LBNFNF_plot.csv", header=T)
SS <- read.csv("smc_SS_plot.csv", header=T)

# Merge
runs <- rbind(IC,AR, SS)

# add unique id for label+run
runs$plot_id <- paste(runs$label, runs$plot_num, sep="_")

# set generation time (and mutation rate for axis label)
gen=17.8
mu=1.53e-8

#group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")
group_colors <- c(IC="#4575B4",ARLBNF="#ABD9E9", SS="#A50026")

# Adjust order to make legend by latitude
runs$label <- factor(runs$label, levels = c("IC", "ARLBNF", "SS"))

# Quick test plot in ggplot
ggplot()+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label))+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  annotation_logticks(sides = "b")+
  theme_classic()+ylab("Ne")+
  xlab(paste("Years ago (g=",gen,", mu=",mu,")", sep=""))

# Extract median run per group
IC <- IC %>% group_by(plot_num) %>% mutate(row=row_number())
IC_med <- aggregate(IC[,c("x","y")], list(IC$row), FUN=median)
colnames(IC_med)[1] <- "row"
IC_med$label <- "IC"

AR <- AR %>% group_by(plot_num) %>% mutate(row=row_number())
AR_med <- aggregate(AR[,c("x","y")], list(AR$row), FUN=median)
colnames(AR_med)[1] <- "row"
AR_med$label <- "ARLBNF"

#LBNF <- LBNF %>% group_by(plot_num) %>% mutate(row=row_number())
#LBNF_med <- aggregate(LBNF[,c("x","y")], list(LBNF$row), FUN=median)
#colnames(LBNF_med)[1] <- "row"
#LBNF_med$label <- "LBNF"

SS <- SS %>% group_by(plot_num) %>% mutate(row=row_number())
SS_med <- aggregate(SS[,c("x","y")], list(SS$row), FUN=median)
colnames(SS_med)[1] <- "row"
SS_med$label <- "SS"

# Merge
MED <- rbind(IC_med,AR_med, SS_med)

# adjust label order by latitude (for legend)
MED$label <- factor(MED$label, levels = c("IC","ARLBNF", "SS"))

# Plot all runs + median in ggplot, also log scaling the y axis
smc<-ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  annotate("rect", xmin = 100, xmax =500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 100, xmax =200, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label), alpha=0.3)+
  geom_step(data=MED, aes(x=(x*gen),y=y, group=label, color=label), lwd=1.5)+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  #scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10')+
  annotation_logticks(sides = "lbt")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA, colour = "black", size=1),text=element_text(family="serif", size=12))+ 
  scale_y_continuous(trans='log10', breaks=c(100,300,1000,3000,10000, 30000), labels=c("100"="0.1","300"="0.3", "1000"="1 ", "3000"="3 ","10000"="10","30000"="30"), limits=c(300,20000), expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(100,300000), expand=c(0,0))+
  xlab(expression(paste("Years ago")))+
  labs(color="Region")
smc

ggsave("SMC_plot_3group_IC_ARLBNF_SS.png", width=7,height=4,dpi=1000)

smc_all / smc
ggsave("SMC_1group_3group_IC_ARLBNF_SS.png", width=5.5,height=6,dpi=1000)
