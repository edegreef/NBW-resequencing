# plot SMC++ models using the .csv outputs from 'smc++ plot' function

library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/demography/smcpp2/min100kb")

# final plots used were from 1) 1-group with all samples and 2) 3-group divided by AR+IC,LB+NF, and SS 
# (starting around line 106 & 222)

#### plotting 5-group (each subregion separated)
AR <- read.csv("5_group_seppop/smc2_AR_plot.csv", header=T)
IC <- read.csv("5_group_seppop/smc2_IC_plot.csv", header=T)
LB <- read.csv("5_group_seppop/smc2_LB_plot.csv", header=T)
NF <- read.csv("5_group_seppop/smc2_NF_plot.csv", header=T)
SS <- read.csv("5_group_seppop/smc2_SS_plot.csv", header=T)

# merge
runs <- rbind(AR, IC, LB, NF, SS)

# add unique id for label+run
runs$plot_id <- paste(runs$label, runs$plot_num, sep="_")

# set generation time (and mutation rate for axis label)
gen=17.8
mu="1.53x10^-8"

# set colors
display.brewer.all(colorblindFriendly = T)
display.brewer.pal(n = 10, name = 'RdYlBu')
brewer.pal(n = 10, name = 'RdYlBu')

group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")

# adjust order to make legend by latitude
runs$label <- factor(runs$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# quick plot in ggplot
ggplot()+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label))+
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

# merge
MED <- rbind(AR_med, IC_med, LB_med, NF_med, SS_med)

# adjust label order by latitude (for legend)
MED$label <- factor(MED$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# plot all runs and median in bold in ggplot
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
#ggsave("SMC_plot_5group_2_noLGM.png", width=7,height=4,dpi=1000)

#expression(Pairwise~italic(r)^{2})

#### next doing 1_group_ALL files (when SMC was estimated with all samples, but I tried 5 runs to use each pop as derived thing so using 5 data files again)

AR <- read.csv("1_group_ALL/smc2_ALL.AR_plot.csv", header=T)
IC <- read.csv("1_group_ALL/smc2_ALL.IC_plot.csv", header=T)
LB <- read.csv("1_group_ALL/smc2_ALL.LB_plot.csv", header=T)
NF <- read.csv("1_group_ALL/smc2_ALL.NF_plot.csv", header=T)
SS <- read.csv("1_group_ALL/smc2_ALL.SS_plot.csv", header=T)

AR$label <- "AR"
IC$label <- "IC"
LB$label <- "LB"
NF$label <- "NF"
SS$label <- "SS"

# merge
runs2 <- rbind(AR, IC, LB, NF, SS)

# add unique id for label+run
runs2$plot_id <- paste(runs2$label, runs2$plot_num, sep="_")

# set generation time (and mutation rate for axis label)
gen=17.8
mu=1.53e-8

group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")

# adjusting order to make legend by latitude
runs2$label <- factor(runs2$label, levels = c("IC", "AR", "LB", "NF", "SS"))

# quick plot in ggplot
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


# extract median run per for all
#MED <- runs %>% group_by(plot_num) %>% mutate(row2=row_number())
AR <- AR %>% group_by(plot_num) %>% mutate(row=row_number())
IC <- IC %>% group_by(plot_num) %>% mutate(row=row_number())
LB <- LB %>% group_by(plot_num) %>% mutate(row=row_number())
NF <- NF %>% group_by(plot_num) %>% mutate(row=row_number())
runs <- rbind(AR, IC, LB, NF, SS)

MED <- runs2 %>% group_by(plot_num) %>% mutate(row2=row_number())
MED_med <- aggregate(MED[,c("x","y")], list(MED$row2), FUN=median)
colnames(MED_med)[1] <- "row"
MED_med$label <- "ALL"

MED <- MED_med

# for MED thing to work on its own need to make sure end result is 101 not 505 (weird plot was based on row num issue)
MED$group <- rep(1:101, times=5)
MED <- aggregate(MED[,c("x","y")], list(MED$group), FUN=median)
colnames(MED)[1] <- "row"
MED$label <- "ALL"

col <- (ALL="#88419D")

# plot all runs and median in bold in ggplot
smc_all<-ggplot()+
    geom_step(data=runs2, aes(x=(x*gen),y=y, group=plot_id),color="white", alpha=0)+
annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  annotate("rect", xmin = 100, xmax =500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  #geom_vline(xintercept = 500, linetype="dashed", color="gray")+
  #geom_vline(xintercept = 200, linetype="dashed", color="gray")+
  geom_step(data=runs2, aes(x=(x*gen),y=y, group=plot_id),color="#8C6BB1", alpha=0.2)+
 geom_step(data=MED, aes(x=(x*gen),y=y, group=label, color=label), lwd=2)+
  scale_color_manual(values=col)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA, colour = "black", size=1))+ 
  annotation_logticks(sides = "lbt")+
  scale_y_continuous(trans='log10', breaks=c(100,300,1000,3000,10000, 30000), labels=c("100"="0.1","300"="0.3", "1000"="1 ", "3000"="3 ","10000"="10","30000"="30"), limits=c(100,20000), expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(100,300000), expand=c(0,0))+
  xlab(expression(paste("Years ago")))+# (g=17.8, mu=1.53x")~10^{-8}~")"))+
  labs(color="Region")
  #ggtitle("Estimated as one population")

smc_all
#ggsave("SMC_plot_1_group_ALL_purp_properlabel.png", width=7,height=4,dpi=1000)

#smc_all / smc
#ggsave("SMC_plots_ALL_5group_2_noLGM.png", width=7,height=7.5,dpi=1000)


#### 3-group (AR+IC, LB+NF, SS)

ARIC <- read.csv("3_group/smc_ARIC.AR_plot.csv", header=T)
LBNF <- read.csv("3_group/smc_LBNF.NF_plot.csv", header=T)
SS <- read.csv("5_group_seppop/smc2_SS_plot.csv", header=T)

# merge
runs <- rbind(ARIC, LBNF, SS)

# add unique id for label+run
runs$plot_id <- paste(runs$label, runs$plot_num, sep="_")

# set generation time (and mutation rate for axis label)
gen=17.8
mu=1.53e-8

group_colors <- c(ARIC="#4575B4", LBNF="#F46D43", SS="#A50026")

# adjusting order to make legend by latitude
runs$label <- factor(runs$label, levels = c("ARIC", "LBNF", "SS"))

# quick plot in ggplot
ggplot()+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label))+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  annotation_logticks(sides = "b")+
  theme_classic()+ylab("Ne")+
  xlab(paste("Years ago (g=",gen,", mu=",mu,")", sep=""))

# extract median run per group
ARIC <- ARIC %>% group_by(plot_num) %>% mutate(row=row_number())
ARIC_med <- aggregate(ARIC[,c("x","y")], list(ARIC$row), FUN=median)
colnames(ARIC_med)[1] <- "row"
ARIC_med$label <- "ARIC"

LBNF <- LBNF %>% group_by(plot_num) %>% mutate(row=row_number())
LBNF_med <- aggregate(LBNF[,c("x","y")], list(LBNF$row), FUN=median)
colnames(LBNF_med)[1] <- "row"
LBNF_med$label <- "LBNF"

SS <- SS %>% group_by(plot_num) %>% mutate(row=row_number())
SS_med <- aggregate(SS[,c("x","y")], list(SS$row), FUN=median)
colnames(SS_med)[1] <- "row"
SS_med$label <- "SS"

# merge
MED <- rbind(ARIC_med, LBNF_med, SS_med)

# adjust label order by latitude (for legend)
MED$label <- factor(MED$label, levels = c("ARIC", "LBNF", "SS"))

# plot all runs + median in ggplot
smc<-ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray30")+
  annotate("rect", xmin = 100, xmax =500, ymin=0, ymax=20000, alpha = 0.2, fill = "gray60")+
  geom_step(data=runs, aes(x=(x*gen),y=y, group=plot_id, color=label), alpha=0.3)+
  geom_step(data=MED, aes(x=(x*gen),y=y, group=label, color=label), lwd=1.5)+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10')+
  #scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10')+
  annotation_logticks(sides = "lbt")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA, colour = "black", size=1))+ 
  scale_y_continuous(trans='log10', breaks=c(100,300,1000,3000,10000, 30000), labels=c("100"="0.1","300"="0.3", "1000"="1 ", "3000"="3 ","10000"="10","30000"="30"), limits=c(100,20000), expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(100,300000), expand=c(0,0))+
  xlab(expression(paste("Years ago")))+
  labs(color="Region")
smc

#ggsave("SMC_plot_3_group.png", width=7,height=4,dpi=1000)

# saving SMC++ plot with all samples as one group, and SMC++ plot with samples divided into 3 groups/regions
smc_all / smc
ggsave("SMC_all_3group_100mark_ROB.png", width=6.3,height=7,dpi=1000)