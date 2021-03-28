# Plotting PSMC data with and without bootstrapping
 
library(ggplot2)
library(tidyverse)

# Without bootstrap
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/demography/psmc_plot")

IC <- read.delim("IC_1.53e-8_g17.8_plot.0.txt", header=FALSE)
AR <- read.delim("AR_1.53e-8_g17.8_plot.0.txt", header=FALSE)
LB <- read.delim("LB_1.53e-8_g17.8_plot.0.txt", header=FALSE)
NF <- read.delim("NF_1.53e-8_g17.8_plot.0.txt", header=FALSE)
SS <- read.delim("SS_1.53e-8_g17.8_plot.0.txt", header=FALSE)

# plotting just 1
ggplot(data=IC, aes(x=V1, y=V2)) + 
    geom_step() + 
    scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10',limits=c(1e4,1e6)) + 
    scale_y_continuous(limits=c(0,2))+ 
    #geom_hline(yintercept = 5) + 
    annotation_logticks(sides = "lb") +
    theme_classic()

# plot all 5 samples
ggplot() + 
  geom_step(data=AR, aes(x=V1, y=V2, color="AR"), lwd=1) + 
  geom_step(data=IC, aes(x=V1, y=V2, color="IC"),lwd=1) + 
  geom_step(data=LB, aes(x=V1, y=V2, color="LB"),lwd=1) + 
  geom_step(data=NF, aes(x=V1, y=V2, color="NF"),lwd=1) + 
  geom_step(data=SS, aes(x=V1, y=V2, color="SS"),lwd=1) + 
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10',limits=c(1e4,1e6)) + 
  scale_y_continuous(limits=c(0,2))+ 
  #geom_hline(yintercept = 5) + 
  annotation_logticks(sides = "b") +
  theme_classic()+
  labs(x="Years ago", y="Effective pop size 10^4", color="Legend")+
  ggtitle("psmc, mutation rate 1.53e-8, generation time 17.8")


# With boostrapping files
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/demography/psmc_plot/boot")

# create list of file names to load (100 per group is too much to load individually)
AR_filenames <- list.files(pattern="AR*")
IC_filenames <- list.files(pattern="IC*")
LB_filenames <- list.files(pattern="LB*")
NF_filenames <- list.files(pattern="NF*")
SS_filenames <- list.files(pattern="SS*")

# import them together (by group)
AR_files <- lapply(AR_filenames,function(i){
  read.delim(i, header=FALSE)
})

IC_files <- lapply(IC_filenames,function(i){
  read.delim(i, header=FALSE)
})

LB_files <- lapply(LB_filenames,function(i){
  read.delim(i, header=FALSE)
})

NF_files <- lapply(NF_filenames,function(i){
  read.delim(i, header=FALSE)
})

SS_files <- lapply(SS_filenames,function(i){
  read.delim(i, header=FALSE)
})


# merge files from a large list format into a dataframe
AR_boot <- do.call(rbind, AR_files)
IC_boot <- do.call(rbind, IC_files)
LB_boot <- do.call(rbind, LB_files)
NF_boot <- do.call(rbind, NF_files)
SS_boot <- do.call(rbind, SS_files)

# every 58 rows is new boot. need column to label this as run number.
AR_boot$run <- rep(c(1:101), each=58)
IC_boot$run <- rep(c(1:101), each=58)
LB_boot$run <- rep(c(1:101), each=58)
NF_boot$run <- rep(c(1:101), each=58)
SS_boot$run <- rep(c(1:101), each=58)

AR_boot$pop <- "AR"
IC_boot$pop <- "IC"
LB_boot$pop <- "LB"
NF_boot$pop <- "NF"
SS_boot$pop <- "SS"

# merge
runs <- rbind(AR_boot, IC_boot, LB_boot, NF_boot, SS_boot)

# add unique id for label+run
runs$plot_id <- paste(runs$pop, runs$run, sep="_")

# set generation time (and mutation rate for axis label)
gen=17.8
mu=1.53e-8

# set colors
group_colors <- c(IC="#4575B4",AR="#ABD9E9", LB="#FEE090", NF="#F46D43", SS="#A50026")

# adjusting order to make legend by latitude
runs$pop <- factor(runs$pop, levels = c("IC", "AR", "LB", "NF", "SS"))

# extract median run per group
AR <- AR_boot %>% group_by(run) %>% mutate(row=row_number())
AR_med <- aggregate(AR[,c("V1","V2")], list(AR$row), FUN=median)
colnames(AR_med)[1] <- "row"
AR_med$pop <- "AR"

IC <- IC_boot %>% group_by(run) %>% mutate(row=row_number())
IC_med <- aggregate(IC[,c("V1","V2")], list(IC$row), FUN=median)
colnames(IC_med)[1] <- "row"
IC_med$pop <- "IC"

LB <- LB_boot %>% group_by(run) %>% mutate(row=row_number())
LB_med <- aggregate(LB[,c("V1","V2")], list(LB$row), FUN=median)
colnames(LB_med)[1] <- "row"
LB_med$pop <- "LB"

NF <- NF_boot %>% group_by(run) %>% mutate(row=row_number())
NF_med <- aggregate(NF[,c("V1","V2")], list(NF$row), FUN=median)
colnames(NF_med)[1] <- "row"
NF_med$pop <- "NF"

SS <- SS_boot %>% group_by(run) %>% mutate(row=row_number())
SS_med <- aggregate(SS[,c("V1","V2")], list(SS$row), FUN=median)
colnames(SS_med)[1] <- "row"
SS_med$pop <- "SS"

# merge
MED <- rbind(AR_med, IC_med, LB_med, NF_med, SS_med)

# adjust label order by latitude (for legend)
MED$pop <- factor(MED$pop, levels = c("IC", "AR", "LB", "NF", "SS"))

# plot all runs & median highlighted in ggplot
# the annotate("rect"..) is adding in the LGM (last glacial maximum) and LGP (last glacial period) timeline

psmc_boot <-ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=100000, alpha = 0.2, fill = "gray60")+
  annotate("rect", xmin = 19000, xmax =26500, ymin=0, ymax=100000, alpha = 0.2, fill = "gray30")+
  geom_step(data=runs, aes(x=V1,y=V2, group=plot_id, color=pop), alpha=0.1)+
  geom_step(data=MED, aes(x=V1,y=V2, group=pop, color=pop), lwd=0.8)+
  scale_color_manual(values=group_colors)+
  scale_x_continuous(trans='log10',breaks=c(1e3,1e4,1e5,1e6))+#,limits = c(3e3,2e6))+
  #scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10')+# limits = c(0.1,50))+
  annotation_logticks(sides = "lb")+
  theme_classic()+ylab("Effective population size (10^4)")+
  xlab(paste("Years ago (g=",gen,", mu=",mu,")", sep=""))

psmc_boot
ggsave("plot/PSMC_boot_editline2_fullY_xlim_plot_colors_LGM.png", width=7,height=4,dpi=1000)
