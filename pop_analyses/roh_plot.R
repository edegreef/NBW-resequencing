# plotting roh using outputs from PLINK
# ran ROH for each population/group separately
# for one group:
# /home/degreefe/programs/plink --allow-extra-chr --bfile snp.scotian_shelf --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 50 --homozyg-snp 50 --out SS

# then plot in R!
library(ggplot2)
setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/roh/")

IC <- read.table("ICtest.hom.indiv", header=T)
ARLBNF <- read.table("ARLBNFtest.hom.indiv", header=T)
SS <- read.table("SStest.hom.indiv", header=T)

IC$region <- "IC"
ARLBNF$region <- "ARLBNF"
SS$region <- "SS"

merge <- rbind(IC, ARLBNF, SS)
#merge <- rbind(IC, SS)

set <- merge[, c("NSEG", "KB", "region")]

set$MB <- (set$KB)/1000

roh <- ggplot(data=set, aes(x=MB,y=NSEG))+
  geom_point(aes(color=region), size=4, alpha=0.7)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (Mb)")+
  ylab("Number of ROHs")+
  theme(legend.position="none", text=element_text(family="serif"))+
  scale_colour_manual(values = c("#ABD9E9", "#4575B4","#A50026"))

roh

bar_NSEG <- ggplot(set, aes(x=reorder(region,NSEG), y=NSEG, fill=region))+
  geom_boxplot(width=0.7, lwd=1)+
  #geom_jitter(width=0.05, colour="gray50", size=2)+
  geom_point(shape = 21, colour = "black", fill = "gray",size=1.5)+
  theme_classic()+
  ylab("Number of ROHs")+
  xlab("Region")+
  theme(legend.position="none", text=element_text(family="serif"))+
  scale_fill_manual(values = c("#ABD9E9", "#4575B4","#A50026"))
bar_NSEG

# adjust site labels IC -> JM, and ARLBNF -> WNA)

labels <- c("JM", "WNA", "SS") 
bar_KB <- ggplot(set, aes(x=reorder(region,MB), y=MB, fill=region))+
  geom_violin(width=0.8, lwd=1)+
  geom_boxplot(width=.1, fill="white", alpha=1)+
  #geom_point(shape = 21, colour = "black", fill = "gray",size=1.5)+
  theme_classic()+
  ylab("Total length of ROHs (Mb)")+
  xlab("Region")+
  theme(legend.position="none", text=element_text(family="serif"))+
  scale_fill_manual(values = c("#ABD9E9", "#4575B4","#A50026"))+
  scale_x_discrete(labels=labels)+
  scale_y_continuous(breaks = round(seq(min(set$KB), max(set$KB), by = 10),1))

  
bar_KB
ggsave("roh_plot_50kb_colour_violin_WNA_JM.png", width=4.5, height=3, dpi=600)

#library(patchwork)
#roh + bar_NSEG / bar_KB
#ggsave("roh_plots_100kb_50snp_5miss_colour.png", width=8, height=5, dpi=300)

#bar_NSEG + bar_KB
#ggsave("roh_plots_100kb_50snp_5miss_colour_barsonly.png", width=7, height=3.5, dpi=400)


