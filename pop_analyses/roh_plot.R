# plotting roh using outputs from PLINK
# ran ROH for each population/group separately
# for one group:
# /home/degreefe/programs/plink --allow-extra-chr --bfile snp.scotian_shelf --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 50 --homozyg-snp 50 --out SS

# then plot in R!
library(ggplot2)

JM <- read.table("IC.hom.indiv", header=T)
WNA <- read.table("ARLBNF.hom.indiv", header=T)
SS <- read.table("SS.hom.indiv", header=T)

JM$region <- "JM"
WNA$region <- "WNA"
SS$region <- "SS"

merge <- rbind(JM, WNA, SS)

set <- merge[, c("NSEG", "KB", "region")]

roh <- ggplot(data=set, aes(x=KB,y=NSEG))+
  geom_point(aes(color=region), size=4, alpha=0.7)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (kb)")+
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

# adjust ARLBNF label to WNA
labels <- c("JM", "WNA", "SS")

bar_KB <- ggplot(set, aes(x=reorder(region,KB), y=KB, fill=region))+
  geom_violin(width=0.8, lwd=1)+
  geom_boxplot(width=.1, fill="white", alpha=1)+
  #geom_point(shape = 21, colour = "black", fill = "gray",size=1.5)+
  theme_classic()+
  ylab("Total length of ROHs (kb)")+
  xlab("Region")+
  theme(legend.position="none", text=element_text(family="serif"))+
  scale_fill_manual(values = c("#ABD9E9", "#4575B4","#A50026"))+
  scale_x_discrete(labels=labels)+
  scale_y_continuous(breaks = round(seq(min(set$KB), max(set$KB), by = 10000),1))
  
bar_KB
