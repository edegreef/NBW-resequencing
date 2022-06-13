# Making list of SNPs to remove that are not in Hardy-Weinberg equilibrium

# First look at hardy statistics through vcftools with --hardy
#vcftools --gzvcf NBW_SNPS_filtered_population_analyses.vcf.gz --hardy --out NBW_SNPS_filtered_population_analyses_hwe
#vcftools --gzvcf NBW_SNPS_filtered_demography_analyses.vcf.gz --hardy --out NBW_SNPS_filtered_demography_analyses_hwe
 

# In R, look at the .hwe output
data<-read.table("NBW_SNPS_filtered_demography_analyses_hwe.hwe",header=T)

split <- data.frame(do.call('rbind', strsplit(as.character(data$OBS.HOM1.HET.HOM2.),'/',fixed=TRUE)))

part <-data[,1:2]
ready <-cbind(part,split)

names(ready)<-c("chr","pos","homo1","het","homo2")
ready$homo1<-as.numeric(as.character(ready$homo1))
ready$homo2<-as.numeric(as.character(ready$homo2))
ready$het<-as.numeric(as.character(ready$het))
ready$obs_het<-0

for (i in nrow(ready)){
  ready$obs_het<-ready$het/(ready$homo1+ready$het+ready$homo2)
}

# look at histogram
hist(ready$obs_het,breaks=60)

na<-subset(ready,ready$obs_het=="NaN")
het<-subset(ready,ready$obs_het>=0.6)
remove<-rbind(na,het)
remove_2col<-remove[,1:2]

# save list of snps
write.table(remove_2col,file="tmp_HWE_filter_demoganalyses",sep="\t",quote=FALSE,row.names=F)


# back in vcftools, use --exclde positions (example:)
#vcftools --gzvcf NBW_SNPS_filtered_population_analyses.vcf.gz --exclude-positions tmp_HWE_filter_demoganalyses --recode --recode-INFO-all --out NBW_SNPS_filtered_demography_analyses.hwe

