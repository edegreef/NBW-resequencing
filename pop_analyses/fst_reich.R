# Pairwise FST using Reich's unbiased estimate to handle unequal and small sample sizes. (from Reich et al. 2009: https://www.nature.com/articles/nature08365)

# 1) Function code to create Reich.Fst from Gaetano et al. 2014: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091237#s6, supplement: https://doi.org/10.1371/journal.pone.0091237.s002)

Reich.Fst <- function(pop1,pop2,call.rate = 0.95, top.number = 10){
  # remove the SNPs that are not in common between the 2 populations
  snp.to.keep <- intersect(row.names(pop1),row.names(pop2))
  if (length(snp.to.keep) == 0){print("Error: no SNP in common");return(NULL)}
  pop1.k <- pop1[snp.to.keep,]
  pop2.k <- pop2[snp.to.keep,]
  # change the reference allele if is not concordant between the 2 populations
  if (sum(pop1.k$A1 == pop2.k$A1) != length(snp.to.keep)){
    idx <- which(pop1.k$A1 != pop2.k$A1)
    idx.rev <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 == pop2.k$A2)
    idx.rm  <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 != pop2.k$A2)
    if(length(idx.rev) > 0){
      provv <- pop1.k$A1[idx.rev]
      pop1.k$A1[idx.rev] <- pop1.k$A2[idx.rev]
      pop1.k$A2[idx.rev] <- provv
      provv <- pop1.k$x0[idx.rev]
      pop1.k$x0[idx.rev] <- pop1.k$x2[idx.rev]
      pop1.k$x2[idx.rev] <- provv}
    if(length(idx.rm) > 0){      
      pop1.k <- pop1.k[-idx.rm,]
      pop2.k <- pop2.k[-idx.rm,]}}
  # remove SNPs with low call rate in one or both populations
  x0 <- pop1.k$x0
  x1 <- pop1.k$x1
  x2 <- pop1.k$x2
  s <- x0 + x1 + x2
  y0 <- pop2.k$x0
  y1 <- pop2.k$x1
  y2 <- pop2.k$x2
  t <- y0 + y1 + y2
  idx.rm.pop1 <- which(s < max(s)*call.rate)
  idx.rm.pop2 <- which(t < max(t)*call.rate)
  idx.rm.all <- union(idx.rm.pop1,idx.rm.pop2)
  x0 <- x0[-idx.rm.all]
  x1 <- x1[-idx.rm.all]
  x2 <- x2[-idx.rm.all]
  s  <- s[-idx.rm.all]
  y0 <- y0[-idx.rm.all]
  y1 <- y1[-idx.rm.all]
  y2 <- y2[-idx.rm.all]
  t  <- t[-idx.rm.all]
  # compute SNP_Fst and global Fst estimators in presence of inbreeding   
  e.x <- ((x1 + 2*x2)/(2*s) - (y1 + 2*y2)/(2*t))^2 + x1/(4*s*s) + y1/(4*t*t)
  e.h1 <- (x0*x2 + (x0 + x2)*x1/2 + x1*(x1-1)/4)/(s*(s-1))
  e.h2 <- (y0*y2 + (y0 + y2)*y1/2 + y1*(y1-1)/4)/(t*(t-1))
  N <- e.x - e.h1/s - e.h2/t
  D <- N + e.h1 + e.h2
  Fst.v <- N/D
  names(Fst.v) <- row.names(pop1.k[-idx.rm.all,])
  Fst.o <- Fst.v[order(Fst.v,decreasing=TRUE)]
  F.global <- sum(N)/sum(D)
  se1 <- sd(N)/sqrt(length(N))
  se2 <- sd(D)/sqrt(length(N))
  se.F <- sqrt(se1*se1 + se2*se2)
  F_L95 <- F.global - 1.96*se.F
  F_U95 <- F.global + 1.96*se.F
  Z <- F.global/se.F
  p <- 2*(1 - pnorm(Z))
  if(p < 2e-16){p <- "less than 2e-16"}
  output <- list()
  output[[1]] <- c(F.global,F_L95,F_U95,p)
  names(output[[1]]) <- c("Reich.Fst","L.95%.CI","U.95%.CI","p.val")
  output[[2]] <- data.frame(Fst.o[1:top.number])
  names(output[[2]]) <- c("Reich.Fst")
  return(output)}


# 2) Set up data to run Reich's fst
# I couldn't get it to work as a combined genind object, so I did each pair separately. Formatting data into pop matrices with ncessary input info.

# From the supplemental page in Gaetano et al. 2014:
# "input data frame pop1 is a N x 5 matrix
# where N is the number of SNPs
# row names correspond to the SNP name
# x0 represent the number of samples with 0 copies of the variant allele
# x1 represent the number of samples with 1 copy of the variant allele
# x2 represent the number of samples with 2 copies of the variant allele
# A1 common allele
# A2 variant allele"


library(adegenet)
library(vcfR)
library(tidyverse)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/fst")

# read in vcfs for each pop 
vcf_arctic <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.arctic.vcf.gz") 
vcf_iceland <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.iceland.vcf.gz") 
vcf_labrador <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.labrador.vcf.gz") 
vcf_newfoundland <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.newfoundland.vcf.gz") 
vcf_scotianshelf <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.scotianshelf.vcf.gz") 

# convert vcf to genlight object
genlight_arctic <- vcfR2genlight(vcf_arctic)
genlight_iceland <- vcfR2genlight(vcf_iceland)
genlight_labrador <- vcfR2genlight(vcf_labrador)
genlight_newfoundland <- vcfR2genlight(vcf_newfoundland)
genlight_scotianshelf <- vcfR2genlight(vcf_scotianshelf)

# reformat data to matrices
mat_arctic <- as.matrix(genlight_arctic)
mat_iceland <- as.matrix(genlight_iceland)
mat_labrador <- as.matrix(genlight_labrador)
mat_newfoundland <- as.matrix(genlight_newfoundland)
mat_scotianshelf <- as.matrix(genlight_scotianshelf)

# transpose matrix so rows=snps, and also convert to dataframe
mat_arctic <- as.data.frame(t(mat_arctic))
mat_iceland <- as.data.frame(t(mat_iceland))
mat_labrador <- as.data.frame(t(mat_labrador))
mat_newfoundland <- as.data.frame(t(mat_newfoundland))
mat_scotianshelf <- as.data.frame(t(mat_scotianshelf))

# make columns of allele counts. The n in 'rowSums(mat[,1:n]...)' is sample size. (i.e. arctic has 11 samples). If the rowSums is not capped to sample size then it may include x0 - x2 columns in counts which will mess up counts.
# do for each pop
mat_arctic$x0 <- rowSums(mat_arctic[,1:11] == 0, na.rm=TRUE)
mat_arctic$x1 <- rowSums(mat_arctic[,1:11] == 1, na.rm=TRUE)
mat_arctic$x2 <- rowSums(mat_arctic[,1:11] == 2, na.rm=TRUE)

mat_iceland$x0 <- rowSums(mat_iceland[,1:7] == 0, na.rm=TRUE)
mat_iceland$x1 <- rowSums(mat_iceland[,1:7] == 1, na.rm=TRUE)
mat_iceland$x2 <- rowSums(mat_iceland[,1:7] == 2, na.rm=TRUE)

mat_labrador$x0 <- rowSums(mat_labrador[,1:3] == 0, na.rm=TRUE)
mat_labrador$x1 <- rowSums(mat_labrador[,1:3] == 1, na.rm=TRUE)
mat_labrador$x2 <- rowSums(mat_labrador[,1:3] == 2, na.rm=TRUE)

mat_newfoundland$x0 <- rowSums(mat_newfoundland[,1:6] == 0, na.rm=TRUE)
mat_newfoundland$x1 <- rowSums(mat_newfoundland[,1:6] == 1, na.rm=TRUE)
mat_newfoundland$x2 <- rowSums(mat_newfoundland[,1:6] == 2, na.rm=TRUE)

mat_scotianshelf$x0 <- rowSums(mat_scotianshelf[,1:9] == 0, na.rm=TRUE)
mat_scotianshelf$x1 <- rowSums(mat_scotianshelf[,1:9] == 1, na.rm=TRUE)
mat_scotianshelf$x2 <- rowSums(mat_scotianshelf[,1:9] == 2, na.rm=TRUE)

# subset matrix to just get snp name and allele counts
allele_count_arctic <- mat_arctic %>% select(c(x0,x1,x2))
allele_count_iceland <- mat_iceland %>% select(c(x0,x1,x2))
allele_count_labrador <- mat_labrador %>% select(c(x0,x1,x2))
allele_count_newfoundland <- mat_newfoundland %>% select(c(x0,x1,x2))
allele_count_scotianshelf <- mat_scotianshelf %>% select(c(x0,x1,x2))

# get the allele infos (A/T/G/C), read vcfs as tables
table_arctic <- read.table("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.arctic.vcf.gz")
table_iceland <- read.table("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.iceland.vcf.gz")
table_labrador <- read.table("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.labrador.vcf.gz")
table_newfoundland <- read.table("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.newfoundland.vcf.gz")
table_scotianshelf <- read.table("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.scotianshelf.vcf.gz")

# extracting the common allele (A1) and variant allele (A2) 
variants_arctic <- table_arctic %>% select(c(V4,V5))
colnames(variants_arctic) <- c("A1", "A2")

variants_iceland <- table_iceland %>% select(c(V4,V5))
colnames(table_iceland) <- c("A1", "A2")

variants_labrador <- table_labrador %>% select(c(V4,V5))
colnames(variants_labrador) <- c("A1", "A2")

variants_newfoundland <- table_newfoundland %>% select(c(V4,V5))
colnames(variants_newfoundland) <- c("A1", "A2")

variants_scotianshelf <- table_scotianshelf %>% select(c(V4,V5))
colnames(variants_scotianshelf) <- c("A1", "A2")


# merge allele counts and types, then it should be in the right format for Reich.Fst 
pop_arctic <- cbind(allele_count_arctic, variants_arctic)
pop_iceland <- cbind(allele_count_iceland, variants_iceland)
pop_labrador <- cbind(allele_count_labrador, variants_labrador)
pop_newfoundland <- cbind(allele_count_newfoundland, variants_newfoundland)
pop_scotianshelf <- cbind(allele_count_scotianshelf, variants_scotianshelf)

# run Reich's fst for each pair
fst_AR_IC <- Reich.Fst(pop_arctic, pop_iceland, call.rate=0.75, top.number=10)
fst_AR_LB <- Reich.Fst(pop_arctic, pop_labrador, call.rate=0.75, top.number=10)
fst_AR_NF <- Reich.Fst(pop_arctic, pop_newfoundland, call.rate=0.75, top.number=10)
fst_AR_SS <- Reich.Fst(pop_arctic, pop_scotianshelf, call.rate=0.75, top.number=10)
fst_IC_LB <- Reich.Fst(pop_iceland, pop_labrador, call.rate=0.75, top.number=10)
fst_IC_NF <- Reich.Fst(pop_iceland, pop_newfoundland, call.rate=0.75, top.number=10)
fst_IC_SS <- Reich.Fst(pop_iceland, pop_scotianshelf, call.rate=0.75, top.number=10)
fst_LB_NF <- Reich.Fst(pop_labrador, pop_newfoundland, call.rate=0.75, top.number=10)
fst_LB_SS <- Reich.Fst(pop_labrador, pop_scotianshelf, call.rate=0.75, top.number=10)
fst_NF_SS <- Reich.Fst(pop_newfoundland, pop_scotianshelf, call.rate=0.75, top.number=10)

# take a look at results
fst_AR_IC #0.0044668074487502
fst_AR_LB #0.00348639995077524
fst_AR_NF #0.00348639995077524
fst_AR_SS #0.00795411534952307
fst_IC_LB #0.00978923539254032
fst_IC_NF #0.00828351161628797
fst_IC_SS #0.0104239348934766
fst_LB_NF #0.00353676945870719
fst_LB_SS #0.0100488019039871
fst_NF_SS #0.00717667192101823

# save results in a thing
a <- data.frame(matrix(unlist(fst_AR_IC[1]), nrow=length(fst_AR_IC[1]), byrow=TRUE))
b <- data.frame(matrix(unlist(fst_AR_LB[1]), nrow=length(fst_AR_LB[1]), byrow=TRUE))
c <- data.frame(matrix(unlist(fst_AR_NF[1]), nrow=length(fst_AR_NF[1]), byrow=TRUE))
d <- data.frame(matrix(unlist(fst_AR_SS[1]), nrow=length(fst_AR_SS[1]), byrow=TRUE))
e <- data.frame(matrix(unlist(fst_IC_LB[1]), nrow=length(fst_IC_LB[1]), byrow=TRUE))
f <- data.frame(matrix(unlist(fst_IC_NF[1]), nrow=length(fst_IC_NF[1]), byrow=TRUE))
g <- data.frame(matrix(unlist(fst_IC_SS[1]), nrow=length(fst_IC_SS[1]), byrow=TRUE))
h <- data.frame(matrix(unlist(fst_LB_NF[1]), nrow=length(fst_LB_NF[1]), byrow=TRUE))
i <- data.frame(matrix(unlist(fst_LB_SS[1]), nrow=length(fst_LB_SS[1]), byrow=TRUE))
j <- data.frame(matrix(unlist(fst_NF_SS[1]), nrow=length(fst_NF_SS[1]), byrow=TRUE))

merged <- rbind(a,b,c,d,e,f,g,h,i,j)

# making column with pair labels to add to merged data (make sure in same order)
pair <- data.frame(c("AR-IC", "AR-LB", "AR-NF", "AR-SS", "IC-LB", "IC-NF", "IC-SS", "LB-NF", "LB-SS", "NF-SS"))
fst_together <- cbind(pair, merged)
colnames(fst_together) <- c("pair", "Reich_fst", "ll_95CI", "ul_95CI", "p_val")

write.csv(fst_together, "Reichs_fst_output_NBW2_SNPS_2M.csv")
