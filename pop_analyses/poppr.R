# creating phylogenetic tree with genetic distances using poppr package

library(vcfR)
library(mmod)
library(poppr)
library(ape)
library(magrittr)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/poppr_tree")

vcf <- read.vcfR("C:/Users/eveli/Dropbox/NBW-me/reseq_newsnps/fst/NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.vcf")

# Reformat vcfR object into a genlight object
#snps <- vcfR2genlight(vcf)

# Reformat to a matrix
#mat <- as.matrix(snps)

#mat[mat==0] <- 11
#mat[mat==1] <- 12
#mat[mat==2] <- 22

# Provesti distance
#PD <- provesti.dist(mat)

#theTree <- PD %>%
#  nj() %>%    # calculate neighbor-joining tree
#  ladderize() # organize branches by clade
#plot(theTree)
#add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

#boot <- aboot(mat, dist = provesti.dist, sample = 100, tree = "nj", cutoff = 50, quiet = TRUE)

# using genind object instead of matrix
my_genind <- vcfR2genind(vcf)
my_genind

# add pop info for now
sample_info <- read.csv("C:/Users/eveli/Dropbox/NBW-me/sample_info_updated_25Jan2021.csv", header=T)
sample_info_subset <- subset(sample_info, remove_indiv!="Y") 

# extract pop info
my_genind@pop <- as.factor(sample_info_subset$region)
my_genind

genind_PD <- provesti.dist(my_genind)

theTree <- genind_PD %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) 

genind_boot <- aboot(my_genind, dist = provesti.dist, sample = 100, tree = "nj", cutoff = 50, quiet = TRUE)

# should save RData at this point b/c the aboot can take a while and don't want to redo later

cols <- c("#ABD9E9", "#4575B4", "#FEE090", "#F46D43", "#A50026")
#plot.phylo(genind_boot)

#### phylogram
png("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/poppr_tree/genind_obj/boot_obj_colorbylat.png", width=2300, height=3000, res=400)
phylo <- plot.phylo(genind_boot, type="phylogram", 
           cex = 0.8, font = 2, adj = 0,
           tip.color = cols[my_genind@pop],
           label.offset = 0.0125)
          # x.lim=0.02)
nodelabels(genind_boot$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)
dev.off()


#### radial
png("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/poppr_tree/genind_obj/boot_radial_colorbylat2.png", width=4000, height=5000, res=400)
ax <- plot.phylo(genind_boot, type="radial", 
                    cex = 0.8, font = 2, adj = 0,
                    tip.color = cols[my_genind@pop],
                    label.offset = 0.0125)
                    #x.lim=0.02)
nodelabels(genind_boot$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
dev.off()

