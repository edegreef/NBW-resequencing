#!/bin/bash

# run admixture

cd /home/degreefe/NBW/reseq_newsnps/samplesizes/admixture

# prep bim file to have integer contig names
sed 's/^.\{,6\}//' NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.bim > NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.bim

# keep bfiles names consistent
cp NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.bed NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.bed
cp NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.fam NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.fam

# run admixture for K 1-5 (-B includes bootstrapping to calculate standard errors)
for K in 1 2 3 4 5
do /home/degreefe/programs/admixture/dist/admixture_linux-1.3.0/admixture -B NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.bed $K
done

# calculate CV error
for K in 1 2 3 4 5
do /home/degreefe/programs/admixture/dist/admixture_linux-1.3.0/admixture --cv NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.ec.bed $K | tee log${K}.out
done