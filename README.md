![Logo](https://github.com/edegreef/NBW-resequencing/blob/main/NBW-cartoon-forgit.JPG)

This is a repository for scripts I used in analyzing northern bottlenose whale (*Hyperoodon ampullatus*) resequencing data. These scripts were executed on computing resources through the University of Manitoba and Compute Canada. Project is a collaboration with the Marine Gene Probe Lab & Whitehead Lab at the Dalhousie University and the University of Manitoba.
<br/>
### Reference genome folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/reference_genome)
* Evaluated genome quality with *QUAST*, *Assemblathon*, *BUSCO*)
* Mapped chromosome alignment (synteny with *Satsuma*)

### Resequencing data folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/resequencing_data)
#### Initial SNP prep
01. Lowered variation in sample coverage, and checked modal coverage with `bam_coverage.sh` (all samples were assess for modal coverage before downsampling, then checked again afterwards).
02. Called variants from bams and reference genome using *Platypus*
03. Removed indels then extracted snp metrics with *bcftools* and *vcftools*
04. Looked at snp metrics in *R*
05. Filtered snps based on quality with *GATK* and *vcftools*
06. Used *bcftools* to edit sample ID label (removing the path in each sample ID name) and to add word "contig" to each scaffold name (in hindsight I probably should have used "scaffold" instead of "contig") so data can be used in downstream analyses that won't take data with an interger as scaffold name.
#### Sex chromosomes
07. Notes on using *DifCover* for coverage comparisons between males and females
08. Examined *DifCover* results in R
09. Used *bedtools* to make annotated bed file with X and Y regions 
10. Finalized list of sex scaffolds in R
#### More SNP prep
11. Filtered X and Y-linked SNPs. Used *bcftools* to create vcfs of X and Y snps, then used CHROM and POS information from those vcfs to filter out X and Y snps in full snp set using *vcftools*. Also made a separate file with LD-pruned snps using `LD_pruning.sh`.
12. Filtered out some individuals with very high missing data and ran pairise kinship estimates using *plink*.
13. Examined kinship results in *R* to identify kin pairs.
14. Removed an individual from each kin pair, and prepped snp files for pop analyses. 

### Population analyses folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/pop_analyses)
* Made site map including ocean depth using a variety of *R* packages
* Evaluated PCA with *pcadapt* and created heatmap matrix of PC1 and PC2 distances
* Evaluated t-SNE with *Rtsne*
* Estimated FST with *hierfstat*, created heatmap, and also analyzed isolation-by-distance
* Looked at SNMF results and plotted pie chart admixture map
* Ran *ADMIXTURE* and *NGSadmix* with shell script
* Compared ancestral admixture analyses together in R
