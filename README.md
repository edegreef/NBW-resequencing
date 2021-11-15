![Logo](https://github.com/edegreef/NBW-resequencing/blob/main/NBW-cartoon-forgit.JPG)

This is a repository for scripts used in analyzing northern bottlenose whale (*Hyperoodon ampullatus*) resequencing data. These scripts were executed on computing resources through the University of Manitoba and Compute Canada. Project is a collaboration with the Marine Gene Probe Lab & Whitehead Lab at the Dalhousie University and the University of Manitoba.
<br/>
### Reference genome folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/reference_genome)
* Evaluated genome quality with *QUAST*, *Assemblathon*, *BUSCO*
* Mapped chromosome alignment through synteny with *Satsuma*
* Annotated genome with *MAKER* and added functional annotations with *BLAST+* and *Interproscan*

### Resequencing data folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/resequencing_data)
### Processing sequencing data
01. Trimmed fastqs with *trimmomatic*
02. Mapped fastqs to reference genome with *bwa*
03. Removed duplicate reads with *picard*
04. Added readgroups with *picard*
05. Checked modal coverage of all the samples
06. Downsampled some bams with *GATK* to maintain a consistent modal coverage before calling snps
07. Merged sequence files for two individuals to increase coverage
#### Initial SNP prep
08. Called variants from bams and reference genome using *Platypus*
09. Removed indels then extracted snp metrics with *bcftools* and *vcftools*
10. Looked at snp metrics in *R*
11. Filtered snps based on quality with *GATK* and *vcftools*
12. Used *bcftools* to edit sample ID label (removing the path in each sample ID name) and to add word "contig" to each scaffold name so data can be used in downstream analyses that won't take data with an interger as scaffold name.
#### Sex chromosomes
13. Notes on using *DifCover* for coverage comparisons between males and females
14. Examined *DifCover* results in R
15. Used *bedtools* to make annotated bed file with X and Y regions 
16. Finalized list of sex scaffolds in R
#### More SNP prep
17. Filtered X and Y-linked SNPs. Used *bcftools* to create vcfs of X and Y snps, then used CHROM and POS information from those vcfs to filter out X and Y snps in full snp set using *vcftools*. Also made a separate file with LD-pruned snps using `LD_pruning.sh`.
18. Note on running *breakdancer* to detect structural variants.
19. Examining structural variant regions in R and creating a list from vcf positions for filtering.
20. a) Filtered out some individuals with very high missing data and ran pairise kinship estimates using *plink*. b) Examined kinship results in *R* to identify kin pairs. c) Removed an individual from each kin pair, and prepped snp files for pop analyses. 

### Population analyses folder: [:file_folder:](https://github.com/edegreef/NBW-resequencing/tree/main/pop_analyses)
* Made site map including ocean depth using a variety of *R* packages
* Evaluated PCA with *pcadapt* and created heatmap matrix of PC1 and PC2 distances
* Evaluated t-SNE with *Rtsne*
* Estimated FST and heterozygosity with *hierfstat*, created FST heatmap, and also analyzed isolation-by-distance
* Identified private and shared snps with bcftools isec
* Ran sNMF with *LEA*, looked at results, and plotted pie chart admixture map
* Ran *ADMIXTURE* with shell script
* Filtered beagle file to match snp sites and then ran *NGSadmix*
* Compared ancestral admixture analyses together in *R*
* Estimated changes in effective population sizes with *SMC++* and *PSMC*

