#!/bin/bash

# prepping bed files to merge info of all 6 DifCover runs into one file to examine X and Y-linked scaffolds

# make bed file from genome fasta or fasta.fai 
# using .fasta
faidx -i bed Northern_bottlenose_whale_051018_shortLabel.fasta > NBW_genome.bed

# using .fasta.fai
#awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' Northern_bottlenose_whale_051018_shortLabel.fasta.fai > NBW_genome.bed

# make 10kb windows
bedtools makewindows -b NBW_genome.bed -w 10000 > NBW_genome.10KBwindows.bed

# next need to annotate bed file with the M/F series run outputs for which scaffolds are sex-linked or not. Order of files here are M1F1, M1F2, M2F1, M2F2, M1M2, F1F2

# will do X-linked first
bedtools annotate -i NBW_genome.10KBwindows.bed -files X_linked_scafwindows.M1F1.bed X_linked_scafwindows.M1F2.bed X_linked_scafwindows.M2F1.bed X_linked_scafwindows.M2F2.bed Control_scafwindows.M1M2.bed Control_scafwindows.F1F2.bed > NBW_genome.10KBwindows.Xlinked.bed

# then Y-linked (control stays the same)
bedtools annotate -i NBW_genome.10KBwindows.bed -files Y_linked_scafwindows.M1F1.bed Y_linked_scafwindows.M1F2.bed Y_linked_scafwindows.M2F1.bed Y_linked_scafwindows.M2F2.bed Control_scafwindows.M1M2.bed Control_scafwindows.F1F2.bed > NBW_genome.10KBwindows.Ylinked.bed

# then load the new bed files into R to sort out the stuffs to make list of sex-linked scaffolds, selecting scaffolds that are consistently enriched in M/F combos, and not enriched in the MM FF controls 