#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=picard_readgroup
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/8.3.0 picard/2.20.6 samtools/1.9

# Create reference genome dictionary (.dict file)
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=/scratch/edegreef/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta O=/scratch/edegreef/ref_genome/Northern_bottlenose_whale_051018_shortLabel.dict 

# Add read group information to each deDup bam
ls /scratch/edegreef/NBW/aligned_bams/*deDup.bam | sed 's/.deDup.bam$//' | parallel --jobs 4 'java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={}.deDup.bam O={}.deDupRG.bam RGID={} RGPL=illumina RGSM={} RGLB=lib1 RGPU=unit1'

# Index bams
ls /scratch/edegreef/NBW/aligned_bams/*deDupRG.bam | sed 's/.deDupRG.bam$//' | parallel --jobs 4 'samtools index {}.deDupRG.bam'
