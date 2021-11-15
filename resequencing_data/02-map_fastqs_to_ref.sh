#!/bin/bash

#SBATCH --time=168:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=bwa_mem
#SBATCH --output=%x-%j.out

# Load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.12

cd /scratch/edegreef/NBW

# Index reference (took about an hour for 2GB genome)
#bwa index /scratch/edegreef/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta
#samtools faidx /scratch/edegreef/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta

# Map paired.fastq's to reference. This step took a long time, example: ~5 days to make one 19GB bam file. Depending on how many samples, may want to increase job time limit or cpus (or spread data over several computing clusters and run multiple jobs to avoid using up all resource allocation in one place).
ls *_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz$//' | parallel --jobs $SLURM_CPUS_PER_TASK 'bwa mem /scratch/edegreef/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz | samtools sort -T {}_tmp -o {}.sorted.bam'

# Index bam file
ls *.sorted.bam | parallel 'samtools index {}'
