#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=merge_bam
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=def-coling
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=6
#SBATCH --output=%x-%j.out
#SBATCH --kill-on-invalid=yes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

module load StdEnv/2020 gcc/9.3.0 samtools/1.13

cd /scratch/edegreef/NBW/aligned_bams/merged_bams

# Merging these two pairs of sorted.bam files (trying to increase coverage for same individual, which were confirmed through kinship analysis on previous snp sets):
#O-BKW1333-717-502 & Ha14-06-716-502
#O-BKW1337-711-502 & Ha14-10-711-502

# Merge bams
samtools merge Ha14-06-716-502-2.sorted.bam /scratch/edegreef/NBW/aligned_bams/O-BKW1333-717-502.sorted.bam /scratch/edegreef/NBW/aligned_bams/intermediate_bams/Ha14-06-716-502.sorted.bam

samtools merge Ha14-10-711-502-2.sorted.bam /scratch/edegreef/NBW/aligned_bams/O-BKW1337-711-502.sorted.bam /scratch/edegreef/NBW/aligned_bams/intermediate_bams/Ha14-10-711-502.sorted.bam

# Index bams
samtools index Ha14-06-716-502-2.sorted.bam
samtools index Ha14-10-711-502-2.sorted.bam

# Check header 
#samtools view -H Ha14-06-716-502-2.sorted.bam > Ha14-06-716-502-2.header.sam
#samtools view -H Ha14-10-711-502-2.sorted.bam > Ha14-10-711-502-2.header.sam

# Deduplicate
module load nixpkgs/16.09 picard/2.20.6

# Remove duplicate reads
ls *.sorted.bam | sed 's/.sorted.bam$//' | parallel --jobs 4 'java -Xmx28000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.sorted.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true'

# Add read group information to each deDup bam
ls *deDup.bam | sed 's/.deDup.bam$//' | parallel --jobs 4 'java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={}.deDup.bam O={}.deDupRG.bam RGID={} RGPL=illumina RGSM={} RGLB=lib1 RGPU=unit1'

# Re-load samtools
module load StdEnv/2020 gcc/9.3.0 samtools/1.13

# Index bams
ls *deDupRG.bam | sed 's/.deDupRG.bam$//' | parallel --jobs 4 'samtools index {}.deDupRG.bam'

# Check coverage
./bam_coverage.sh Ha14-06-716-502-2.deDupRG.bam
./bam_coverage.sh Ha14-10-711-502-2.deDupRG.bam