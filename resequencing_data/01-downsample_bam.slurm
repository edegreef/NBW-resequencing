#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=24
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=downsample
#SBATCH --output=%x-%j.out

# Downsampling some of the bam files to a lower coverage, aiming for ~5x for the higher coverage samples
# Proportion to 'keep' value based on target coverage divided by original coverage

# In this case just need to do for 9 samples:
#6x	- NBW-2019-B17-716-503 (keep 0.833)
#6x	- NBW1904-718-503 (keep 0.833)
#7x	- Ha14-07-718-502 (keep 0.714)
#9x	- NBW19-03-712-504 (keep 0.556)
#10x - NBW-2019-B11-718-505 (keep 0.5)
#10x - NBW-2019-B13-711-505 (keep 0.5)
#10x - Ha14-13-716-505 (keep 0.5) ##for some reason the first time this over-downsampled to 2, so increased keep to 0.6 and re-ran that one individually
#11x - NBW-2019-B14-712-505 (keep 0.455)
#11x - Ha14-18-717-505 (keep 0.455)


bam_ID1=NBW-2019-B17-716-503
bam_ID2=NBW1904-718-503
bam_ID3=Ha14-07-718-502
bam_ID4=NBW19-03-712-504
bam_ID5=NBW-2019-B11-718-505
bam_ID6=NBW-2019-B13-711-505
bam_ID7=Ha14-13-716-505
bam_ID8=NBW-2019-B14-712-505
bam_ID9=Ha14-18-717-505

keep1=0.833
keep2=0.833
keep3=0.714
keep4=0.556
keep5=0.5
keep6=0.5
keep7=0.5
keep8=0.455
keep9=0.455

cd /scratch/edegreef/NBW/aligned_bams/dedupRG/downsampled_bams
module load nixpkgs/16.09 gatk/4.1.2.0

gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID1.deDupRG.bam -O $bam_ID1.deDupRG.downsampled.bam -P $keep1
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID2.deDupRG.bam -O $bam_ID2.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID3.deDupRG.bam -O $bam_ID3.deDupRG.downsampled.bam -P $keep3
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID4.deDupRG.bam -O $bam_ID4.deDupRG.downsampled.bam -P $keep4
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID5.deDupRG.bam -O $bam_ID5.deDupRG.downsampled.bam -P $keep5
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID6.deDupRG.bam -O $bam_ID6.deDupRG.downsampled.bam -P $keep6
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID7.deDupRG.bam -O $bam_ID7.deDupRG.downsampled.bam -P $keep7
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID8.deDupRG.bam -O $bam_ID8.deDupRG.downsampled.bam -P $keep8
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW/aligned_bams/dedupRG/$bam_ID9.deDupRG.bam -O $bam_ID9.deDupRG.downsampled.bam -P $keep9

# Check new modal coverage
module load StdEnv/2020 samtools/1.11

./bam_coverage.sh $bam_ID1.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID2.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID3.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID4.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID5.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID6.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID7.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID8.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID9.deDupRG.downsampled.bam

# Index all new bams at the end
module load StdEnv/2020 samtools/1.11

for i in *.bam
do
echo "Indexing: "$i        
samtools index $i
done
