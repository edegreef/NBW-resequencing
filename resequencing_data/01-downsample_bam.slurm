#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=downsample_batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=24
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

# downsampling some of the bam files to a lower coverage, aiming for ~5x for the higher coverage samples, running script in batches with 4 samples each (20 total). Proportion to 'keep' value based on target coverage divided by original coverage

bam_ID=NBW-2016-14
bam_ID2=O-BKW1305-710-504
bam_ID3=O-BKW1312-711-504
bam_ID4=O-BKW1340-717-505
keep=0.833

cd /scratch/edegreef/NBW_downsampled_bams
module load nixpkgs/16.09 gatk/4.1.2.0

gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW_dedup/$bam_ID.deDupRG.bam -O $bam_ID.deDupRG.downsampled.bam -P $keep
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW_dedup/$bam_ID2.deDupRG.bam -O $bam_ID2.deDupRG.downsampled.bam -P $keep
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW_dedup/$bam_ID3.deDupRG.bam -O $bam_ID3.deDupRG.downsampled.bam -P $keep
gatk --java-options "-Xmx24G" DownsampleSam -I /scratch/edegreef/NBW_dedup/$bam_ID4.deDupRG.bam -O $bam_ID4.deDupRG.downsampled.bam -P $keep


# check new modal coverage
module load StdEnv/2020 samtools/1.11

./bam_coverage.sh $bam_ID.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID2.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID3.deDupRG.downsampled.bam
./bam_coverage.sh $bam_ID4.deDupRG.downsampled.bam

# index all new bams at the end (run this only on last batch)
module load StdEnv/2020 samtools/1.11

for i in *.bam
do
echo "Indexing: "$i        
samtools index $i
done
