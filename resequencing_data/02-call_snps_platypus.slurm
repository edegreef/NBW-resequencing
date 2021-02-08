#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=platypus
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=24
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

module load nixpkgs/16.09  gcc/5.4.0
module load nixpkgs/16.09  gcc/7.3.0
module load platypus/0.8.1

# run platypus to call variants from bams & reference genome 
# make sure all bam files are indexed beforehand. ALso couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=/scratch/edegreef/NBW_downsampled_bams/HamLB03-1.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/HamLB03-2.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/HamLB03-3.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/Hyam-2016-06-pooled.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-12-710-509.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-13-702-507.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-15-708-507.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-16-705-509.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2016-02-706-509.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-2016-05-715-506.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2016-06.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2016-09.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-2016-10-706-505.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2016-11.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-2016-12.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2016-13-707-505.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2016-14.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-2016-17.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2017-01.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/NBW-2017-03.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/NBW-2017-06.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2017-2-705-505.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2017-4-704-505.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2019-01-709-505.deDupRG.bam,/scratch/edegreef/NBW_dedup/NBW-2019-02-714-509.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW0957-712-505.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW0958-716-503.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/O-BKW0960-717-503.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW0962-718-503.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1305-710-504.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1312-711-504.deDupRG.downsampled.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1313-718-505.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/O-BKW1314-712-504.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1315-716-504.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1316-717-504.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1317-710-503.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1318-711-503.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1319-712-503.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/O-BKW1328-711-505.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1330-716-505.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/O-BKW1331-712-502.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1332-716-502.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1333-717-502.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1334-718-502.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1335-710-505.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1337-711-502.deDupRG.bam,/scratch/edegreef/NBW_downsampled_bams/O-BKW1340-717-505.deDupRG.downsampled.bam,/scratch/edegreef/NBW_dedup/O-BKW1343-710-502.deDupRG.bam,/scratch/edegreef/NBW_dedup/O-BKW1345-718-504.deDupRG.bam \
--refFile=/scratch/edegreef/Northern_bottlenose_whale_051018_shortLabel.fasta \
--output=/scratch/edegreef/platypus/NBW_platypus_allvariantcalls.vcf \
--nCPU=24
