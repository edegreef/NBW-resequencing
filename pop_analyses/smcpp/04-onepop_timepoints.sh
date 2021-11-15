#!/bin/bash
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-20
#SBATCH --job-name=estimate_SS

# 20 array runs/iterations of 'smc++ estimate' for a population (this script is for one pop, run script for each pop separately)
# also remember to make the folders first ("smc_out_timepoints" a folder for each pop ID)

pop=SS

source /home/edegreef/smcpp/bin/activate

mkdir /scratch/edegreef/SMC_min100kb/smc_out_timepoints/$pop/run_${SLURM_ARRAY_TASK_ID}/

echo smc_out_timepoints/$pop/run_${SLURM_ARRAY_TASK_ID}/

smc++ estimate 1.53e-8 -o /scratch/edegreef/SMC_min100kb/smc_out_timepoints/$pop/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/SMC_min100kb/masked_data2/$pop.*
