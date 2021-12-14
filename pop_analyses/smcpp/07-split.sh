#!/bin/bash

#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-20
#SBATCH --job-name=IC_ARLBNF_split

source /home/edegreef/smcpp/bin/activate

pop=IC_SS

echo run_${SLURM_ARRAY_TASK_ID}

smc++ split -o /scratch/edegreef/smcpp/split_estimate/$pop/combined_run${SLURM_ARRAY_TASK_ID} /scratch/edegreef/smcpp/smc_out_timepoints/IC/run_${SLURM_ARRAY_TASK_ID}/model.final.json /scratch/edegreef/smcpp/smc_out_timepoints/SS/run_${SLURM_ARRAY_TASK_ID}/model.final.json /scratch/edegreef/smcpp/masked_data_run_ICSS/*.smc
