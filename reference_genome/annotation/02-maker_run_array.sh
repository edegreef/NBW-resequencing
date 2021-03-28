#!/bin/bash
#SBATCH --job-name=maker_run_array
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50000
#SBATCH --cpus-per-task=20
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 
#SBATCH --array=1-21

cd /home/edegreef/maker

# load maker 2.31.10 version
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 maker/2.31.10

# set up array info/list to run all chunks in one job submission
list=NBW_fasta_20chunk_list.txt
string="sed -n "${SLURM_ARRAY_TASK_ID}"p ${list}"
sample=$($string)
echo ${sample}

# start maker things
echo "Starting run at: `date`"

# copy maker files into a TMPDIR and set to TMPDIR to save file space
cp maker*.ctl ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

# run maker on reference fasta chunk
maker -g /home/edegreef/maker/reference_chunks/${sample}.fasta maker_opts.ctl maker_bopts.ctl maker_exe.ctl -cpus ${SLURM_CPUS_PER_TASK}

# merge output fastas and gffs
fasta_merge -d ${sample}.maker.output/${sample}_master_datastore_index.log
gff3_merge -d ${sample}.maker.output/${sample}_master_datastore_index.log

# copy merged output files to specific directory since everything in TMPDIR will go away
cp ${sample}.maker.output/${sample}_master_datastore_index.log ${SLURM_SUBMIT_DIR}
cp *.all.maker.*.fasta ${SLURM_SUBMIT_DIR}
cp *.all.gff ${SLURM_SUBMIT_DIR}

echo "Program finished with exit code $? at: `date`"