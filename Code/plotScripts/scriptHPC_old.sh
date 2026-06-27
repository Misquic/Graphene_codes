#!/bin/bash -l
#SBATCH -J Bilayer
#SBATCH --array=1-227%100
#SBATCH --cpus-per-task=3
#SBATCH --mem 8000
#SBATCH --time=70:59:00
#SBATCH --account=plgmrenca2025-cpu
#SBATCH -p plgrid
#SBATCH --output=./results/slurm_%A_%a.out
#SBATCH --error=./results/slurm_%A_%a.err

module add intel/2023b
module add Miniconda3
eval "$(conda shell.bash hook)"
conda activate normal
cd $SLURM_SUBMIT_DIR

export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

START=$(( (SLURM_ARRAY_TASK_ID-1)*100 + 1 ))
END=$(( START + 99 ))

TOTAL=$(wc -l < ./results/resultsTest4/commands.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

sed -n "${START},${END}p" ./results/resultsTest4/commands.txt | while read CMD
do
    eval "$CMD"
done
