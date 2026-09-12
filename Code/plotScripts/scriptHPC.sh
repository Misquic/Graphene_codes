#!/bin/bash -l
#SBATCH -J Bilayer_ni
#SBATCH --cpus-per-task=2
#SBATCH --mem 7000
#SBATCH --time=2:30:00
#SBATCH --account=plgmrenca2025-cpu
#SBATCH -p plgrid

DIR=$1
RUNS_PER_JOB=$2
TOTAL=$(wc -l < "$DIR/commands.txt")

module load intel/2023b
module load Miniconda3
eval "$(conda shell.bash hook)"
conda activate normal
cd "$SLURM_SUBMIT_DIR"

export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

START=$(( (SLURM_ARRAY_TASK_ID-1)*RUNS_PER_JOB + 1 ))
END=$(( START + RUNS_PER_JOB - 1 ))

if [ $END -gt $TOTAL ]; then
  END=$TOTAL
fi

sed -n "${START},${END}p" "${DIR}/commands.txt" |
while IFS= read -r CMD
do
  eval "$CMD"
done
