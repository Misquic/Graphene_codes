#!/bin/bash -l
set -u

module load Miniconda3
eval "$(conda shell.bash hook)"
conda activate normal

echo $#
if [ $# -ne 1 ]; then
  echo "illegal number of parameters"
  exit 1
fi

SF=8
EXECUTABLE=Transport2DI_rand
DIR="$1"_sf"$SF"

echo "preparing commands and directories for $DIR"

pythonCmd="python $SCRIPTS/wholeSim.py 1 8 -60 40 \
--dB=0.1   \
--dVb=1 \
--allResultsDir=$DIR \
--clearDir=1 \
--prepCmdsOnly=1 \
--saveStdout=1 \
--saveCurrents=1 \
--sf=$SF \
--Executable=\"./$EXECUTABLE\""

echo "$pythonCmd" >> commandGen.txt
eval "$pythonCmd"
mv ./commandGen.txt "./$DIR/commandGen.txt"

mkdir -p "$DIR/outs"

echo "running supercomputer"

TOTAL=$(wc -l < "$DIR/commands.txt")
if [ "$TOTAL" -eq 0 ]; then
    echo "$DIR/commands.txt is empty"
    exit 1
fi

RUNS_PER_JOB=15
NJOBS=$(( (TOTAL + RUNS_PER_JOB - 1) / RUNS_PER_JOB ))
echo "NJOBS $NJOBS"
echo "RUNS_PER_JOB $RUNS_PER_JOB"
if [ "$NJOBS" -gt 20000 ]; then
  echo "Submitted too much jobs: $NJOBS"
  exit 1
fi

# Create the sbatch command
sbatchCmd="sbatch \
--array=1-${NJOBS} \
--output=\"${DIR}/outs/slurm_%A_%a.out\" \
--error=\"${DIR}/outs/slurm_%A_%a.err\" \
-J Bilayer_\"$EXECUTABLE\"_sf\"$SF\" \
$SCRIPTS/scriptHPC.sh \"$DIR\" $RUNS_PER_JOB"

# Log it
echo "$sbatchCmd" >> lastSbatch.txt

# Execute it
eval "$sbatchCmd"

### python wholeSim.py --plotAll=1 --allResultsDir="$DIR"