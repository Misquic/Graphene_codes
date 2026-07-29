#!/bin/bash -l
set -u

module load Miniconda3
eval "$(conda shell.bash hook)"
conda activate normal

echo $#
if [ $# -ne 1 ] && [ $# -ne 2 ] && [ $# -ne 3 ]; then
  echo "illegal number of parameters"
  exit 1
fi

SF=8
SEED=12345
WAIT=""

if [ $# -ge 2 ]; then
  if [ $2 == "1" ] || [ $2 == "T" ] || [ $2 == "t" ] || [ $2 == "true" ] || [ $2 == "True" ]; then
    WAIT="--wait"
  fi
fi

if [ $# -eq 3 ]; then
  SEED=$3
fi

EXECUTABLE=Transport2D_seed
DIR="$1"_sf"$SF"_S"$SEED/"


echo "preparing commands and directories for $DIR"

pythonCmd="python $PLOT_SCRIPTS/wholeSim.py 1 8 -60 40 \
--dB=0.1   \
--dVb=1 \
--allResultsDir=$DIR \
--clearDir=1 \
--prepCmdsOnly=1 \
--saveStdout=1 \
--saveCurrents=1 \
--sf=$SF \
--seed=$SEED \
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

RUNS_PER_JOB=30
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
-J Bilayer_\"$EXECUTABLE\"_sf\"$SF\"_$SEED \
$WAIT \
$PLOT_SCRIPTS/scriptHPC.sh \"$DIR\" $RUNS_PER_JOB"

# Log it
echo "$sbatchCmd" >> lastSbatch.txt

echo "Post process cmd: $PLOT_SCRIPTS/wholeSim.py --plotAll=1 --allResultsDir=\"$DIR\" --leadInfo=\"3 4 1 2\""

# Execute it
eval "$sbatchCmd"

if [ $WAIT == "--wait" ]; then
  echo "$PLOT_SCRIPTS/wholeSim.py --plotAll=1 --allResultsDir=\"$DIR\" --leadInfo=\"3 4 1 2\""
  python $PLOT_SCRIPTS/wholeSim.py --plotAll=1 --allResultsDir="$DIR" --leadInfo="3 4 1 2"
fi

### python wholeSim.py --plotAll=1 --allResultsDir="$DIR"
