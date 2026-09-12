#!/bin/bash -l
set -u

echo $#
if [ $# -ne 1 ] && [ $# -ne 2 ]; then
  echo "illegal number of parameters"
  exit 1
fi

DIR=$1
mkdir $DIR
for SEED in 1 12 1234 12345 112 1123 11234 54321
do
  COMMAND="$PLOT_SCRIPTS/run.sh $DIR/ T $SEED > $DIR/$SEED.txt 2>&1 &"
  echo $COMMAND
  eval $COMMAND
done