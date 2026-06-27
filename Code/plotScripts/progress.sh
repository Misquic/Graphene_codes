#!/bin/bash -l
set -u

echo $#
if [ $# -ne 1 ]; then
  echo "illegal number of parameters"
  exit 1
fi

DIR=$1

TOTAL=$(( $(wc -l < "$DIR/commands.txt") - 1 ))

COMPLETED=$(ls "$DIR"dirs/*/Transmissions.csv | wc -l)
# COMPLETED_COUNT=$(echo "$COMPLETED" | wc -l)

PERCENT=$(( $COMPLETED * 100 / $TOTAL))
echo "Completed: $COMPLETED / $TOTAL = $PERCENT %"
