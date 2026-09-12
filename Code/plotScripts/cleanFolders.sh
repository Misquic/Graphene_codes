#!/bin/bash -l

numFolderToDelete=0

echo "Folders to clean:"
for dir in "$1"/*/; do
  if [ -e "$dir/dirs/" ]; then
    echo "$dir/dirs/"
    numFolderToDelete=$(( $numFolderToDelete + 1 ))
  fi
  if [ -e "$dir/outs/" ]; then
    echo "$dir/outs/"
    numFolderToDelete=$(( $numFolderToDelete + 1 ))
  fi
done

if [ $numFolderToDelete -eq 0 ]; then
  echo "Nothing to delete, exiting"
  exit 0
fi

echo "Proceed? [y/n]"
read PROCEED

if [ $PROCEED != "y" ] && [ $PROCEED != "Y" ] && [ $PROCEED != "1" ]; then
  echo "exiting"
  exit 0
fi

echo "deleting"
for dir in "$1"/*/; do
  rm -rf "$dir/dirs/" "$dir/outs/" &
done