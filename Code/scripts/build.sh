#!/bin/bash

if [ -e "main" ]; then 
    rm main
fi;

startDir=$(pwd)
project=$(basename $startDir)

opt=$(echo "$1" | tr '[:upper:]' '[:lower:]')

cd "../../"
case "$opt" in
    -d|d|debug)
        echo "building Debug"
        cmake --build build/Debug --target $project -- -j 8 # -- passes next arguments to make
        ;;
    -r|r|release)
        echo "building Release"
        cmake --build build/Release --target $project -- -j 8 # -- passes next arguments to make
        ;;
    *)
        echo "Invalid build option 1: <release/debug>"
        cd $startDir
        exit 1
        ;;
esac

cd $startDir
