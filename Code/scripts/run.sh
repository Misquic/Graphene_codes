#!/bin/bash

source build.sh $@

if [ -e "main" ]; then  
    ./main
else 
    echo "build failed"
    exit 1
fi