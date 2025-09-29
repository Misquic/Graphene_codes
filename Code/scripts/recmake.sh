#!/bin/bash

if [ -d "./build" ]; then
    rm -r ./build/
fi
mkdir build

mkdir build/Debug
mkdir build/Release

cmake -S . -B ./build/Release -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake -S . -B ./build/Debug -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
