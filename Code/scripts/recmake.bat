@echo off

rm -r ./build/
mkdir build
cd ./build/
mkdir Debug
mkdir Release
cd ..

cmake -S . -B ./build/Release -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake -S . -B ./build/Debug -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug
