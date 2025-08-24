@echo off

call build.bat %1

if exist ".\main.exe" ( 
    main.exe
) else (
    echo "build failed"
)