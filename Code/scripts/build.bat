@echo off
if exist ".\main.exe" ( 
    rm main.exe
)

Set startDir=%cd%

@REM cutting path so only name of project is saved in Project variable
for %%A in (%startDir%) do Set "Project=%%~nA"
cd ../../

if "%1" == "-Debug" (
    GOTO DEBUG
) else if "%1" == "-debug" (
    GOTO DEBUG
) else if "%1" == "-d" (
    GOTO DEBUG
) else if "%1" == "-D" (
    GOTO DEBUG
) else if "%1" == "Debug" (
    GOTO DEBUG
) else if "%1" == "debug" (
    GOTO DEBUG
) else if "%1" == "d" (
    GOTO DEBUG
) else if "%1" == "D" (
    GOTO DEBUG
) else if "%1" == "-Release" (
    GOTO RELEASE
) else if "%1" == "-release" (
    GOTO RELEASE
) else if "%1" == "-r" (
    GOTO RELEASE
) else if "%1" == "-R" (
    GOTO RELEASE
) else if "%1" == "release" (
    GOTO RELEASE
) else if "%1" == "Release" (
    GOTO RELEASE
) else if "%1" == "r" (
    GOTO RELEASE
) else if "%1" == "R" (
    GOTO RELEASE
) 
echo "Invalid build option 1: <release/debug>"
exit /b 1

:DEBUG
    echo "building Debug"
    cmake --build build/Debug --target %Project%  -- -j 8
    GOTO CHECKERROR
:RELEASE
    echo "building Release"
    cmake --build build/Release --target %Project% -- -j 8
    GOTO CHECKERROR

:CHECKERROR
if %ERRORLEVEL% neq 0 (
    echo Build failed. Exiting...
    exit /b 1
)

cd /d %startDir%
