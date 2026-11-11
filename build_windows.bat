@echo off
setlocal
cl /nologo /std:c++17 /EHsc /O2 /utf-8 main.cpp /Fe:retaining_wall.exe
if errorlevel 1 (
    echo Build failed.
    exit /b 1
)
echo Build succeeded: retaining_wall.exe
