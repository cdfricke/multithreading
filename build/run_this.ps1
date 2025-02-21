# Build Script for multithreading project
# by Connor Fricke
# Build the project with this script from the project directory with 
#    powershell -File .\build\run_this.ps1

g++.exe -g .\cpp\rejection_sampling.cpp .\cpp-library\Calculus.cpp -o .\build\rejection_sampling.exe -I .\cpp-library

# to add arguments to the program, append $arguments to the end of the executable
$arguments = 1, 0.05, 100000    # SHAPE, DMIN, VOL
.\build\rejection_sampling.exe $arguments