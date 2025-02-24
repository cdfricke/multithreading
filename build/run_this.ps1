# Build Script for multithreading project
# by Connor Fricke
# Build the project with this script from the project directory with 
#    powershell -File .\build\run_this.ps1

# * COMPILING *
# compile with g++, show all warnings and elevate all warnings to errors
# add my library to include path, and add the .cpp files for that library to the compilation files list
# the -D option specifies a define, by removing it, the program won't write data to any files
$cppFiles = "cpp/rejection_sampling.cpp", "cpp-library/Calculus.cpp", "cpp-library/Stats.cpp"
$includeDir = "cpp-library"
$executable = ".\build\rejection_sampling.exe"
g++.exe -Wall -Werror -O2 `
-I $includeDir -D WRITE_DATA `
$cppFiles `
-o $executable

# * VARS *
$SPHERE = 0
#$CUBE = 1
#$CUBOID_A = 2
#$CUBOID_B = 3
$CE3 = 0
#$CE4 = 1

# * RUN PROGRAM *
# to use custom arguments to the program, append $arguments to the end of the executable
# and likewise remove it to use default parameters
$arguments = $SPHERE, $CE3, 0.05, 3.0, 10.0, 2500.0    # SHAPE, DIST, DMIN, DMAX, DEPTH, AREA
.\build\rejection_sampling.exe $arguments