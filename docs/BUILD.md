# Building CardioMechanics from source

First, you have to make sure all requirements for building CardioMechanics are met.
Afterwards, CardioMechanics can be build using CMake.
Building from source was tested on Linux and Intel based Mac systems. 
ARM based Macs may require some changes.

## Requirements

The following requirements have to be installed before trying to build CardioMechanics from source.
We recommend using a package manager (we use [macports](https://www.macports.org) on our macOSX systems) whenever possible.
* C and C++ compilers (e.g. gcc/g++ or clang/clang++)
* [CMake](https://cmake.org)
* [zlib](https://zlib.net)
* [gfortran](https://gcc.gnu.org/fortran/)
* [git](https://git-scm.com)
* [make](https://git-scm.com)
* [PETSc](https://www.mcs.anl.gov/petsc/)
* [VTK](https://vtk.org)
* [Open MPI](https://www.open-mpi.org)
* [Python3](https://www.python.org) (optional, if you want to use some of the provided tools)

With [installRequirements.sh](/installRequirements.sh) we provide a script to compile [Open MPI](https://www.open-mpi.org), [PETSc](https://www.mcs.anl.gov/petsc/), and [VTK](https://vtk.org) from source with the most recently tested versions to ensure compatibility.
Building with CMake as described in the next step requires the location and version of the tools as set in the script.
If you do want to use alternative locations/versions you have to link them as required.
By default, the script [installRequirements.sh](/installRequirements.sh) compiles all components using a single process, which takes a considerable amount of time.
If you want to speed up this process, use the command `make` with the option `-j X` where `X` is the number of processes you want to use.
Additionally, adjust `export AUTOMAKE_JOBS=X` specifically for [Open MPI](https://www.open-mpi.org).

## Building using CMake

Before continuing with compiling CardioMechanics using [CMake](https://cmake.org), add the following environmental variables to your systems configuration file ( e.g. .bashrc or .zshrc)

If you installed the requirements with homebrew add:
```
export kaRootDir=$HOME/CardioMechanics
export THIRDPARTY_HOME=$kaRootDir/thirdparty
export PETSC_DIR=/opt/homebrew/opt/petsc
```
If you are using [installRequirements.sh](/installRequirements.sh) add:
```
export kaRootDir=$HOME/CardioMechanics
export THIRDPARTY_HOME=$kaRootDir/thirdparty
export PETSC_DIR=$THIRDPARTY_HOME/macosx
export PETSC_ARCH=petsc-v3.19.1
```
In case you did not use either option to install the requirements add the paths `kaRootDir` and `THIRDPARTY_HOME` in the same way as proposed above and edit the path for `PETSC_DIR` according to your actual path.

Next add the location of the executables to your PATH variable (replace macosx with linux if you are on a linux machine).
```
PATH="$PATH:$THIRDPARTY_HOME/macosx/openMPI-64bit/bin"
PATH="$PATH:$kaRootDir/_build/bin/macosx"
PATH="$PATH:$kaRootDir/tools/python"
export PATH
```
We assume here that you copied the repository to your `$HOME` directory.
If you chose a different root directory, adjust the `kaRootDir` variable accordingly.
If you used a different PETSc version in [installRequirements.sh](/installRequirements.sh), adjust the `PETSC_ARCH` variable.
Now run 
```
cmake -S . -B _build
```
to create the `_build` folder and compile the code using
```
cmake --build _build
```

## Troubleshooting

If you experience the following error 
```
nlohmann/json.hpp is not found (fatal error)
```
you will have to add the nlohmann-json manually. 
If you used homebrew to install the requirements you will have to include the following line in your .zshrc:
```
export CPATH=$CPATH:/opt/homebrew/opt/nlohmann-json/include
```