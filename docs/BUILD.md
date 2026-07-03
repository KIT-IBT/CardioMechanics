# Building CardioMechanics from source

Building from source was tested on Linux and Intel based Mac systems.
ARM based Macs may require some changes.

## Requirements

All build dependencies except PETSc are available from the system package manager.
We recommend using one whenever possible.

* C and C++ compilers (e.g. gcc/g++ or clang/clang++)
* [CMake](https://cmake.org) (>= 3.21 for the preset build)
* [zlib](https://zlib.net)
* [gfortran](https://gcc.gnu.org/fortran/)
* [git](https://git-scm.com)
* [make](https://www.gnu.org/software/make/)
* [VTK](https://vtk.org)
* [Open MPI](https://www.open-mpi.org)
* [PETSc](https://www.mcs.anl.gov/petsc/) — built from source, see below
* [Python3](https://www.python.org) (optional, for the provided tools)

On Debian/Ubuntu:
```sh
sudo apt install build-essential gfortran cmake git zlib1g-dev \
                 libopenmpi-dev openmpi-bin libvtk9-dev qtbase5-dev
```
`qtbase5-dev` is needed because Ubuntu's VTK is built with Qt modules, so
`find_package(VTK)` resolves Qt5 even though CardioMechanics uses none of VTK's GUI
modules. On macOS the equivalent packages come from Homebrew:
```sh
brew install cmake open-mpi vtk gfortran nlohmann-json
```

## Building PETSc from source

PETSc is the only dependency that has to be built from source: the packaged builds do
not ship all the solvers CardioMechanics needs. The configuration below has most recently
been tested with PETSc `v3.24.0`.

CardioMechanics uses the MUMPS and SuperLU direct solvers, so PETSc is configured with
`--download-mumps` / `--download-superlu` / `--download-superlu_dist`. The remaining
downloads are their dependencies: `--download-scalapack` (required by MUMPS),
`--download-metis` / `--download-parmetis` (parallel ordering), `--download-hypre` (a
runtime-selectable preconditioner), and `--download-fblaslapack` (reference BLAS/LAPACK).
`--download-cmake` builds a private CMake for those package builds: some of them now
require CMake >= 3.26, newer than some distributions ship (Ubuntu 22.04 has 3.22).

Pick an install prefix (`$HOME/software` below), then configure, build, install, and check:
```sh
PETSC_VERSION=v3.24.0
PETSC_PREFIX=$HOME/software/petsc-$PETSC_VERSION

git clone --depth 1 --branch $PETSC_VERSION https://gitlab.com/petsc/petsc.git petsc-$PETSC_VERSION
cd petsc-$PETSC_VERSION
unset PETSC_DIR PETSC_ARCH
./configure --prefix=$PETSC_PREFIX \
    --download-cmake \
    --download-fblaslapack --download-mumps --download-scalapack \
    --download-superlu --download-superlu_dist \
    --download-metis --download-parmetis --download-hypre \
    --with-shared-libraries=0 --with-x=0 --with-debugging=0 \
    COPTFLAGS=-O2 CXXOPTFLAGS=-O2 FOPTFLAGS=-O2
    # add e.g. --with-mpi-dir=$(brew --prefix open-mpi) if configure does not
    # find your MPI on PATH
# configure prints the exact make command, including the in-tree PETSC_ARCH
# (arch-darwin-c-opt / arch-linux-c-opt). Run it, then install and check:
make PETSC_DIR=$PWD PETSC_ARCH=arch-darwin-c-opt all
make PETSC_DIR=$PWD PETSC_ARCH=arch-darwin-c-opt install
make PETSC_DIR=$PETSC_PREFIX PETSC_ARCH= check
```
The resulting `--prefix` install is a single-arch tree, so at build time you point
`PETSC_DIR` at `$PETSC_PREFIX` and leave `PETSC_ARCH` empty (see below).

## Environment variables

Add the following to your shell configuration (`.bashrc` / `.zshrc`):
```sh
export kaRootDir=$HOME/CardioMechanics            # repo root; also read at runtime
export PETSC_DIR=$HOME/software/petsc-v3.24.0     # the --prefix you installed PETSc into
export PETSC_ARCH=                                # empty for a --prefix install
```
`kaRootDir` is required regardless of how you build — the binaries also read it **at
runtime** to locate their bundled data files. We assume here that you cloned the repository
into your `$HOME` directory; adjust `kaRootDir` if you chose a different location.

Add the built binaries and the Python tools to your `PATH` (replace `macosx` with `linux`
on a Linux machine):
```sh
export PATH="$PATH:$kaRootDir/_build/release/bin/macosx"
export PATH="$PATH:$kaRootDir/tools/python"
```
Open MPI is already on your `PATH` through the package manager.

## Building with CMake presets

The repository ships a `CMakePresets.json` with `release` and `debug` presets. They set
`kaRootDir` to the source directory automatically and pick up `PETSC_DIR` / `PETSC_ARCH`
from your environment, so with the variables above already exported you can build straight
away:
```sh
cmake --preset release
cmake --build --preset release
```
Use the `debug` preset in the same way for a debug build. Preset builds place binaries in
`_build/<release|debug>/bin/macosx` (or `.../linux`).

If you would rather keep the machine-specific PETSc paths out of your shell environment,
provide them in an uncommitted `CMakeUserPresets.json` next to `CMakePresets.json` that
inherits a base preset instead:
```json
{
  "version": 3,
  "configurePresets": [
    {
      "name": "local-release",
      "inherits": "release",
      "cacheVariables": {
        "PETSC_DIR": "/home/you/software/petsc-v3.24.0",
        "PETSC_ARCH": ""
      }
    }
  ],
  "buildPresets": [
    { "name": "local-release", "configurePreset": "local-release" }
  ]
}
```
Then build with `cmake --preset local-release` and `cmake --build --preset local-release`.
