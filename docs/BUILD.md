# Building CardioMechanics from source

First, you have to make sure all requirements for building CardioMechanics are met.
Afterwards, CardioMechanics can be built using CMake.
Building from source was tested on Linux and macOS (Intel and Apple Silicon).

## Requirements

The following requirements have to be installed before trying to build CardioMechanics from source.
We recommend using a package manager (e.g. [Homebrew](https://brew.sh) on macOS or `apt` on Ubuntu) whenever possible.
* C and C++ compilers (e.g. gcc/g++ or clang/clang++)
* [CMake](https://cmake.org) ≥ 3.20
* [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
* [zlib](https://zlib.net)
* [gfortran](https://gcc.gnu.org/fortran/)
* [git](https://git-scm.com)
* [PETSc](https://petsc.org/) (built with `--prefix` so that a `PETSc.pc` pkg-config file is generated)
* [VTK](https://vtk.org) ≥ 9
* [Open MPI](https://www.open-mpi.org)
* [Python3](https://www.python.org) (optional, for the tools in `tools/python/`)

With [installRequirements.sh](/installRequirements.sh) we provide a script to compile [Open MPI](https://www.open-mpi.org), [PETSc](https://petsc.org/), and [VTK](https://vtk.org) from source with tested versions.
By default the script uses a single process, which is slow.
Speed it up by passing `-j X` to `make` and setting `export AUTOMAKE_JOBS=X` for Open MPI.

## Environment variables

Only `PETSC_DIR` (and optionally `PETSC_ARCH`) are required so that CMake can locate the PETSc pkg-config file.

```sh
export PETSC_DIR=/path/to/petsc          # prefix install directory
export PETSC_ARCH=                        # leave empty if PETSc was installed with --prefix
```

If `PETSC_ARCH` is non-empty, the pkg-config file is expected at `$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/PETSc.pc`.
If `PETSC_ARCH` is empty, it is expected at `$PETSC_DIR/lib/pkgconfig/PETSc.pc`.

VTK and Open MPI are found via CMake's standard search paths.
On macOS with Homebrew, no extra variables are needed — Homebrew's prefix (`/opt/homebrew` or `/usr/local`) is searched automatically.

### PETSc build example

PETSc must be configured with `--prefix` to generate the pkg-config file.  The options below
enable the solver packages used by CardioMechanics:

```sh
cd /path/to/petsc-source
./configure \
    --prefix=/path/to/petsc/install \
    --download-superlu --download-superlu_dist \
    --download-mumps --download-dmumps \
    --download-bison --download-ptscotch \
    --download-scalapack --download-blacs \
    --download-metis --download-parmetis \
    --download-hypre \
    --with-shared-libraries=0 --with-x=0 \
    --with-debugging=0 \
    COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
make all && make install
```

## Building

Configure and build using presets (recommended):

```sh
cmake --preset default   # RelWithDebInfo, build dir _build
cmake --build --preset default -j
```

Or without presets:

```sh
cmake -S . -B _build
cmake --build _build -j
```

Binaries are placed in `_build/bin/`.

Add the binaries and Python tools to your PATH:

```sh
export PATH="$PATH:/path/to/CardioMechanics/_build/bin"
export PATH="$PATH:/path/to/CardioMechanics/tools/python"
```

Available build presets: `default` (RelWithDebInfo), `debug`, `release`.

## Troubleshooting

**`nlohmann/json.hpp` not found**

Install the nlohmann-json package and expose its headers:
```sh
# Homebrew
brew install nlohmann-json
export CPATH=$CPATH:/opt/homebrew/opt/nlohmann-json/include
```
