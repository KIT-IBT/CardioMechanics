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

The `docker/` directory contains Dockerfiles that build Open MPI, PETSc, and VTK from source with tested versions and can serve as a reference for building dependencies manually.

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
    --download-cmake \
    --download-fblaslapack --download-mumps --download-scalapack \
    --download-superlu --download-superlu_dist \
    --download-metis --download-parmetis --download-hypre \
    --with-shared-libraries=0 --with-x=0 \
    --with-debugging=0 \
    COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
make all && make install
```

`--download-fblaslapack` builds PETSc's own reference BLAS/LAPACK: linking a system OpenBLAS can
make MUMPS report the tangent stiffness as numerically singular on some CPUs.  `--download-cmake`
lets PETSc build a recent enough CMake for those solver packages that require one newer than the
system provides.

## Building

Configure and build using presets (recommended):

```sh
cmake --preset release   # Release build, output in _build/release/
cmake --build --preset release -j

cmake --preset debug     # Debug build, output in _build/debug/
cmake --build --preset debug -j
```

Or without presets:

```sh
cmake -S . -B _build/release -DCMAKE_BUILD_TYPE=Release
cmake --build _build/release -j
```

Binaries are placed in `_build/<preset>/bin/`.

Add the binaries and Python tools to your PATH:

```sh
export PATH="$PATH:/path/to/CardioMechanics/_build/release/bin"
export PATH="$PATH:/path/to/CardioMechanics/tools/python"
```

### Pinning PETSC_DIR per preset

If you maintain multiple PETSc builds (e.g. optimized and debug), create a local
`CMakeUserPresets.json` at the repo root to pin `PETSC_DIR` per preset without
modifying the committed `CMakePresets.json`:

```json
{
  "version": 3,
  "configurePresets": [
    {
      "name": "local-release",
      "inherits": "release",
      "binaryDir": "${sourceDir}/_build/release",
      "environment": { "PETSC_DIR": "/path/to/petsc-opt", "PETSC_ARCH": "" }
    },
    {
      "name": "local-debug",
      "inherits": "debug",
      "binaryDir": "${sourceDir}/_build/debug",
      "environment": { "PETSC_DIR": "/path/to/petsc-deb", "PETSC_ARCH": "" }
    }
  ],
  "buildPresets": [
    { "name": "local-release", "configurePreset": "local-release", "configuration": "Release" },
    { "name": "local-debug",   "configurePreset": "local-debug",   "configuration": "Debug" }
  ]
}
```

Add `CMakeUserPresets.json` to `.gitignore` to keep machine-local paths out of version control.

## Troubleshooting

**`nlohmann/json.hpp` not found**

Install the nlohmann-json package and expose its headers:
```sh
# Homebrew
brew install nlohmann-json
export CPATH=$CPATH:/opt/homebrew/opt/nlohmann-json/include
```
