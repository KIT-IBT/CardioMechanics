#!/bin/bash

MPI_VERSION=v4.1.5
VTK_VERSION=v9.2.6
PETSC_VERSION=v3.19.1

THIRDPARTY_HOME=$(pwd)/thirdparty
mkdir -p ${THIRDPARTY_HOME}/src
prefixPath=${THIRDPARTY_HOME}/install

# compile OpenMPI
export AUTOMAKE_JOBS=1
cd ${THIRDPARTY_HOME}/src
git clone --depth 1 --branch ${MPI_VERSION} https://github.com/open-mpi/ompi.git openmpi-${MPI_VERSION}
cd openmpi-${MPI_VERSION}
./autogen.pl
./configure --prefix=${prefixPath}/openMPI-64bit
make
make install

# compile VTK
cd ${THIRDPARTY_HOME}/src
git clone --depth 1 --branch ${VTK_VERSION} https://gitlab.kitware.com/vtk/vtk.git vtk-${VTK_VERSION}
cd vtk-${VTK_VERSION}
# rename a variable to avoid a redeclaration under Linux
echo "Replacing HZ with H_Z in VTK..."
grep -rl HZ . | xargs sed -i 's/HZ/H_Z/g'
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${prefixPath}/vtk-${VTK_VERSION}
make
make install

# compile PETSc
cd ${THIRDPARTY_HOME}/src
git clone --depth 1 --branch ${PETSC_VERSION} https://gitlab.com/petsc/petsc.git petsc-${PETSC_VERSION}
cd petsc-${PETSC_VERSION}
unset PETSC_DIR
unset PETSC_ARCH
./configure \
    --prefix=${prefixPath}/petsc-${PETSC_VERSION} \
    --with-cmake=1 \
    --with-mpi-dir=${prefixPath}/openMPI-64bit \
    --download-superlu --download-superlu_dist \
    --download-mumps --download-dmumps \
    --download-metis --download-parmetis \
    --download-bison --download-ptscotch \
    --download-scalapack --download-blacs \
    --download-hypre \
    --with-shared-libraries=0 --with-x=0 \
    --with-debugging=no \
    COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3
make all
make install

echo ""
echo "Done. Add the following to your shell configuration (.zshrc / .bashrc):"
echo ""
echo "  export PETSC_DIR=${prefixPath}/petsc-${PETSC_VERSION}"
echo "  export PETSC_ARCH="
echo "  export PATH=\"\$PATH:${prefixPath}/openMPI-64bit/bin\""
echo "  export VTK_DIR=${prefixPath}/vtk-${VTK_VERSION}"
