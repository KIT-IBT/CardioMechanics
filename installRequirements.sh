#!/bin/bash

MPI_VERSION=v4.1.5
VTK_VERSION=v9.2.6
PETSC_VERSION=v3.19.1

# env variables for CardioMechanics
export kaRootDir=`pwd`
export THIRDPARTY_HOME=$kaRootDir/thirdparty

kernel=`uname -s`
if [ "$kernel" = "Linux" ]; then
    KAPREFIX="linux"
elif [ "$kernel" = "Darwin" ]; then
    KAPREFIX="macosx"
else
    KAPREFIX=""
fi
export KAPREFIX

mkdir -p $THIRDPARTY_HOME/src
prefixPath=${THIRDPARTY_HOME}/${KAPREFIX}

# compile OpenMPI
FC=gfortran-mp-12
export AUTOMAKE_JOBS=1
cd $THIRDPARTY_HOME/src
git clone --depth 1 --branch ${MPI_VERSION} https://github.com/open-mpi/ompi.git openmpi-${MPI_VERSION}
cd openmpi-${MPI_VERSION}
./autogen.pl
./configure --prefix=${prefixPath}/openMPI-64bit FC=${FC}
make 
make install

# compile vtk
cd $THIRDPARTY_HOME/src
git clone --depth 1 --branch ${VTK_VERSION} https://gitlab.kitware.com/vtk/vtk.git vtk-${VTK_VERSION}
cd vtk-${VTK_VERSION}
# we have to rename a variable to avoid a redeclaration under Linux
echo "Replacing HZ with H_Z in VTK..."
grep -rl HZ . | xargs sed -i 's/HZ/H_Z/g'
mkdir build
cd build
cmake .. -D CMAKE_INSTALL_PREFIX=${prefixPath}/vtk-${VTK_VERSION}
make
make install

# compile PETSc
cd $THIRDPARTY_HOME/src
git clone --depth 1 --branch ${PETSC_VERSION} https://gitlab.com/petsc/petsc.git petsc-${PETSC_VERSION}
cd petsc-${PETSC_VERSION}
unset PETSC_DIR
unset PETSC_ARCH
./configure --prefix=${prefixPath}/petsc-${PETSC_VERSION} --with-cmake=1 --with-mpi-dir=${prefixPath}/openMPI-64bit --download-superlu --download-superlu_dist --download-mumps --download-dmumps  --download-metis --download-parmetis --download-bison --download-ptscotch --download-scalapack --download-blacs --download-hypre --with-shared-libraries=0 --with-x=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 --with-debugging=no
make PETSC_DIR=${THIRDPARTY_HOME}/src/petsc-${PETSC_VERSION} PETSC_ARCH=arch-darwin-c-opt all
make PETSC_DIR=${THIRDPARTY_HOME}/src/petsc-${PETSC_VERSION} PETSC_ARCH=arch-darwin-c-opt install
make PETSC_DIR=${prefixPath}/petsc-${PETSC_VERSION} PETSC_ARCH="" check