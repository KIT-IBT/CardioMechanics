ARG PETSC_IMAGE=cardiomechanics/thirdparty-petsc

# Build on top of the PETSc image so OpenMPI and PETSc come from a single source;
# this rules out any MPI version drift between the two images.
FROM ${PETSC_IMAGE}

# VTK and the remaining build dependencies come from the package manager. VTK is not
# linked into PETSc, so it has no cross-image consistency requirement.
RUN apt-get update && apt-get install -y \
        clang libvtk9-dev pkg-config zlib1g-dev freeglut3-dev mesa-common-dev \
    && rm -rf /var/lib/apt/lists/*

ENV CC=clang
ENV CXX=clang++
ENV kaRootDir=/CardioMechanics

WORKDIR /CardioMechanics
COPY . .

RUN cmake -S . -B _build
RUN cmake --build _build --parallel 4

ENV PATH=$PATH:/CardioMechanics/_build/bin/linux
WORKDIR /CardioMechanics/_build/bin/linux
