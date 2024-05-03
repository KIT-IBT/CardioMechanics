ARG VTK_IMAGE=cardiomechanics/thirdparty-vtk
ARG OPENMPI_PETSC_IMAGE=cardiomechanics/thirdparty-openmpi-petsc

FROM ${VTK_IMAGE} as vtkImage
FROM ${OPENMPI_PETSC_IMAGE} as openmpiPetscImage

FROM ubuntu:22.04

WORKDIR /CardioMechanics
COPY . .

RUN apt update && apt install -y build-essential gfortran python3 git cmake autoconf libtool
ENV MPI_VERSION=v4.1.5
ENV VTK_VERSION=v9.2.6
ENV PETSC_VERSION=v3.19.1

# env variables for CardioMechanics
ENV kaRootDir=/CardioMechanics
ENV THIRDPARTY_HOME=$kaRootDir/thirdparty
RUN mkdir -p ${THIRDPARTY_HOME}/src
ENV prefixPath=${THIRDPARTY_HOME}/install

COPY --from=vtkImage ${prefixPath}/vtk-${VTK_VERSION} ${prefixPath}/vtk-${VTK_VERSION}
COPY --from=openmpiPetscImage ${prefixPath}/openMPI-64bit ${prefixPath}/openMPI-64bit
COPY --from=openmpiPetscImage ${prefixPath}/petsc-${PETSC_VERSION} ${prefixPath}/petsc-${PETSC_VERSION}

ENV VTK_DIR=${prefixPath}/vtk-${VTK_VERSION}
ENV PETSC_DIR=${prefixPath}/petsc-${PETSC_VERSION}
ENV PATH=$PATH:${prefixPath}/openMPI-64bit/bin
ENV LD_LIBRARY_PATH=${prefixPath}/openMPI-64bit/lib:${prefixPath}/vtk-${VTK_VERSION}/lib

RUN apt install -y pkg-config mesa-common-dev libopenblas-dev zlib1g-dev freeglut3-dev
RUN cmake -S. -B_build
RUN cmake --build _build --parallel 4