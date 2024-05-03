ARG VTK_IMAGE=cardiomechanics/thirdparty-vtk
ARG OPENMPI_PETSC_IMAGE=cardiomechanics/thirdparty-openmpi-petsc

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

COPY --from=$VTK_IMAGE ${prefixPath}/vtk-${VTK_VERSION} ${prefixPath}/vtk-${VTK_VERSION}
COPY --from=$OPENMPI_PETSC_IMAGE ${prefixPath}/openMPI-64bit ${prefixPath}/${prefixPath}/openMPI-64bit
COPY --from=$OPENMPI_PETSC_IMAGE ${prefixPath}/petsc-${PETSC_VERSION} ${prefixPath}/petsc-${PETSC_VERSION}

