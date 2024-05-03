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

COPY --from=cardiomechanics/thirdparty-vtk ${prefixPath}/vtk-${VTK_VERSION} ${prefixPath}/vtk-${VTK_VERSION}
COPY --from=cardiomechanics/thirdparty-openmpi-petsc ${prefixPath}/openMPI-64bit ${prefixPath}/${prefixPath}/openMPI-64bit
COPY --from=cardiomechanics/thirdparty-openmpi-petsc ${prefixPath}/petsc-${PETSC_VERSION} ${prefixPath}/petsc-${PETSC_VERSION}

