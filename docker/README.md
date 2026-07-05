# Create docker image for CardioMechanics

The main image is based on Ubuntu 22.04. VTK and Open MPI come from the package manager;
only PETSc is built from source, in a separate cached image. The main image builds `FROM`
the PETSc image so both share the same Open MPI.

From the current directory, build the PETSc image and then the main image:

```
docker build -t cardiomechanics/thirdparty-petsc -f Dockerfile-thirdparty-petsc .
docker build -t cardiomechanics/cardiomechanics -f Dockerfile ..
```
