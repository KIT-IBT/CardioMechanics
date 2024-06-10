# Create docker image for CardioMechanics

From the current directory, build the Docker images for dependencies and the main Docker image:

```
docker build -t cardiomechanics/thirdparty-vtk -f Dockerfile-thirdparty-vtk .
docker build -t cardiomechanics/thirdparty-openmpi-petsc -f Dockerfile-thirdparty-openmpi-petsc .
docker build -t cardiomechanics/cardiomechanics -f Dockerfile ..
```

The main image `cardiomechanics/cardiomechanics` is based on Ubuntu 22.04.
