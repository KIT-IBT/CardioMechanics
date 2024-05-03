# Create docker image for CardioMechanics

```
docker build -t cardiomechanics/thirdparty-vtk -f Dockerfile-thirdparty-vtk .
docker build -t cardiomechanics/thirdparty-openmpi-petsc -f Dockerfile-thirdparty-openmpi-petsc .
docker build -t cardiomechanics/cardiomechanics -f Dockerfile ..
```

