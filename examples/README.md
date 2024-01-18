# Examples

Some examples for the use of CardioMechanics can be found here.

## benchmark2015

[benchmark2015](./benchmark2015) contains the three benchmark problems defined in the publication by [Land et al. (2015)](https://royalsocietypublishing.org/doi/10.1098/rspa.2015.0641).
The setup presented here uses P2 elements with a quasi-incompressible formulation in the constitutive model.
To run the problems, simply use the commands
```
CardioMechanics -settings Problem1.xml
CardioMechanics -settings Problem2.xml
CardioMechanics -settings Problem3.xml
```
Optionally, you can add `-debug` to get more output in the terminal or change specific PETSc options at runtime using e.g. `-ksp_atol <atol>`, `-ksp_rtol <rtol>`, `-ksp_max_it <its>` or `-ksp_monitor` to change and monitor convergence metrics of the underlying linear system solver.
Naturally, you can run the problems on multiple processes as well by preceding the above commands with `mpirun -np <n>`, where `n` is the number of processes you want to use.

> [!IMPORTANT]
> If you want accurate results according to the paper, you have to comment lines 45-46 and uncomment line 49 in file [CBTensionModel.h](/mechanics/src/CBTensionModel/CBTensionModel.h). Then recompile the code. Otherwise the active stress tensor is not the same as defined in the paper. This only concerns Problem3.

## EM01

This example realizes a simple electromechanical problem on a cuboid geometry with dimensions $20.0 \times 16.0 \times 12.0$ mm.
You can use this to familiarize yourself with the setup of all files required for electromechanical simulations.

The electrophysiology setup directly follows the problem definition of the [N-version benchmark by Niederer et al. (2011)](https://royalsocietypublishing.org/doi/full/10.1098/rsta.2011.0139) with an external stimulus provided in an $1.5 \times 1.5 \times 1.5$ mm area from a corner of the cube.
A sensor is placed at the location of P8 of the N-version benchmark to record $V_m$ and intracellular calcium, meaning you can directly compare to the results of the pure EP benchmark if you wish to do so.
For passive mechanics, the Holzapfel-Odgen material with $a = 330$ Pa, $a_{ff} = 18535$ Pa, $a_{ss} = 2564$ Pa, $a_{fs} = 417$ Pa, $b = 9.242$, $b_{ff} = 15.972$, $b_{ss} = 10.446$, $b_{fs} = 11.602$, and $\kappa = 10^6$ Pa was used.
Active stress is calculated by the Land17 model and zero normal displacement boundary conditions are placed at the boundaries where the stimulus is applied.

Before you can run the simulation itself, you first have to assemble the mass and stiffness matrix as well as the material vector for the monodomain equation by using the `BidomainMatrixGenerator`:
```
cd EM01/geoFiles
BidomainMatrixGenerator cube_0.25mm.vtu -o cube_0.25mm -mono ../materialFiles/materialIntra.def
```
Then run the simulation using
```
mpirun -np <n> CardioMechanics -settings M_1mm.xml
```
assuming you are in the EM01/settings directory.

## JPhys

This folder contains the files required to reproduce the electromechanically coupled whole heart simulations shown in the latest publication by [Gerach and Loewe (2024)](https://doi.org/10.1113/JP285022).
Due to file size limitations on GitHub, the geometry has to be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.10526554) first.
Make sure to download version 1.1 of the files.

Unzip the files and put the file `EP.vtu` into the directory `JPhys/geoFiles` and first assemble the FEM matrices using
```
cd JPhys/geoFiles
BidomainMatrixGenerator EP.vtu -o heart -v -mono ../materialFiles/materialIntra.def -tissue 1,2,3,72,73,74,75,76,77,78,79,80,32,33
```
followed by
```
cd ../settings
mpirun -np <n> CardioMechanics -settings scenario_I.xml
mpirun -np <n> CardioMechanics -settings scenario_II.xml
mpirun -np <n> CardioMechanics -settings scenario_III.xml
mpirun -np <n> CardioMechanics -settings scenario_IV.xml
mpirun -np <n> CardioMechanics -settings scenario_V.xml
```
to run the simulations.
To give you an idea on expected simulation times, simulating a single heart beat of 0.8s took around 4 hours on a Mac Pro equipped with a 28-core 2,5GHz Intel Xeon W processor.

> [!IMPORTANT]
> For simulation scenarios I and III you need to deactivate the troponin C feedback formulation in the source code of the ionic models by commenting the line containing `#define TRPN` in the files [OHaraRudyParameters.h](/electrophysiology/src/CellModel/OHaraRudyParameters.h) and [CourtemancheParameters.h](/electrophysiology/src/CellModel/CourtemancheParameters.h).