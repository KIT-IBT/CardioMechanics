# CardioMechanics

CardioMechanics is a simulation environment designed for (but not limited to) cardiac tissue electromechanics problems.
As of 09/2023, active development on this software is no longer persued.
Therefore, support is limited.
Even though this repository includes a standalone cardiac electrophysiology solver, we recommend using [openCARP](https://opencarp.org) instead if you are only interested in cardiac electrophysiology.

## Features
The following features are available:
* Single cell electrophysiology model test suite.
* Electrophysiology simulations using Monodomain and Bidomain.
* Passive and active tissue mechanics.
* Support for 1-, 2-, and 4-chamber heart geometries.
* Recovery of the reference configuration.
* Passive parameter optimization.
* Lumped parameter model of the circulatory system featuring heart valve dynamics.
* Electromechanically coupled simulations.
* XML-based simulation setup.

## Building from Source
See docs/BUILD.md

## Examples
See docs/EXAMPLES.md (coming soon)

## User manual
See docs/MANUAL.pdf (coming soon)

## License
This project is licensed under the GNU General Public License 3 - see the LICENSE.md file for details.

## Acknowledgments
CardioMechanics builds on
* [CMake](https://cmake.org)
* [PETSc](https://www.mcs.anl.gov/petsc/)
* [VTK](https://vtk.org)
* [Open MPI](https://www.open-mpi.org)

and the contributions of many former (PhD-)students at IBT.