# Created on May 5 2015 by ew095

cmake_minimum_required( VERSION 2.8 )

include( ../cmake/IBTDefault.cmake )

project( MathBasics )

IBTProject()

# Tell compiler where to look for libraries and header files
include_directories(
    .
    ../Lattice
    ../LatticeTools
)

# List of all libraries to be linked
kaAddLibs (
    Lattice
    z
)

kaMakeObject( kaTabularizedFunction kaTabularizedFunction.cpp )

kaMakeTools(
    MatrixNOpMatrixN
    MatrixNSet
    FormFunctionsTest,${kaTabularizedFunction}
    Tensor2Test,${kaTabularizedFunction}
    Tensor4Test,${kaTabularizedFunction}
)
