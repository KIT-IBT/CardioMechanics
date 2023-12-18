cmake_minimum_required( VERSION 2.8 )

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )

set ( OLD_OPT_LOCAL_VALUE ${IBT_USE_OPT_LOCAL} )
set ( IBT_USE_OPT_LOCAL OFF CACHE BOOL "abc" FORCE )

project( acCELLerate )
IBTProject()

add_definitions( -Wno-unused-private-field -Wno-unused-variable -Wno-reorder )

# Tell compiler where to look for header files
include_directories(
    .
    ../Lattice
    ../MathBasics
    ../CellModel
)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
string( REPLACE ";" ",-l" VTK_LIBS "-l${VTK_LIBRARIES}" )

set ( CELLMODEL_LIBS -lCellModel )

find_package(PETSC REQUIRED)
kaAddLibdirs(${PETSC_LIBDIRS})
include_directories(${PETSC_INCLUDES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# List of all libraries to be linked
kaAddLibs(
    ${PETSC_LIBRARIES}
    Lattice
    z
)

kaAddObjects( ./PETScLSE.cpp )
kaAddObjects( ../MathBasics/kaTabularizedFunction.cpp )

kaMakeTools(
    acCELLerate.cpp,acltCellModel.cpp,${CELLMODEL_LIBS}
)

set ( IBT_USE_OPT_LOCAL ${OLD_OPT_LOCAL_VALUE} CACHE BOOL "abc" FORCE )