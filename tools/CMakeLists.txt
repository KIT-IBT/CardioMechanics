cmake_minimum_required( VERSION 2.8 )

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )

set ( OLD_OPT_LOCAL_VALUE ${IBT_USE_OPT_LOCAL} )
set ( IBT_USE_OPT_LOCAL OFF CACHE BOOL "abc" FORCE )

project(tools)
IBTProject()

# Tell compiler where to look for header files
include_directories (
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/MathBasics
    ../electrophysiology/src/acCELLerate
    ../mechanics/src/MathPack
    ../mechanics/src/CBMaterial
    ./ElectrophysiologyUtility
    ./MechanicsUtility
)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
string( REPLACE ";" ",-l" VTK_LIBS "-l${VTK_LIBRARIES}" )

find_package(PETSC REQUIRED)
kaAddLibdirs(${PETSC_LIBDIRS})
include_directories(${PETSC_INCLUDES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

kaBuildLinkLibrariesString(PETSC_LIBRARIES_STRING "${PETSC_LIBRARIES}")
kaBuildLinkLibrariesString(VTK_LIBRARIES_STRING "${VTK_LIBRARIES}")
kaBuildLinkLibrariesString(MPI_LIBRARIES_STRING "${MPI_LIBRARIES}")
kaBuildLinkLibrariesString(CommonLibs "${PETSC_LIBRARIES};${VTK_LIBRARIES};z")

kaMakeTool( PETScVec2VTK
	./ElectrophysiologyUtility/PETScVec2VTK.cpp
	${VTK_LIBRARIES_STRING}
    ${PETSC_LIBRARIES_STRING}
)

kaMakeTool( BidomainMatrixGenerator
	./ElectrophysiologyUtility/BidomainMatrixGenerator.cpp
    ./ElectrophysiologyUtility/Verbosity.cpp
	${CommonLibs}
)

kaMakeTool( VTPExtractSurfaceFromMesh
    ./MechanicsUtility/VTPExtractSurfaceFromMesh.cpp
    ${VTK_LIBRARIES_STRING}
)

kaMakeTool( FixNodes
    ./MechanicsUtility/FixNodes.cpp
)

kaMakeTool( ConvertT4toT10
    ./MechanicsUtility/ConvertT4toT10.cpp
)
