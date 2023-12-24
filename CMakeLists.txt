cmake_minimum_required(VERSION 2.8)

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )

project( CardioMechanics )
set ( USE_LIBCXX ON )

add_subdirectory( electrophysiology/src/Lattice )
add_subdirectory( electrophysiology/src/CellModel )
add_subdirectory( electrophysiology/src/acCELLerate )
add_subdirectory( mechanics )
add_subdirectory( tools )

if(CMAKE_BUILD_TYPE)
    message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )
endif(CMAKE_BUILD_TYPE)