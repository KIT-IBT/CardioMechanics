# Created on Oct 16 2014 by Silvio Bauer

cmake_minimum_required( VERSION 3.5 )

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )

project( Lattice )
IBTProject()

# Special flags for libLattice
add_definitions( -DKALATTICEDEBUG=0 )

# I don't assume someone wants to fix these:
add_definitions( -Wno-unused-private-field -Wno-sign-compare -Wno-reorder )
add_definitions( -Wno-unused-variable )

# Tell compiler where to look for header files
include_directories( "../MathBasics" "." )

# Target for Lattice library
kaMakeLib( Lattice kaMachineOS.cpp kaBasicIO.cpp )

# Include library to the rest of the project
kaAddLibs( Lattice z )
