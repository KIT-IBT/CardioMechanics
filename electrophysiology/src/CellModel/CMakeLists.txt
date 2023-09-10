# Created on Feb 4 2014 by ew095

cmake_minimum_required( VERSION 2.8 )

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )

project( CellModel )
IBTProject()

add_definitions( -Wno-unused-private-field -Wno-unused-variable )
add_definitions( -Wno-unneeded-internal-declaration -Wno-unused-function )

# Tell compiler where to look for header files
include_directories(. ../Lattice)

kaMakeObject( ParameterHandling
    Parameter.cpp ParameterSwitch.cpp
    ParameterLoader.cpp
)

# While it is possible to use GLOB *Parameters.cpp here, please don't.
# If a code file changes, the icmake build system would have no way of knowing
# that it should re-generate.
kaMakeObject ( Parameters
    CourtemancheParameters.cpp
    TenTusscherEtAlParameters.cpp
    TenTusscher2Parameters.cpp
    SachseEtAlParameters.cpp
    MaleckarEtAlParameters.cpp
    KoivumaekiEtAlParameters.cpp
    GrandiEtAlVentricleParameters.cpp
    GrandiEtAlAtriumParameters.cpp
    BeelerReuterParameters.cpp
    HybridModelParameters.cpp
    FabbriParameters.cpp
    Land17Parameters.cpp
    HimenoParameters.cpp
    OHaraRudyParameters.cpp
    MitchellSchaefferParameters.cpp
    FitzhughNagumoParameters.cpp
)

kaMakeObject( CellModels
    CourtemancheEtAl.cpp
    TenTusscher.cpp
    TenTusscher2.cpp
    SachseEtAl.cpp
    MaleckarEtAl.cpp
    KoivumaekiEtAl.cpp
    GrandiEtAlVentricle.cpp
    GrandiEtAlAtrium.cpp
    BeelerReuter.cpp
    HybridModel.cpp
    FabbriEtAl.cpp
    Land17.cpp
    Himeno.cpp
    OHaraRudy.cpp
    MitchellSchaeffer.cpp
    FitzhughNagumo.cpp
    )

kaMakeObject( kaVersion ../Lattice/kaVersion.cpp )
kaMakeObject( kaParser ../Lattice/kaParser.cpp )
kaMakeObject( nrutil nrutil.cpp )

kaMakeLib( CellModel ${Parameters} ${CellModels} ${ParameterHandling})
kaAddLibs(
    ${CellModel}
    Lattice
    z
)

kaMakeTools(
    ElphyModelTest.cpp,-DUSE_EMT,${kaVersion},../Lattice/kaMachineOS.cpp
    CellModelTest.cpp,-DUSE_EMT,-DUSE_FMT,${kaVersion},../Lattice/kaMachineOS.cpp
    ForceModelTest.cpp,-DUSE_FMT,${kaVersion},../Lattice/kaMachineOS.cpp
    APAnalysis.cpp,${nrutil}
)
