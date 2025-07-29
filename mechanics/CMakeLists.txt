cmake_minimum_required(VERSION 3.5)

if( DEFINED ENV{kaRootDir} )
    set( kaRootDir $ENV{kaRootDir} )
else ( )
    message( WARNING "Environment variable kaRootDir should be defined.")
endif( DEFINED ENV{kaRootDir} )

include( ${kaRootDir}/cmake/IBTDefault.cmake )
set ( OLD_OPT_LOCAL_VALUE ${IBT_USE_OPT_LOCAL} )
set ( IBT_USE_OPT_LOCAL OFF CACHE BOOL "abc" FORCE )

project(CardioMechanics)
IBTProject()

add_definitions( -DUSE_MUMPS -Wno-c++11-narrowing -DACLT_NO_MAIN)


find_package(VTK REQUIRED)
if (${VTK_MAJOR_VERSION} GREATER 8)
    include_directories(${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
else ()
    include(${VTK_USE_FILE})
endif ()

find_package(PETSC REQUIRED)
kaAddLibdirs(${PETSC_LIBDIRS})

find_package(MPI REQUIRED)

include_directories(
    ./src/CBCavity
    ./src/CBConstitutiveModel
    ./src/CBData
    ./src/CBElement
    ./src/CBTensionModel
    ./src/CBFormulation
    ./src/CBMaterial
    ./src/CBModel
    ./src/CBSolver
    ./src/CBSolverPlugins
    ./src/CBStatus
    ./src/CardioMechanics
    ./src/DCTK
    ./src/Frizzle
    ./src/MathPack
    ./src/ParameterMap
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/MathBasics
    ../electrophysiology/src/CellModel
    ../electrophysiology/src/acCELLerate
    /opt/local/include/
    ${MPI_INCLUDE_PATH}
    ${PETSC_INCLUDES}
)

kaBuildLinkLibrariesString(PETSC_LIBRARIES_STRING "${PETSC_LIBRARIES}")
kaBuildLinkLibrariesString(VTK_LIBRARIES_STRING "${VTK_LIBRARIES}")
kaBuildLinkLibrariesString(MPI_LIBRARIES_STRING "${MPI_LIBRARIES}")


SET( SRCS
    src/CardioMechanics/CardioMechanics.cpp
    src/CBCavity/CBCircLeftHeart.cpp    
    src/CBCavity/CBCircOneVentricle.cpp
    src/CBCavity/CBCircIsovolumetric.cpp
    src/CBCavity/CBCircTwoVentricles.cpp
    src/CBCavity/CBCircWholeHeart.cpp
    src/CBCavity/CBCircWholeHeartValves.cpp
    src/CBCavity/CBCircWholeHeartValvesDynamic.cpp
    src/CBCavity/CBCircWindkessel3.cpp
    src/CBCavity/CBCircWindkessel3AsymptoticPressure.cpp
    src/CBCavity/CBCirculationCavity.cpp
    src/CBConstitutiveModel/CBConstitutiveModel.cpp
    src/CBConstitutiveModel/CBConstitutiveModelFactory.cpp
    src/CBConstitutiveModel/CBConstitutiveModelGuccione.cpp
    src/CBConstitutiveModel/CBConstitutiveModelUsyk.cpp
    src/CBConstitutiveModel/CBConstitutiveModelNeoHooke.cpp
    src/CBConstitutiveModel/CBConstitutiveModelHolzapfel.cpp
    src/CBConstitutiveModel/CBConstitutiveModelMooneyRivlin.cpp
    src/CBData/CBDataCtrl.cpp
    src/CBData/CBDataFromFile.cpp
    src/CBData/CBDataFromFunction.cpp
    src/CBData/CBDataPerMaterial.cpp
    src/CBData/CBFileManager.cpp
    src/CBElement/CBElement.cpp
    src/CBElement/CBElementAdapter.cpp
    src/CBElement/CBElementCavity.cpp
    src/CBElement/CBElementContactMaster.cpp
    src/CBElement/CBElementContactSlave.cpp
    src/CBElement/CBElementContactRobin.cpp
    src/CBElement/CBElementFactory.cpp
    src/CBElement/CBElementSolid.cpp
    src/CBElement/CBElementSolidT10.cpp
    src/CBElement/CBElementSolidT4.cpp
    src/CBElement/CBElementSurfaceT3.cpp
    src/CBElement/CBElementSurfaceT6.cpp
    src/CBTensionModel/CBTensionFactory.cpp
    src/CBTensionModel/CBTensionFromFunction.cpp
    src/CBTensionModel/CBTensionModelLumens.cpp
    src/CBTensionModel/CBTensionModelBestel.cpp
    src/CBTensionModel/CBTensionModelLand17.cpp
    src/CBTensionModel/CBTensionModelTanH.cpp
    src/CBFormulation/CBFormulationTotalLagrangian.cpp
    src/CBMaterial/CBMaterial.cpp
    src/CBMaterial/CBMaterialFactory.cpp
    src/CBModel/CBModel.cpp
    src/CBModel/CBModelExporter.cpp
    src/CBModel/CBModelExporterVTK.cpp
    src/CBModel/CBModelLoaderTetgen.cpp
    src/CBSolver/CBSolver.cpp
    src/CBSolver/CBTiming.cpp
    src/CBSolver/CBSolverEquilibrium.cpp
    src/CBSolver/CBSolverNewmarkBeta.cpp
    src/CBSolverPlugins/CBApplyPressure.cpp
    src/CBSolverPlugins/CBApplyPressureFromFunction.cpp
    src/CBSolverPlugins/CBApplyPressureFromFunctionNodeExport.cpp
    src/CBSolverPlugins/CBReferenceRecovery.cpp
    src/CBSolverPlugins/CBCirculation.cpp
    src/CBSolverPlugins/CBContactHandling.cpp
    src/CBSolverPlugins/CBRobinBoundary.cpp
    src/CBSolverPlugins/CBRobinBoundaryGeneral.cpp
    src/CBSolverPlugins/CBLoadUnloadedState.cpp
    src/CBSolverPlugins/CBSolverPluginCavities.cpp
    src/CBSolverPlugins/CBSolverPluginFactory.cpp
    src/CBSolverPlugins/CBacCELLerate.cpp
    src/CBStatus/CBStatus.cpp
    src/DCTK/DCCtrl.cpp
    src/DCTK/DCCtrlPETSc.cpp
    src/DCTK/DCMatrix.cpp
    src/DCTK/DCMatrixPETSc.cpp
    src/DCTK/DCVector.cpp
    src/DCTK/DCVectorPETSc.cpp
    src/Frizzle/filesystem.cpp
    src/Frizzle/stringtools.cpp
    src/ParameterMap/ParameterMap.cpp
    ../electrophysiology/src/acCELLerate/acCELLerate.cpp
    ../electrophysiology/src/acCELLerate/acltConditions.h
    ../electrophysiology/src/acCELLerate/acltProjectFile.h
    ../electrophysiology/src/acCELLerate/acltSensors.h
    ../electrophysiology/src/acCELLerate/PETScConditions.h
    ../electrophysiology/src/acCELLerate/PETScLSE.cpp
    ../electrophysiology/src/acCELLerate/timer.h
    ../electrophysiology/src/acCELLerate/acltCellModel.cpp
    ../electrophysiology/src/Lattice/kaMachineOS.cpp
    ../electrophysiology/src/Lattice/kaBasicIO.cpp
)


if (NOT ${CMAKE_GENERATOR} STREQUAL "Xcode")
    kaMakeObject( Common ${SRCS} )
else (NOT ${CMAKE_GENERATOR} STREQUAL "Xcode")
    file(GLOB_RECURSE all_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
    kaMakeObject( Common ${SRCS} ${all_headers} )
    message(STATUS " generate Xcode " )
endif (NOT ${CMAKE_GENERATOR} STREQUAL "Xcode")

kaBuildLinkLibrariesString( CommonLibs "${PETSC_LIBRARIES};${VTK_LIBRARIES};z" )

kaBuildLinkLibrariesString( CELLMODEL_LIBS "CellModel" )

kaMakeTool( CardioMechanics
    ./src/CardioMechanics/main.cpp
    ${Common}
    ${CommonLibs}
    ${CELLMODEL_LIBS}
)

set ( IBT_USE_OPT_LOCAL ${OLD_OPT_LOCAL_VALUE} CACHE BOOL "abc" FORCE )
