set(SRCS
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
    ../electrophysiology/src/acCELLerate/PETScLSE.cpp
    ../electrophysiology/src/acCELLerate/acltCellModel.cpp
    ../electrophysiology/src/Lattice/kaMachineOS.cpp
    ../electrophysiology/src/Lattice/kaBasicIO.cpp
)

add_library(CardioMechanicsCore STATIC ${SRCS})
target_include_directories(CardioMechanicsCore PUBLIC
    src/CBCavity
    src/CBConstitutiveModel
    src/CBData
    src/CBElement
    src/CBTensionModel
    src/CBFormulation
    src/CBMaterial
    src/CBModel
    src/CBSolver
    src/CBSolverPlugins
    src/CBStatus
    src/CardioMechanics
    src/DCTK
    src/Frizzle
    src/MathPack
    src/ParameterMap
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/MathBasics
    ../electrophysiology/src/CellModel
    ../electrophysiology/src/acCELLerate
)
# USE_MUMPS: required by PETSc MUMPS solver integration (compiled into PETSc via --download-mumps).
# ACLT_NO_MAIN: suppress acCELLerate's main() since its sources are compiled into this library.
target_compile_definitions(CardioMechanicsCore PRIVATE USE_MUMPS ACLT_NO_MAIN)
target_compile_options(CardioMechanicsCore PRIVATE -Wno-c++11-narrowing)
target_link_libraries(CardioMechanicsCore PUBLIC
    CellModel
    Lattice
    PkgConfig::PETSc
    MPI::MPI_CXX
    ${VTK_LIBRARIES}
    z
    cm_common_flags
)

add_executable(CardioMechanics src/CardioMechanics/main.cpp)
target_link_libraries(CardioMechanics PRIVATE CardioMechanicsCore)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS CardioMechanics MODULES ${VTK_LIBRARIES})
endif()
