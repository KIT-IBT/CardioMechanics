add_executable(PETScVec2VTK
    ElectrophysiologyUtility/PETScVec2VTK.cpp
    ../electrophysiology/src/acCELLerate/PETScLSE.cpp
)
target_include_directories(PETScVec2VTK PRIVATE
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/acCELLerate
)
target_link_libraries(PETScVec2VTK PRIVATE
    Lattice
    PkgConfig::PETSc
    MPI::MPI_CXX
    ${VTK_LIBRARIES}
    z
    cm_common_flags
)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS PETScVec2VTK MODULES ${VTK_LIBRARIES})
endif()

add_executable(BidomainMatrixGenerator
    ElectrophysiologyUtility/BidomainMatrixGenerator.cpp
    ElectrophysiologyUtility/Verbosity.cpp
)
target_include_directories(BidomainMatrixGenerator PRIVATE
    ElectrophysiologyUtility
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/MathBasics
    ../electrophysiology/src/acCELLerate
    ../mechanics/src/CBMaterial
    ../mechanics/src/MathPack
)
target_link_libraries(BidomainMatrixGenerator PRIVATE
    Lattice
    PkgConfig::PETSc
    MPI::MPI_CXX
    ${VTK_LIBRARIES}
    z
    cm_common_flags
)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS BidomainMatrixGenerator MODULES ${VTK_LIBRARIES})
endif()

add_executable(VTPExtractSurfaceFromMesh
    MechanicsUtility/VTPExtractSurfaceFromMesh.cpp
)
target_include_directories(VTPExtractSurfaceFromMesh PRIVATE
    ../electrophysiology/src/Lattice
    ../electrophysiology/src/MathBasics
)
target_link_libraries(VTPExtractSurfaceFromMesh PRIVATE
    Lattice
    ${VTK_LIBRARIES}
    cm_common_flags
)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS VTPExtractSurfaceFromMesh MODULES ${VTK_LIBRARIES})
endif()

add_executable(ExtractSurfaceNodesFromVTU MechanicsUtility/ExtractSurfaceNodesFromVTU.cpp)
target_link_libraries(ExtractSurfaceNodesFromVTU PRIVATE
    ${VTK_LIBRARIES}
    cm_common_flags
)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS ExtractSurfaceNodesFromVTU MODULES ${VTK_LIBRARIES})
endif()

add_executable(FixNodes MechanicsUtility/FixNodes.cpp)
target_link_libraries(FixNodes PRIVATE cm_common_flags)

add_executable(ConvertT4toT10 MechanicsUtility/ConvertT4toT10.cpp)
target_link_libraries(ConvertT4toT10 PRIVATE cm_common_flags)

install(TARGETS
    PETScVec2VTK
    BidomainMatrixGenerator
    VTPExtractSurfaceFromMesh
    ExtractSurfaceNodesFromVTU
    FixNodes
    ConvertT4toT10
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
