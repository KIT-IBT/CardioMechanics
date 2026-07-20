add_executable(acCELLerate
    acCELLerate.cpp
    acltCellModel.cpp
    PETScLSE.cpp
    ../MathBasics/kaTabularizedFunction.cpp
)
target_include_directories(acCELLerate PRIVATE . ../Lattice ../MathBasics ../CellModel)
target_compile_options(acCELLerate PRIVATE
    -Wno-unused-private-field -Wno-unused-variable -Wno-reorder
)
target_link_libraries(acCELLerate PRIVATE
    CellModel
    Lattice
    PkgConfig::PETSc
    MPI::MPI_CXX
    ${VTK_LIBRARIES}
    z
    cm_common_flags
)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    vtk_module_autoinit(TARGETS acCELLerate MODULES ${VTK_LIBRARIES})
endif()

install(TARGETS acCELLerate RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
