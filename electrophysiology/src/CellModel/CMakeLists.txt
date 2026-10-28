add_library(CellModel STATIC
    Parameter.cpp
    ParameterSwitch.cpp
    ParameterLoader.cpp
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
    OHaraRudyIsoParameters.cpp
    TomekParameters.cpp
    TWorldParameters.cpp
    MitchellSchaefferParameters.cpp
    FitzhughNagumoParameters.cpp
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
    OHaraRudyIso.cpp
    Tomek.cpp
    TWorld.cpp
    MitchellSchaeffer.cpp
    FitzhughNagumo.cpp
)
target_include_directories(CellModel PUBLIC . ../Lattice)
target_compile_options(CellModel PRIVATE
    -Wno-unused-private-field -Wno-unused-variable
    -Wno-unneeded-internal-declaration -Wno-unused-function
)
target_link_libraries(CellModel PUBLIC Lattice z cm_common_flags)

add_executable(ElphyModelTest
    ElphyModelTest.cpp
    ../Lattice/kaVersion.cpp
    ../Lattice/kaMachineOS.cpp
)
target_compile_definitions(ElphyModelTest PRIVATE USE_EMT)
target_link_libraries(ElphyModelTest PRIVATE CellModel)

add_executable(CellModelTest
    CellModelTest.cpp
    ../Lattice/kaVersion.cpp
    ../Lattice/kaMachineOS.cpp
)
target_compile_definitions(CellModelTest PRIVATE USE_EMT USE_FMT)
target_link_libraries(CellModelTest PRIVATE CellModel)

add_executable(ForceModelTest
    ForceModelTest.cpp
    ../Lattice/kaVersion.cpp
    ../Lattice/kaMachineOS.cpp
)
target_compile_definitions(ForceModelTest PRIVATE USE_FMT)
target_link_libraries(ForceModelTest PRIVATE CellModel)

install(TARGETS
    ElphyModelTest
    CellModelTest
    ForceModelTest
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
