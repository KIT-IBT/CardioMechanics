add_library(Lattice STATIC kaMachineOS.cpp kaBasicIO.cpp)
target_include_directories(Lattice PUBLIC . ../MathBasics)
target_compile_definitions(Lattice PUBLIC KALATTICEDEBUG=0)
target_compile_options(Lattice PRIVATE
    -Wno-unused-private-field -Wno-sign-compare -Wno-reorder -Wno-unused-variable
)
target_link_libraries(Lattice PUBLIC cm_common_flags z)
