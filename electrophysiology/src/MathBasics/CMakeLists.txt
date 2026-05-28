cmake_minimum_required(VERSION 3.20)
project(MathBasics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT TARGET Lattice)
    add_subdirectory(../Lattice ${CMAKE_CURRENT_BINARY_DIR}/Lattice)
endif()

include_directories(. ../Lattice ../LatticeTools)

add_library(kaTabularizedFunction OBJECT kaTabularizedFunction.cpp)

foreach(tool MatrixNOpMatrixN MatrixNSet)
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} PRIVATE Lattice z)
endforeach()

foreach(tool FormFunctionsTest Tensor2Test Tensor4Test)
    add_executable(${tool} ${tool}.cpp $<TARGET_OBJECTS:kaTabularizedFunction>)
    target_link_libraries(${tool} PRIVATE Lattice z)
endforeach()
