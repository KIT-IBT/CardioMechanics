cmake_minimum_required(VERSION 3.20)
project(CardioMechanics LANGUAGES CXX)

# VTK 9+ requires C++17; set globally to avoid flag conflicts when targets propagate their standards.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Homebrew on macOS (Apple Silicon: /opt/homebrew, Intel: /usr/local)
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH /opt/homebrew /usr/local)
endif()

# ccache: auto-enable when found
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    message(STATUS "ccache found: ${CCACHE_PROGRAM}")
endif()

find_package(MPI REQUIRED)
find_package(VTK REQUIRED)
if(VTK_VERSION VERSION_LESS "9.0")
    include(${VTK_USE_FILE})
endif()

find_package(PkgConfig REQUIRED)
# Locate PETSc pkg-config file via PETSC_DIR env var if set
if(DEFINED ENV{PETSC_DIR})
    if(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
        set(_petsc_pc "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig")
    else()
        set(_petsc_pc "$ENV{PETSC_DIR}/lib/pkgconfig")
    endif()
    set(ENV{PKG_CONFIG_PATH} "${_petsc_pc}:$ENV{PKG_CONFIG_PATH}")
endif()
pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET PETSc)

# Interface target that propagates compiler flags and debug macros to all targets
add_library(cm_common_flags INTERFACE)
target_compile_options(cm_common_flags INTERFACE
    -pedantic -Wall -Werror=vla -Wno-extra-semi
)
target_compile_definitions(cm_common_flags INTERFACE
    $<$<CONFIG:Debug>:KADEBUG=1;KAMATERIALDEBUG=1>
    $<$<NOT:$<CONFIG:Debug>>:KADEBUG=0;KAMATERIALDEBUG=0;NDEBUG>
)
if(APPLE)
    target_compile_definitions(cm_common_flags INTERFACE __DARWIN__)
elseif(UNIX)
    target_compile_definitions(cm_common_flags INTERFACE _DEFAULT_SOURCE)
endif()

add_subdirectory(electrophysiology/src/Lattice)
add_subdirectory(electrophysiology/src/CellModel)
add_subdirectory(electrophysiology/src/acCELLerate)
add_subdirectory(mechanics)
add_subdirectory(tools)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
