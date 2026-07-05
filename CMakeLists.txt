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
# Request only the VTK modules whose headers are used, so find_package(VTK) does
# not resolve unrelated dependencies (e.g. Qt via Ubuntu's GUI-enabled VTK).
find_package(VTK REQUIRED COMPONENTS
    CommonCore CommonMath CommonDataModel
    FiltersCore FiltersGeometry FiltersModeling FiltersVerdict
    IOXML IOLegacy IOGeometry
    RenderingCore RenderingOpenGL2
)

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
# PETSc's public Libs is only -lpetsc; its solver dependencies (MUMPS, ScaLAPACK,
# BLAS/LAPACK, MPI-Fortran) live in Libs.private. Linking only -lpetsc lets those
# symbols bind to other providers pulled in transitively (e.g. VTK's Accelerate
# BLAS), which corrupts MUMPS's parallel distributed-RHS solve and segfaults.
# Link the full static closure so they resolve against PETSc's own libraries.
set_property(TARGET PkgConfig::PETSc APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES "${PETSc_STATIC_LDFLAGS}")
# Some of those dependencies (e.g. HYPRE) have @rpath install names, so PETSc's
# library directory must be on the runtime search path to resolve them.
list(APPEND CMAKE_BUILD_RPATH ${PETSc_LIBRARY_DIRS})
list(APPEND CMAKE_INSTALL_RPATH ${PETSc_LIBRARY_DIRS})
if(APPLE)
    list(APPEND CMAKE_BUILD_RPATH /opt/homebrew/lib)
    list(APPEND CMAKE_INSTALL_RPATH /opt/homebrew/lib)
endif()

# Interface target that propagates compiler flags and debug macros to all targets
add_library(cm_common_flags INTERFACE)
target_compile_options(cm_common_flags INTERFACE
    -pedantic -Wall -Werror=vla -Wno-extra-semi
)
target_compile_definitions(cm_common_flags INTERFACE
    $<$<CONFIG:Debug>:KADEBUG=1;KAMATERIALDEBUG=1>
    $<$<NOT:$<CONFIG:Debug>>:KADEBUG=0;KAMATERIALDEBUG=0;NDEBUG>
)
# Default location of the bundled model parameter files (electrophysiology/data),
# so binaries resolve them without a runtime kaRootDir env var. Setting kaRootDir
# still overrides it.
target_compile_definitions(cm_common_flags INTERFACE CM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
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
