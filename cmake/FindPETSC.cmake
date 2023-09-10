# Look for PETSc install location respecting PETSC_ARCH

include ( kaPrefix )

set ( PETSC_ARCH $ENV{PETSC_ARCH} CACHE STRING "Subdir of PETSC_DIR if applicable")

if ( PETSC_ARCH )
    set ( petsc_arch_sl "${PETSC_ARCH}/" )
endif ( PETSC_ARCH )

find_path ( PETSC_DIR "${petsc_arch_sl}include/petsc.h" HINTS ENV PETSC_DIR
    PATHS $ENV{THIRDPARTY_HOME} DOC "Directory under which PETSc is installed.")

if ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}conf/PETScConfig.cmake" )

    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set ( PETSC_BLAS_LIB "-framework Accelerate" CACHE STRING "Blas for PETSc")
        set ( PETSC_LAPACK_LIB "-framework Accelerate" CACHE STRING "Lapack for PETSc" )
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    include ( "${PETSC_DIR}/${petsc_arch_sl}conf/PETScConfig.cmake" )

    set ( PETSC_LIBRARIES "-L${PETSC_DIR}/${petsc_arch_sl}lib" "-lpetsc" ${PETSC_PACKAGE_LIBS})
    set ( PETSC_INCLUDES "${PETSC_DIR}/${petsc_arch_sl}include" ${PETSC_PACKAGE_INCLUDES} )
    set ( PETSC_FOUND YES)

    # newer versions of petsc (>=3.6) have the petscconfig.cmake in a different directory
elseif ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScBuildInternal.cmake" OR EXISTS "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScConfig.cmake" )
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set ( PETSC_BLAS_LIB "-framework Accelerate" CACHE STRING "Blas for PETSc")
        set ( PETSC_LAPACK_LIB "-framework Accelerate" CACHE STRING "Lapack for PETSc" )
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    if ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScBuildInternal.cmake" )
        include ( "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScBuildInternal.cmake" )
    else ()
        include ( "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScConfig.cmake" )
    endif ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}lib/petsc/conf/PETScBuildInternal.cmake" )

    set ( PETSC_LIBRARIES "-L${PETSC_DIR}/${petsc_arch_sl}lib" "-lpetsc" ${PETSC_PACKAGE_LIBS})
    set ( PETSC_INCLUDES "${PETSC_DIR}/${petsc_arch_sl}include" ${PETSC_PACKAGE_INCLUDES} )
    set ( PETSC_FOUND YES)

elseif ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}include/petsc.h" )
    #if ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}include/petsc.h" )

    set ( PETSC_USE_STATIC ON CACHE BOOL "Statically link PETSc." )

    set ( PETSC_FOUND YES)
    set ( PETSC_INCLUDES "${PETSC_DIR}/${petsc_arch_sl}include")

    set( PETSC_LIBDIR "${PETSC_DIR}/${petsc_arch_sl}lib")

    set ( PKGCFG_EXTRA_ARGS "" )
    if ( ${PETSC_USE_STATIC} )
        set ( PKGCFG_EXTRA_ARGS "--static" )
    endif ( ${PETSC_USE_STATIC} )

    set ( ENV{PKG_CONFIG_PATH} "${PETSC_LIBDIR}/pkgconfig" )
    find_program(PKGCONFIG_EXECUTABLE NAMES pkg-config )

    execute_process( COMMAND ${PKGCONFIG_EXECUTABLE} PETSc ${PKGCFG_EXTRA_ARGS} --cflags-only-I OUTPUT_VARIABLE tmpinc )
    execute_process( COMMAND ${PKGCONFIG_EXECUTABLE} PETSc ${PKGCFG_EXTRA_ARGS} --libs OUTPUT_VARIABLE tmplib )

    string ( REPLACE "-I" "" tmpinc ${tmpinc} )
    string ( REPLACE " " ";" PETSC_INCLUDES ${tmpinc} )
    string ( STRIP "${tmplib}" tmplib )
    string ( REPLACE " " ";" tmplib ${tmplib} )

    foreach ( libentry ${tmplib} )
        if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            if ( ${libentry} MATCHES "blas" )
                set ( libentry "Accelerate" )
            elseif ( ${libentry} MATCHES "lapack" )
                set ( libentry "Accelerate" )
            elseif ( ${libentry} MATCHES "libpthread" )
                set ( libentry "pthread" )
            endif ( ${libentry} MATCHES "blas" )
        endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

        if ( ${libentry} MATCHES "^-L.+" ) # Library dir
            string ( SUBSTRING ${libentry} 2 -1 libentry )
            if ( EXISTS ${libentry} )
                list ( APPEND PETSC_LIBDIRS ${libentry} )
            else ( )
                message ( WARNING "PETSc Linkdir ${libentry} doesn't exist." )
            endif ( EXISTS ${libentry} )
        elseif ( ${libentry} MATCHES "^-l.+" ) #lib name
            string ( SUBSTRING ${libentry} 2 -1 libentry )
            list ( APPEND PETSC_LIBRARIES ${libentry} )
        elseif (IS_ABSOLUTE ${libentry}) # Absolute path to library
            if ( EXISTS ${libentry} )
                list ( APPEND PETSC_LIBRARIES ${libentry} )
            else ( ) # Add library name w/o path so that it can be set from cache
                get_filename_component ( libentry "${libentry}" NAME )
                string ( REGEX REPLACE "^lib" "" libentry "${libentry}" )
                string ( REGEX REPLACE "\\.a$" "" libentry "${libentry}" )
                string ( REGEX REPLACE "\\.dylib$" "" libentry "${libentry}" )
                string ( REGEX REPLACE "\\.so$" "" libentry "${libentry}" )
                list ( APPEND PETSC_LIBRARIES ${libentry} )
            endif ( EXISTS ${libentry} )
        else ( ) # Just a library
            list ( APPEND PETSC_LIBRARIES ${libentry} )
        endif ( ${libentry} MATCHES "^-L.+" )
    endforeach ( libentry ${tmplib} )

    list ( APPEND PETSC_LIBRARIES "scalapack" )
    list ( REMOVE_DUPLICATES PETSC_LIBDIRS )
    list ( REMOVE_DUPLICATES PETSC_LIBRARIES )

    message ( STATUS "PETSc libs: ${PETSC_LIBRARIES}" )
    message ( STATUS "PETSc libdirs: ${PETSC_LIBDIRS}" )

    #pkg_search_module( PETSC ${PETSC_FIND_REQUIRED} ${PETSC_FIND_QUIET} PETSc libpetsc petsc PETSC )

else (  )

    set ( PETSC_FOUND NO )

#endif ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}include/petsc.h" )
endif ( EXISTS "${PETSC_DIR}/${petsc_arch_sl}conf/PETScConfig.cmake" )

include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args ( PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDES )
