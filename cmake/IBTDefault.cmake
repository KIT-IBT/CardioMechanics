# CMake IBT include file
# Created Feb 04, 2015 by ew095
#
# Exports the following functions:
#
# kaAddLibs ( <list> )
# -- The libraries specified in <list> will be linked to all targets in the
#    current project.
#
# kaAddObjects ( <list of sourcefiles> )
# -- The objects created from each sourcefile will be linked to all targets
#    in the current project.
#
# kaMakeLib ( <name> <srcs...> )
# -- Create static library named <name> from <srcs...>
#
# kaMakeSharedLib ( <name> <srcs...> )
# -- Create shared library named <name> from <srcs...>
#
# kaMakeModuleLib ( <name> <srcs...> )
# -- Create .bundle library named <name> from <srcs...>
#
# kaMakeObject ( <name> <srcs...> )
# -- Create object files for <srcs...> that can be linked. The resulting object
#    must be referenced (e.g. in `kaMakeTools()`) by `${name}` (without quotes).
#
# kaMakeTools ( <name.cpp>[,<extra-srcs>]* ... )
# -- Create executables for tools from a list of sources.
#    A tool will be named the basename of its's first sourcefile, e.g. `name`
#    for `name.cpp`. Each list item  may contain multiple source files, objects
#    (from other kaMake... methods) and flags. All of these must be separated
#    by commas.
#    Example:
#        kaMakeObject ( obj obj_src1.cpp obj_src2.cpp )
#        kaMakeTools ( tool1.cpp,-DDEBUG tool2.cpp,extra.cpp tool3.cpp,${obj} )
#
# kaBuildLinkLibrariesString ( <output> <list> )
# -- Join the values in `list` using "-l" and store the resulting string in `output` in order to pass
#    it to e.g. kaMakeTools
#
#
# The following macro is defined:
#
# IBTOverrideCmakeFlag ( flag value type help )
# -- Create a CMakeCache entry named `IBT_flag` that has default value `value`
#    and is of type `type`. The cache description is `help`. The value may be
#    modified/overridden by the user.
#    IMPORTANT: This flag will override the flag CMAKE_flag.
#    Motivation: Some CMAKE_* flags are set to a default value in the cache.
#    There is no way to overwrite them without using `FORCE`, but then the user
#    can not modify the value.
#
# JOIN(<list> <glue> <output>)
# -- Joins the values in `list` using `glue` and stores the result in `output`



cmake_minimum_required( VERSION 3.5 )

# Set the preferred C and C++ compilers

#if( DEFINED CMAKE_SYSTEM_NAME )
#    message( FATAL_ERROR "This file (IBTDefault.cmake) must be included "
#        "BEFORE the project(...) or IBTProject(...) directive in your "
#        "CMakeLists.txt. Prefer IBTProject(...) " )
#endif( DEFINED CMAKE_SYSTEM_NAME )

find_program ( IBT_C_COMPILER NAMES $ENV{CC} "clang" "cc" "gcc" DOC "C compiler to use")

if ( EXISTS ${IBT_C_COMPILER} )
    set ( CMAKE_C_COMPILER ${IBT_C_COMPILER} )
endif ( EXISTS ${IBT_C_COMPILER} )

find_program ( IBT_CXX_COMPILER NAMES $ENV{CXX} "clang++" "c++" "g++" DOC "C++ compiler to use")

if ( EXISTS ${IBT_CXX_COMPILER} )
    set ( CMAKE_CXX_COMPILER ${IBT_CXX_COMPILER} )
endif ( EXISTS ${IBT_CXX_COMPILER} )

# This macro creates IBT_flag with a default but changeable value.
# CMAKE_flag will always be set to IBT_flag.
# This way, we can override a CMAKE default without FORCE, so that the user
# can change it.
macro ( IBTOverrideCmakeFlag flag value type help)
    set ( IBT_${flag} ${value} CACHE ${type} ${help} )
    set ( CMAKE_${flag} ${IBT_${flag}} CACHE INTERNAL ${help} FORCE )
endmacro ( IBTOverrideCmakeFlag )

# Use this macro instead after the project(..) directive
macro ( IBTProject )

    if( NOT DEFINED CMAKE_SYSTEM_NAME )
        message( FATAL_ERROR "Use the IBTProject() macro after the project(...) directive." )
    endif( NOT DEFINED CMAKE_SYSTEM_NAME )

    if( DEFINED ENV{kaRootDir} )
        set( kaRootDir $ENV{kaRootDir} )
    else ( )
        message( WARNING "Environment variable kaRootDir should be defined.")
    endif( DEFINED ENV{kaRootDir} )

    if( DEFINED $ENV{THIRDPARTY_HOME} )
        set( THIRDPARTY_HOME $ENV{THIRDPARTY_HOME} )
    else ( )
        set( THIRDPARTY_HOME ${kaRootDir}/thirdparty )
    endif( DEFINED $ENV{THIRDPARTY_HOME} )

    set( CMAKE_MODULE_PATH "${kaRootDir}/cmake" )

    include ( kaPrefix ) # sets /macosx/ or /linuxELF/ etc

    # Pre-set the MPI compiler path so that FindMPI finds our MPI
    set ( MPI_C_COMPILER "${THIRDPARTY_HOME}/${KAPREFIX}/openMPI-64bit/bin/mpicc" )
    if ( NOT EXISTS ${MPI_C_COMPILER} )
        unset ( MPI_C_COMPILER )
    else ( )
        set ( MPI_CXX_COMPILER "${THIRDPARTY_HOME}/${KAPREFIX}/openMPI-64bit/bin/mpicxx" )
    endif ( NOT EXISTS ${MPI_C_COMPILER} )
    
    option ( IBT_USE_CCACHE "Use ccache to cache build results" OFF )
    if ( IBT_USE_CCACHE )
      find_program(CCACHE_FOUND ccache)
      if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
      endif(CCACHE_FOUND)
    endif ( IBT_USE_CCACHE )
    
    

    option ( IBT_USE_LIBCXX "Use clang's libc++ instead of the GNU stdlibc++." ON)
    if ( NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
        set ( IBT_USE_LIBCXX OFF CACHE BOOL "Can only use libc++ with clang++" FORCE )
    endif ( NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )

    IBTOverrideCmakeFlag ( INSTALL_PREFIX "${kaRootDir}" PATH "Install under this directory." )
    IBTOverrideCmakeFlag ( C_FLAGS "-std=c11 -pedantic -Wall" STRING "C compiler flags." )
    IBTOverrideCmakeFlag ( CXX_FLAGS "-std=c++11 -pedantic -Wall -Werror=vla -Wno-extra-semi" STRING "C++ compiler flags." )
    
    if ( IBT_USE_LIBCXX AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
        add_compile_options( -stdlib=libc++ )
        set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi" )
    endif ( IBT_USE_LIBCXX AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )

    set ( DEBUG_FLAGS "-g -DKADEBUG=1 -DMATERIAL_PROPERTY_DEBUG=1" )
    set ( RELEASE_FLAGS "-DNDEBUG -DKADEBUG=0 -DMATERIAL_PROPERTY_DEBUG=0" )

    IBTOverrideCmakeFlag ( C_FLAGS_DEBUG "-O0 ${DEBUG_FLAGS}" STRING "Flags for configuration 'Debug'." )
    IBTOverrideCmakeFlag ( C_FLAGS_MINSIZEREL "-Os ${RELEASE_FLAGS}" STRING "Flags for configuration 'MinSizeRel'." )
    IBTOverrideCmakeFlag ( C_FLAGS_RELEASE "-O3 ${RELEASE_FLAGS}" STRING "Flags for configuration 'Release'." )
    IBTOverrideCmakeFlag ( C_FLAGS_RELWITHDEBINFO "-O3 -g ${RELEASE_FLAGS}" STRING "Flags for configuration 'RelWithDebInfo'." )

    IBTOverrideCmakeFlag ( CXX_FLAGS_DEBUG "-O0 ${DEBUG_FLAGS}" STRING "Flags for configuration 'Debug'." )
    IBTOverrideCmakeFlag ( CXX_FLAGS_MINSIZEREL "-Os ${RELEASE_FLAGS}" STRING "Flags for configuration 'MinSizeRel'." )
    IBTOverrideCmakeFlag ( CXX_FLAGS_RELEASE "-O3 ${RELEASE_FLAGS}" STRING "Flags for configuration 'Release'." )
    IBTOverrideCmakeFlag ( CXX_FLAGS_RELWITHDEBINFO "-O3 -g ${RELEASE_FLAGS}" STRING "Flags for configuration 'RelWithDebInfo'." )

    # Using rpath makes no sense in our setup. Might lead to executables that work
    # after building but not after setup.
    IBTOverrideCmakeFlag ( SKIP_RPATH OFF BOOL "Do not set rpath to linking path when building." )
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
    #SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

    if( NOT CMAKE_BUILD_TYPE )
        set( CMAKE_BUILD_TYPE RelWithDebInfo )
    endif( NOT CMAKE_BUILD_TYPE )
    
    set( IBT_BIN_PREFIX bin/${KAPREFIX} CACHE PATH "Installation directory for executables" )
    set( IBT_LIB_PREFIX lib/${KAPREFIX} CACHE PATH "Installation directory for libraries" )
    # set( IBT_INC_PREFIX include/${KAPREFIX} CACHE PATH "Installation directory for headers (currently unused)" )
    # set( IBT_DATA_DIR share/${KAPREFIX}/${PROJECT_NAME}  CACHE PATH "Installation directory for data files (currently unused)" )

    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${IBT_BIN_PREFIX} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${IBT_LIB_PREFIX} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${IBT_LIB_PREFIX} )

    set ( CMAKE_PREFIX_PATH "${kaRootDir}" )
    #set ( CMAKE_LIBRARY_ARCHITECTURE ${KAPREFIX} )
    #foreach ( prefix ${CMAKE_PREFIX_PATH} )
    #    link_directories( ${prefix}/lib/${CMAKE_LIBRARY_ARCHITECTURE} )
    #    link_directories( ${prefix}/${CMAKE_LIBRARY_ARCHITECTURE}/lib )
    #endforeach ( prefix ${CMAKE_PREFIX_PATH} )

    set ( IBT_USE_OPT_LOCAL ON CACHE BOOL "Scan /opt/local/ for thirdparty libraries." )
    if ( IBT_USE_OPT_LOCAL )
        set ( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/opt/local" )
        if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
            set ( CMAKE_FRAMEWORK_PATH ${CMAKE_FRAMEWORK_PATH} "/opt/local/Library/Frameworks" )
        endif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    endif ( IBT_USE_OPT_LOCAL )

    set ( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${THIRDPARTY_HOME}" "${THIRDPARTY_HOME}/${KAPREFIX}" )

    define_property ( DIRECTORY PROPERTY KA_LIBRARIES
        BRIEF_DOCS "List of libraries for ${PROJECT_NAME}"
        FULL_DOCS "These libraries will be linked to all targets in project ${PROJECT_NAME}." )

    define_property ( DIRECTORY PROPERTY KA_LINKDIRS
        BRIEF_DOCS "List of link directories for ${PROJECT_NAME}"
        FULL_DOCS "Libraries for ${PROJECT_NAME} will be searched for in these directories." )

endmacro ( IBTProject )

function( kaMakeObject objname )
    set ( ctr 0 )
    set ( originalObjname ${objname} )
    while ( TARGET ${objname} )
        # Object of this name has already been defined.
        # We just give it a new name. objname is overwritten anyway
        # and must always be used as ${objname} by the caller.

        math ( EXPR ctr ${ctr}+1 )
        set ( objname ${originalObjname}${ctr} )
    endwhile ( TARGET ${objname} )

    add_library( ${objname} OBJECT ${ARGN} )
    set( ${originalObjname} "$<TARGET_OBJECTS:${objname}>" PARENT_SCOPE )
endfunction( kaMakeObject )

macro ( get_doc_libname var lib )
    get_filename_component ( docname "${lib}" NAME_WE )
    if ( NOT ${docname} MATCHES "^lib.+" )
        set ( docname "lib${docname}" )
    endif ( NOT ${docname} MATCHES "^lib.+" )

    set ( ${var} "${docname}" )
endmacro ( get_doc_libname )

macro ( add_library_helper libname libdirs libraries )
    if ( ${${libname}} MATCHES "^-[lL].+" ) # strip prefixes
        string ( SUBSTRING ${${libname}} 2 -1 ${libname} )
    endif ( ${${libname}} MATCHES "^-[lL].+" )


    string( SUBSTRING ${${libname}} 0 1 firstChar)
    string( COMPARE EQUAL ${firstChar} "/" isCanonical )

    if ( TARGET "${${libname}}" )
        list ( APPEND ${libraries} "${${libname}}" )
    elseif ( (NOT ${isCanonical}) AND IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/${${libname}}" )
        list ( APPEND ${libdirs} "${CMAKE_CURRENT_LIST_DIR}/${${libname}}" )
    elseif ( ${isCanonical} AND IS_DIRECTORY "${${libname}}" )
        list ( APPEND ${libdirs} "${${libname}}" )
    elseif ( (NOT ${isCanonical}) AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/${${libname}}" ) # is a file
        list ( APPEND ${libraries} "${CMAKE_CURRENT_LIST_DIR}/${${libname}}" )
    elseif ( ${isCanonical} AND EXISTS "${${libname}}" ) # is a file
        list ( APPEND ${libraries} "${${libname}}" )
    else ( )
        if ( NOT ${${libname}} MATCHES "^-.*" )
            get_doc_libname( docname "${${libname}}" )
            find_library( LIB_${${libname}} "${${libname}}" HINTS ${${libdirs}} DOC
                "Location of ${docname} (needed for ${PROJECT_NAME})" )
            list ( APPEND ${libraries} ${LIB_${${libname}}} )
        else ( ) #append linker flag -...
            list ( APPEND ${libraries} "${${libname}}" )
        endif ( NOT ${${libname}} MATCHES "^-.*" )
    endif ( TARGET "${${libname}}" )
endmacro ( add_library_helper )

function ( kaAddLibs )

    get_directory_property ( libdirs KA_LINKDIRS )

    unset ( qualifier )
    foreach ( libname ${ARGV} )
        if ( ${qualifier} ) # a previous iteration found a qualifier
            add_library_helper ( "${qualifier} ${libname}" libdirs libraries )
            unset ( qualifier )
        else ( )
            if ( libname MATCHES "^(debug|optimized|general)$" )
                set( qualifier "${libname}" ) # set for next iteration
            else ( ) # the current argument is not a qualifier
                # -> add the current argument as a library
               add_library_helper ( libname libdirs libraries )
            endif ( libname MATCHES "^(debug|optimized|general)$" )
        endif ( ${qualifier} )
    endforeach ( libname ${ARGV} )

    foreach ( libdir ${libdirs} )
        string( SUBSTRING ${libdir} 0 1 firstChar)
        string( COMPARE EQUAL ${firstChar} "/" isCanonical )
        if( ${isCanonical} )
            link_directories ( ${libdir} )
        else( ${isCanonical} )
            link_directories ( ${CMAKE_CURRENT_LIST_DIR}/${libdir} )
        endif( ${isCanonical} )
    endforeach( libdir ${libdirs} )


    set_property ( DIRECTORY APPEND PROPERTY KA_LIBRARIES ${libraries} )
    set_property ( DIRECTORY PROPERTY KA_LINKDIRS ${CMAKE_CURRENT_LIST_DIR}/${libdirs} )

endfunction ( kaAddLibs )

function ( kaAddLibdirs )
    foreach ( dir ${ARGV} )
        if ( "${dir}" MATCHES "^-L.+" )
            string ( SUBSTRING ${dir} 2 -1 dir )
        endif ( "${dir}" MATCHES "^-L.+" )

        if ( NOT IS_DIRECTORY "${dir}" )
            message ( WARNING "Can not search directory ${dir}" )
        else ( )
            set_property ( DIRECTORY APPEND PROPERTY KA_LINKDIRS ${dir} )
        endif ( NOT IS_DIRECTORY "${dir}" )
    endforeach ( dir ${ARGV} )

endfunction ( kaAddLibdirs )

function ( kaMakeTool tool )

    set ( extra_includedirs "" )
    set ( extra_sources "" )
    set ( extra_cflags "" )
    set ( extra_ldflags "")
    set ( extra_libs "" )

    get_directory_property ( libdirs KA_LINKDIRS )

    get_filename_component( toolname ${tool} NAME_WE )
    if ( EXISTS "${CMAKE_CURRENT_LIST_DIR}/${tool}" ) #Toolname is also a source (e.g. .cpp) filename
        list( APPEND extra_sources ${CMAKE_CURRENT_LIST_DIR}/${tool})
    endif ( EXISTS "${CMAKE_CURRENT_LIST_DIR}/${tool}" )

    foreach ( arg ${ARGN} )
        if ( "${arg}" MATCHES "^-[DU].+" )
            list( APPEND extra_cflags "${arg}" )
        elseif ( "${arg}" MATCHES "^-L.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search link directory ${dir}" )
            else ( )
                list ( APPEND libdirs "${dir}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( "${arg}" MATCHES "^-I.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search include directory ${dir}" )
            else ( )
                list ( APPEND extra_cflags "${arg}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( ${arg} MATCHES "^-l.+" OR ${arg} MATCHES "\\.(a|so|dylib|dll|lib)$" )
            add_library_helper ( arg libdirs extra_libs )
        elseif ( ${arg} MATCHES "^-.+" )
            list ( APPEND extra_ldflags ${arg} )
        else (  )
            list ( APPEND extra_sources ${arg} )
        endif ( "${arg}" MATCHES "^-[DU].+" )
    endforeach ( arg ${ARGN} )

    add_executable( ${toolname} ${sourcefile} ${extra_sources} ${objs_${CMAKE_PROJECT_NAME}} )

    get_directory_property ( libraries KA_LIBRARIES )

    target_link_libraries( ${toolname} ${extra_libs} )
    target_link_libraries( ${toolname} ${libraries} )

    install( TARGETS ${toolname} RUNTIME DESTINATION ${IBT_BIN_PREFIX} )

    append_flags( ${toolname} COMPILE_FLAGS ${extra_cflags} )
    append_flags( ${toolname} LINK_FLAGS ${extra_ldflags} )

endfunction( kaMakeTool )

function ( kaMakeTools )
    foreach( sources ${ARGV} )
        string( REPLACE "," ";" sourcelist ${sources} )
        list( GET sourcelist 0 sourcefile )
        list( REMOVE_AT sourcelist 0 )

        kaMakeTool( ${sourcefile} "${sourcelist}" )

        get_filename_component( toolname ${sourcefile} NAME_WE )

    endforeach( sources )
endfunction ( kaMakeTools )

function ( kaMakeLib lib )

    add_library( ${lib} STATIC ${ARGN} )
    install( TARGETS ${lib}
    ARCHIVE DESTINATION ${IBT_LIB_PREFIX}
    LIBRARY DESTINATION ${IBT_LIB_PREFIX} )

    # Create an extra target lib<name> (to be able to `make lib<name>`)
    add_custom_target ( lib${lib} DEPENDS ${lib} )

    # BEGIN OBSOLETE:
    # create an alias to the lib we're building that can be linked against
    # as well. this is so that the other targets dont link against a potentially
    # old library that was installed earlier.
    # Caller must use ${lib} when referring to this object collection.
    #add_library ( ${lib}-my ALIAS ${lib} )
    # END OBSOLETE

    # To be compatible with the above behaviour:
    set ( ${lib} "${lib}" PARENT_SCOPE )
endfunction ( kaMakeLib )

function ( kaMakeSharedLib lib )

    set ( extra_includedirs "" )
    set ( extra_sources "" )
    set ( extra_cflags "" )
    set ( extra_ldflags "")
    set ( extra_libs "" )

    get_directory_property ( libdirs KA_LINKDIRS )

    foreach ( arg ${ARGN} )
        if ( "${arg}" MATCHES "^-[DU].+" )
            list( APPEND extra_cflags "${arg}" )
        elseif ( "${arg}" MATCHES "^-L.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search link directory ${dir}" )
            else ( )
                list ( APPEND libdirs "${dir}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( "${arg}" MATCHES "^-I.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search include directory ${dir}" )
            else ( )
                list ( APPEND extra_cflags "${arg}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( ${arg} MATCHES "^-l.+" )
            add_library_helper ( arg libdirs extra_libs )
        elseif ( ${arg} MATCHES "^-.+" )
            list ( APPEND extra_ldflags ${arg} )
        else (  )
            list ( APPEND extra_sources ${arg} )
        endif ( "${arg}" MATCHES "^-[DU].+" )
    endforeach ( arg ${ARGN} )

    add_library( ${lib} SHARED ${extra_sources} )
    target_link_libraries( ${lib} ${extra_libs} )
    get_directory_property ( libraries KA_LIBRARIES )
    target_link_libraries( ${lib} ${libraries} )

    set_target_properties(${lib} PROPERTIES CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    install( TARGETS ${lib}
      ARCHIVE DESTINATION ${IBT_LIB_PREFIX}
      LIBRARY DESTINATION ${IBT_LIB_PREFIX} )

    # Create an extra target lib<name> (to be able to `make lib<name>`)
    add_custom_target ( lib${lib} DEPENDS ${lib} )


    # To be compatible with the above behaviour:
    set ( ${lib} "${lib}" PARENT_SCOPE )
endfunction ( kaMakeSharedLib )

function ( kaMakeModuleLib lib )

    set ( extra_includedirs "" )
    set ( extra_sources "" )
    set ( extra_cflags "" )
    set ( extra_ldflags "")
    set ( extra_libs "" )

    get_directory_property ( libdirs KA_LINKDIRS )

    if( EXISTS ${lib}.cpp )
        set( extra_sources ${lib}.cpp)
    endif( EXISTS ${lib}.cpp )

    foreach ( arg ${ARGN} )
        if ( "${arg}" MATCHES "^-[DU].+" )
            list( APPEND extra_cflags "${arg}" )
        elseif ( "${arg}" MATCHES "^-L.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search link directory ${dir}" )
            else ( )
                list ( APPEND libdirs "${dir}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( "${arg}" MATCHES "^-I.+" )
            string ( SUBSTRING "${arg}" 2 -1 dir )
            if ( NOT IS_DIRECTORY "${dir}" )
                message ( WARNING "Can not search include directory ${dir}" )
            else ( )
                list ( APPEND extra_cflags "${arg}" )
            endif ( NOT IS_DIRECTORY "${dir}" )
        elseif ( ${arg} MATCHES "^-l.+" )
            add_library_helper ( arg libdirs extra_libs )
        elseif ( ${arg} MATCHES "^-.+" )
            list ( APPEND extra_ldflags ${arg} )
        else (  )
            list ( APPEND extra_sources ${arg} )
        endif ( "${arg}" MATCHES "^-[DU].+" )
    endforeach ( arg ${ARGN} )

    add_library( ${lib} MODULE ${extra_sources} )
    target_link_libraries( ${lib} ${extra_libs} )
    get_directory_property ( libraries KA_LIBRARIES )
    target_link_libraries( ${lib} ${libraries} )

    set_target_properties(${lib} PROPERTIES CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    set_target_properties(${lib} PROPERTIES BUNDLE TRUE)

    install( TARGETS ${lib} LIBRARY DESTINATION ${IBT_BIN_PREFIX}/bundles )

    # Create an extra target lib<name> (to be able to `make lib<name>`)
    add_custom_target ( lib${lib} DEPENDS ${lib} )

    # To be compatible with the above behaviour:
    set ( ${lib} "${lib}" PARENT_SCOPE )
endfunction ( kaMakeModuleLib )

# This creates one OBJECT target that is linked into every executable target
function ( kaAddObjects )

    #Create a composite name for the object (file1-file2-file3)
    set ( objname "" )
    foreach( fn ${ARGV} )
        get_filename_component ( fn ${fn} NAME_WE )
        list ( APPEND objname ${fn} )
    endforeach( fn ${ARGV} )
    string ( REPLACE ";" "-" objname ${objname} )

    kaMakeObject( ${objname} ${ARGV} )

    # Add this to a list that is linked to every executable
    # TODO does the scpe work with subdirectories?
    set ( objs_${CMAKE_PROJECT_NAME}
        ${objs_${CMAKE_PROJECT_NAME}} ${${objname}} PARENT_SCOPE )

endfunction ( kaAddObjects )

function ( set_target_stdlib target stdlib )
    get_target_property(cflags ${target} COMPILE_FLAGS )
    set_target_properties( ${target} PROPERTIES COMPILE_FLAGS "${cflags} -stdlib=${stdlib}" )
    get_target_property(ldflags ${target} LINK_FLAGS )
    set_target_properties( ${target} PROPERTIES LINK_FLAGS "${ldflags} -stdlib=${stdlib}" )
endfunction ( set_target_stdlib )

function ( append_flags target flagtype )
    get_target_property(temp ${target} ${flagtype} )
    if (NOT ${temp})
        set(temp "")
    endif (NOT ${temp})
    string( REPLACE ";" " " flags "${ARGN}" )
    set_target_properties( ${target} PROPERTIES ${flagtype} "${temp} ${flags}" )
endfunction ( append_flags )

function( JOIN VALUES GLUE OUTPUT )
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction( JOIN )

function( kaBuildLinkLibrariesString output listIn)
  JOIN("${listIn}" ";-l" _TMP_STR)
  set (${output} "-l${_TMP_STR}" PARENT_SCOPE)
endfunction( kaBuildLinkLibrariesString )

