if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    add_definitions( -D__DARWIN__ )
    set( KAPREFIX "macosx" )
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions( -Dlinux -D_DEFAULT_SOURCE )
    set( KAPREFIX "linux" )
elseif( ${WIN32})
    add_definitions( -DWIN32 )
    set( KAPREFIX "win32" )
elseif( ${CYGWIN})
    add_definitions( -Dlinux -D_DEFAULT_SOURCE -Di386 )
    set( KAPREFIX "win32" )
endif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
