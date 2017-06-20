# use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
find_package( PkgConfig )
if( PKG_CONFIG_FOUND )
    pkg_check_modules( Matheval libmatheval QUIET )
endif()

set( Matheval_DEFINITIONS ${Matheval_CFLAGS_OTHER} )

find_path( Matheval_INCLUDE_DIRS
        NAMES matheval.h
        HINTS
        ${Matheval_INCLUDEDIR}
        ${Matheval_INCLUDE_DIRS}
        )

find_library( Matheval_LIBRARIES
        NAMES matheval
        HINTS
        ${Matheval_LIBDIR}
        ${Matheval_LIBRARY_DIRS}
        )

# handle the QUIETLY and REQUIRED arguments and set Matheval_FOUND to TRUE if
# all listed variables are TRUE
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( Matheval DEFAULT_MSG
        Matheval_LIBRARIES Matheval_INCLUDE_DIRS)

mark_as_advanced(Matheval_INCLUDE_DIRS Matheval_LIBRARIES)