include(FindPackageHandleStandardArgs)

if (NOT WIN32)
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
         pkg_check_modules(PKG_LIBGETOPT libgetopt)
    endif ()
endif (NOT WIN32)

find_path(GETOPT_INCLUDE_DIR getopt.h
    ${PKG_LIBGETOPT_INCLUDE_DIRS}
    /usr/include
    /usr/local/include
)

find_library(GETOPT_LIBRARY
    NAMES
    getopt
    PATHS
    ${PKG_LIBGETOPT_LIBRARY_DIRS}
    /usr/lib
    /usr/local/lib
)

find_package_handle_standard_args(GetOpt DEFAULT_MSG GETOPT_LIBRARY GETOPT_INCLUDE_DIR)
