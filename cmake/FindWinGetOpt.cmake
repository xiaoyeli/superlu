include(FindPackageHandleStandardArgs)

find_path(WinGetOpt_INCLUDE_DIR getopt.h
    ${PKG_LIBGETOPT_INCLUDE_DIRS}
    /usr/include
    /usr/local/include
)

find_library(WinGetOpt_LIBRARY
    NAMES
    getopt
    PATHS
    ${PKG_LIBGETOPT_LIBRARY_DIRS}
    /usr/lib
    /usr/local/lib
)

find_package_handle_standard_args(WinGetOpt DEFAULT_MSG WinGetOpt_LIBRARY WinGetOpt_INCLUDE_DIR)
