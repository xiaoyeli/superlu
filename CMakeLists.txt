######################################################################
#
# CMakeLists.txt for SUPERLU
#
######################################################################

# Required version
cmake_minimum_required(VERSION 2.8.12)

# Project Version
project(SuperLU NONE)
set(VERSION_MAJOR "5")
set(VERSION_MINOR "2")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# IDEAS: xSDK standards module
MESSAGE("\nProcess XSDK defaults ...")
# SET(USE_XSDK_DEFAULTS_DEFAULT TRUE) # Set to false if desired
INCLUDE("cmake/XSDKDefaults.cmake")
INCLUDE(CTest)

######################################################################

######################################################################
#
# Usual initialization stuff
#
######################################################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

#---- For shared library

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#----

enable_language(C)
message("CMAKE_SHARED_LIBRARY_PREFIX - ${CMAKE_SHARED_LIBRARY_PREFIX}")
message("CMAKE_SHARED_LIBRARY_SUFFIX - ${CMAKE_SHARED_LIBRARY_SUFFIX}")
message("CMAKE_IMPORT_LIBRARY_PREFIX - ${CMAKE_IMPORT_LIBRARY_PREFIX}")
message("CMAKE_IMPORT_LIBRARY_SUFFIX - ${CMAKE_IMPORT_LIBRARY_SUFFIX}")
message("CMAKE_STATIC_LIBRARY_PREFIX - ${CMAKE_STATIC_LIBRARY_PREFIX}")
message("CMAKE_STATIC_LIBRARY_SUFFIX - ${CMAKE_STATIC_LIBRARY_SUFFIX}")

if (BUILD_SHARED_LIBS)
  message("-- SuperLU will be built as a shared library.")
  if (WIN32 OR CYGWIN OR MSYS)
    # In Windows, you link to an associated import library, not directly to the shared library runtime (.dll)
    set(PROJECT_NAME_LIB_EXPORT ${CMAKE_IMPORT_LIBRARY_PREFIX}superlu${CMAKE_IMPORT_LIBRARY_SUFFIX})
  else()
    set(PROJECT_NAME_LIB_EXPORT ${CMAKE_SHARED_LIBRARY_PREFIX}superlu${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()
else()
  message("-- SuperLU will be built as a static library.")
  set(PROJECT_NAME_LIB_EXPORT ${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()


if (XSDK_ENABLE_Fortran)
  enable_language (Fortran)
  set(NOFORTRAN FALSE)
endif()
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")

# The XSDK standard does not allow using internally built BLAS
if (USE_XSDK_DEFAULTS)
  set(enable_blaslib_DEFAULT OFF)
else()
  set(enable_blaslib_DEFAULT ON)
endif()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

# setup options
option(enable_blaslib   "Build the CBLAS library (requires f2c)" ${enable_blaslib_DEFAULT})
option(enable_matlabmex "Build the Matlab mex library" OFF)
option(enable_tests     "Build tests" ON)
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_doc_html  "Build HTML documentation"    OFF)
option(enable_doc_man   "Build Unix man pages"        OFF)
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)
# option(enable_examples  "Build examples" ON)

# setup required compiler defines and options.
## add_definitions(-DPRNTlevel=0 -DAdd_)
## get_directory_property( DirDefs COMPILE_DEFINITIONS )
set(CMAKE_C_FLAGS "-DPRNTlevel=0 -DDEBUGlevel=0 ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "")

######################################################################
#
# Find packages
#
######################################################################
#
#--------------------- BLAS ---------------------
if(NOT enable_blaslib)
  if (TPL_BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
  else()
    find_package(BLAS)
    if (BLAS_FOUND)
      set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
        "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
    endif()
  endif()
endif()

if(BLAS_FOUND)
  message("-- Using TPL_BLAS_LIBRARIES='${TPL_BLAS_LIBRARIES}'")
  set(CMAKE_C_FLAGS "-DUSE_VENDOR_BLAS ${CMAKE_C_FLAGS}")
  set(BLAS_LIB ${TPL_BLAS_LIBRARIES})
  # fix up BLAS library name
  string (REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
  set(BLAS_LIB_EXPORT ${BLAS_LIB_STR})
else()
  message("-- Did not find or specify BLAS so configure to build internal CBLAS ...")
  add_subdirectory(CBLAS)
  set(BLAS_LIB blas)
  if (BUILD_SHARED_LIBS)  # export to be referenced by downstream makefile
  #This must always be a static library.
      set(BLAS_LIB_EXPORT ${CMAKE_INSTALL_PREFIX}/CBLAS/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
endif()

######################################################################
#
# TMGLIB_EXPORT value
#
######################################################################
  if (BUILD_SHARED_LIBS)  # export to be referenced by downstream makefile
      if (WIN32 OR CYGWIN OR MSYS)
      # In Windows, you link to an associated import library, not directly to the shared library runtime (.dll)
         set(TMGLIB_EXPORT ${CMAKE_IMPORT_LIBRARY_PREFIX}tmglib${CMAKE_IMPORT_LIBRARY_SUFFIX})
      else()
         set(TMGLIB_EXPORT ${CMAKE_SHARED_LIBRARY_PREFIX}tmglib${CMAKE_SHARED_LIBRARY_SUFFIX})
      endif()
  else()
      set(TMGLIB_EXPORT ${CMAKE_STATIC_LIBRARY_PREFIX}tmglib${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

######################################################################
#
# Include directories
#
######################################################################

include_directories(${CMAKE_SOURCE_DIR}/SRC)

######################################################################
#
# Add subdirectories
#
######################################################################

add_subdirectory(SRC)

if(enable_matlabmex)
  add_subdirectory(MATLAB)
endif()

if(enable_tests)
  enable_testing()
  add_subdirectory(TESTING)
endif()

if(enable_doc)
#  message(FATAL_ERROR "Documentation build requested but not implemented.")
  #implement doxygen

# check if Doxygen is installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    include(GNUInstallDirs)
    message(STATUS "doxygen                           ${DOXYGEN_EXECUTABLE}")
    message(STATUS "doxygen version                   ${DOXYGEN_VERSION}")
    #based from: https://vicrucann.github.io/tutorials/quick-cmake-doxygen/
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    if(DOXYGEN_DOT_FOUND)
      set(HAVE_DOT "YES")
      set(DOT_EXECUTABLE ${DOXYGEN_DOT_EXECUTABLE})
      message(STATUS "Graphviz    (dot)                 ${DOXYGEN_DOT_EXECUTABLE}")
    else()
      set(HAVE_DOT "NO")
      unset(DOT_EXECUTABLE)
      message(STATUS "Graphviz    (dot)                 Not Found")
    endif(DOXYGEN_DOT_FOUND)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Optional - use mscgen
    find_program(MSCGEN mscgen)
    if(MSCGEN)
       message(STATUS "mscgen                            ${MSCGEN}")
       # Doxygen will use this dir to find mscgen
       get_filename_component(MSCGEN_PATH ${MSCGEN} DIRECTORY)
    else()
       message(STATUS "mscgen                            Not Found")
       unset(MSCGEN_PATH)
    endif(MSCGEN)
    # Optional - use dia
    find_program(DIA dia)
    if(DIA)
       message(STATUS "dia                               ${DIA}")
       # Doxygen will use this dir to find mscgen
       get_filename_component(DIA_PATH ${DIA} DIRECTORY)
    else()
       message(STATUS "dia                               Not Found")
       unset(DIA_PATH)
    endif(DIA)
    # note the option ALL which allows to build the docs together with the application
    # We want to set up the docs in HTML format and the option for it in these conditionals
    if(enable_doc_html)
	   set(d_enable_doc_html "YES")
           install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
    else()
           set(d_enable_doc_html "NO")
    endif(enable_doc_html)
    if(enable_doc_man)
	   set(d_enable_doc_man "YES")
	   install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/man/ DESTINATION ${CMAKE_INSTALL_MANDIR})
    else()
	   set(d_enable_doc_man "NO")
    endif(enable_doc_man)
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
  else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()

# file(WRITE "make.defs" "# can be exposed to users"
#  ${CMAKE_C_COMPILER}  )
# configure_file(${CMAKE_SOURCE_DIR}/make.inc.in ${CMAKE_SOURCE_DIR}/make.inc)
configure_file(${SuperLU_SOURCE_DIR}/make.inc.in ${SuperLU_SOURCE_DIR}/make.inc)
# This is just here in case people compile the same code for different sitatuions.
file(COPY ${SuperLU_SOURCE_DIR}/make.inc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/superlu.pc.in ${CMAKE_CURRENT_BINARY_DIR}/superlu.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/superlu.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
