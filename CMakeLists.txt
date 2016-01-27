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
set(VERSION_MINOR "1")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# IDEAS: xSDK standards module
MESSAGE("\nProcess XSDK defaults ...")
# SET(USE_XSDK_DEFAULTS_DEFAULT TRUE) # Set to false if desired
INCLUDE("XSDKDefaults.cmake")
######################################################################

######################################################################
#
# Usual initialization stuff
#
######################################################################
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????

enable_language(C)
# enable_language(C++)
enable_language (Fortran)
# set(NOFORTRAN TRUE)
set(NOFORTRAN FALSE)
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")

# setup options
option(enable_blaslib   "Build the CBLAS library" OFF)
option(enable_matlabmex "Build the Matlab mex library" OFF)
option(enable_tests     "Build tests" ON)
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)
option(enable_examples  "Build examples" ON)

# setup required compiler defines and options.
## add_definitions(-DPRNTlevel=0 -DAdd_)
## get_directory_property( DirDefs COMPILE_DEFINITIONS )
set(CMAKE_C_FLAGS "-DPRNTlevel=0 -DAdd_ ${CMAKE_C_FLAGS}" CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE "-O3 ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")

######################################################################
#
# Find packages
#
######################################################################

if(NOT enable_blaslib)
  set(TPL_BLAS_LIBRARIES "" CACHE FILEPATH
    "Override of list of absolute path to libs for BLAS.")
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
  message(STATUS "Using TPL_BLAS_LIBRARIES='${TPL_BLAS_LIBRARIES}'")
  add_definitions(-DPRNTlevel=0 -DUSE_VENDOR_BLAS)
  set(BLAS_LIB ${TPL_BLAS_LIBRARIES})
else()
  message(STATUS "Did not find or specify BLAS so configure to build internal CBLAS ...")
  add_subdirectory(CBLAS)
  set(BLAS_LIB blas)
endif()

if (enable_tests AND BLAS_FOUND)

  set(TPL_LAPACK_LIBRARIES "" CACHE FILEPATH
    "Override of list of absolute path to libs for LAPACK.")
  if (TPL_LAPACK_LIBRARIES)
    set(LAPACK_FOUND TRUE)
  else()
    find_package(LAPACK)
    if (LAPACK_FOUND)
      set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" CACHE FILEPATH
        "Set from FindLAPACK.cmake LAPACK_LIBRARIES." FORCE)
    endif()
  endif()

  if(LAPACK_FOUND)
    message(STATUS "Using TPL_LAPACK_LIBRARIES='${TPL_LAPACK_LIBRARIES}'")
    add_definitions(-DPRNTlevel=0 -DUSE_VENDOR_BLAS)
    set(LAPACK_LIB ${TPL_LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "Error, must find or specify LAPACK when enable_testing=TRUE"
      " and not building internal CBLAS!")
  endif()

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
  message(FATAL_ERROR "Documentation build requested but not implemented.")
  #implement doxygen
endif()

if(enable_examples)
  add_subdirectory(EXAMPLE)
endif()


# file(WRITE "make.defs" "# can be exposed to users"
#  ${CMAKE_C_COMPILER}  )
#string (REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
# configure_file(${CMAKE_SOURCE_DIR}/make.inc.in ${CMAKE_BINARY_DIR}/make.inc)
#configure_file(${SuperLU_SOURCE_DIR}/make.inc.in ${SuperLU_SOURCE_DIR}/make.inc)
 
