######################################################################
#
# CMakeLists.txt for SUPERLU
#
######################################################################

# Required version
cmake_minimum_required(VERSION 2.8.12)

# Project Version	
project(SuperLU C)
set(VERSION_MAJOR "5")
set(VERSION_MINOR "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR})

######################################################################
#
# IDEAS: xSDK standards module
MESSAGE("\nProcess XSDK defaults ...")
INCLUDE("XSDKDefaults.cmake")
######################################################################

######################################################################
#
# Usual initialization stuff
#
######################################################################
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????

set(NOFORTRAN TRUE)
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
set(CMAKE_C_FLAGS "-DPRNTlevel=0 -DAdd_")
set(CMAKE_C_FLAGS_RELEASE "-O3")

######################################################################
#
# Find packages
#
######################################################################

if(NOT enable_blaslib)
  find_package(BLAS)
endif()

if(BLAS_FOUND)
  add_definitions(-DPRNTlevel=0 -DUSE_VENDOR_BLAS)
  set(BLAS_LIB ${BLAS_LIBRARIES})
else()
  add_subdirectory(CBLAS)
  set(BLAS_LIB blas)
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
string (REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
# configure_file(${CMAKE_SOURCE_DIR}/make.inc.in ${CMAKE_BINARY_DIR}/make.inc)
configure_file(${SuperLU_SOURCE_DIR}/make.inc.in ${SuperLU_SOURCE_DIR}/make.inc)
 
