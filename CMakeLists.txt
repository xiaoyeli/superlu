# Required version
cmake_minimum_required(VERSION 3.12)

# Project Version
project(SuperLU C)
set(VERSION_MAJOR "7")
set(VERSION_MINOR "0")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

# Compatibility to xSDK standard
# (Extreme-scale Scientific Software Development Kit)
#
#MESSAGE("\nProcess XSDK defaults ...")
#SET(USE_XSDK_DEFAULTS_DEFAULT TRUE) # Set to false if desired
#INCLUDE("cmake/XSDKDefaults.cmake")
set(USE_XSDK_DEFAULTS TRUE)

# The XSDK standard does not allow using internally built BLAS
if (NOT "${enable_internal_blaslib}" STREQUAL "")
  if (USE_XSDK_DEFAULTS)
    set(enable_blaslib_xSDK OFF)
  else()
    set(enable_blaslib_xSDK ON)
  endif()
else()
  set(enable_blaslib_xSDK ${enable_internal_blaslib})
endif()

if (NOT "${enable_fortran}" STREQUAL "")
  if (XSDK_ENABLE_Fortran)
    set(enable_fortran_xSDK ON)
  else()
    set(enable_fortran_xSDK OFF)
  endif()
else()
  set(enable_fortran_xSDK ${enable_fortran})
endif()

# set up options
option(enable_internal_blaslib  "Build the CBLAS library" ${enable_blaslib_xSDK})
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable single precision complex library" ON)
option(enable_complex16 "Enable double precision complex (complex16) library" ON)
option(enable_matlabmex "Build the Matlab mex library" OFF)
option(enable_doc       "Add target 'doc' to build Doxygen documentation" OFF)
option(enable_examples  "Build examples" ON)
option(enable_fortran   "Build Fortran interface" ${enable_fortran_xSDK})
option(enable_tests     "Build tests" ON)
option(enable_compatibility_complex "Provide typdef 'complex' for compatibility with older SuperLU version" OFF)

include(CTest)
include(GNUInstallDirs)
include(FeatureSummary)

######################################################################
#
# Usual initialization stuff
#
######################################################################

# set C standard to C99 without GNU extensions, if not specified by the user
if ("${CMAKE_C_STANDARD}" STREQUAL "")
  set(CMAKE_C_STANDARD 99)
endif()
if ("${CMAKE_C_EXTENSIONS}" STREQUAL "")
  set(CMAKE_C_EXTENSIONS FALSE)
endif()

#---- set RPATH handling
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#---- For shared library
if (BUILD_SHARED_LIBS)
  message("-- SuperLU will be built as a shared library.")
  set(PROJECT_NAME_LIB_EXPORT libsuperlu.so)
else()
  message("-- SuperLU will be built as a static library.")
  set(PROJECT_NAME_LIB_EXPORT libsuperlu.a)
endif()

enable_language(C)
if (enable_fortran)
  enable_language(Fortran)
endif()
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")

if(enable_compatibility_complex)
  add_compile_definitions(SUPERLU_TYPEDEF_COMPLEX)
endif()


#-- BLAS
option(TPL_ENABLE_INTERNAL_BLASLIB  "Build the CBLAS library" ${enable_internal_blaslib})
option(TPL_BLAS_LIBRARIES "List of absolute paths to blas libraries [].")
#-- METIS
option(TPL_ENABLE_METISLIB   "Build the METIS library" OFF)
option(TPL_METIS_LIBRARIES "List of absolute paths to METIS link libraries [].")
option(TPL_METIS_INCLUDE_DIRS "List of absolute paths to METIS include directories [].")

# setup required compiler defines and options.
## get_directory_property( DirDefs COMPILE_DEFINITIONS )
#set(CMAKE_C_FLAGS "-DPRNTlevel=0 -DDEBUGlevel=0 ${CMAKE_C_FLAGS}")

set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "")

if(XSDK_INDEX_SIZE EQUAL 64)
    message("-- Using 64 bit integer for index size.")
endif()

if(MSVC)
  add_compile_definitions(NOMINMAX _COMPLEX_DEFINED)
  # suppress compiler warnings that sprintf and similar functions are unsafe (e.g. C4996)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

######################################################################
#
# Find packages
#
######################################################################
#
#--------------------- BLAS ---------------------
if(NOT enable_internal_blaslib)
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
      set(BLAS_LIB_EXPORT ${CMAKE_INSTALL_PREFIX}/CBLAS/libblas.so)
  else()
      set(BLAS_LIB_EXPORT ${CMAKE_INSTALL_PREFIX}/CBLAS/libblas.a)
  endif()
endif()

#--------------------- METIS ---------------------
if (TPL_ENABLE_METISLIB)   ## want to use metis
  if (NOT TPL_METIS_LIBRARIES)
    message(FATAL_ERROR "TPL_METIS_LIBRARIES option should be set for METIS support to be enabled.")
  endif()

  if (NOT TPL_METIS_INCLUDE_DIRS)
    message(FATAL_ERROR "TPL_METIS_INCLUDE_DIRS option be set for METIS support to be enabled.")
  endif()
  foreach(dir ${TPL_METIS_INCLUDE_DIRS})
    if (NOT MINGW AND NOT EXISTS ${dir})
      message(FATAL_ERROR "METIS include directory not found: ${dir}")
    endif()
    set(CMAKE_C_FLAGS "-I${dir} ${CMAKE_C_FLAGS}")
  endforeach()

  message("-- Enabled support for METIS.")
  set(METIS_FOUND TRUE)

  set(METIS_LIB ${TPL_METIS_LIBRARIES})
  # fix up METIS library names
  string (REPLACE ";" " " METIS_LIB_STR "${METIS_LIB}")
  set(METIS_LIB_EXPORT ${METIS_LIB_STR})
else()
  message("-- Will not link with METIS.")
endif()

if(TPL_ENABLE_METISLIB AND NOT METIS_FOUND)
  find_package(METIS)
  if(METIS_FOUND)
    set(METIS_LIB METIS::METIS)
    set(TPL_METIS_INCLUDE_DIRS "")
  endif()
endif()

if(METIS_FOUND)
  set(HAVE_METIS TRUE)
endif()

if(MSVC)
  find_package(WinGetOpt)
endif()


######################################################################
#
# Include directories
#
######################################################################
include_directories(${CMAKE_BINARY_DIR}/SRC) # For superlu_config.h
if (TPL_METIS_INCLUDE_DIRS)
  include_directories(${TPL_METIS_INCLUDE_DIRS})  ## metis
endif ()

# Generate various configure files with proper definitions

# file(WRITE "make.defs" "# can be exposed to users"
#  ${CMAKE_C_COMPILER}  )
# configure_file(${CMAKE_SOURCE_DIR}/make.inc.in ${CMAKE_SOURCE_DIR}/make.inc)
configure_file(${SuperLU_SOURCE_DIR}/make.inc.in ${SuperLU_SOURCE_DIR}/make.inc)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/superlu.pc.in ${CMAKE_CURRENT_BINARY_DIR}/superlu.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/superlu.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

configure_file(${SuperLU_SOURCE_DIR}/SRC/superlu_config.h.in ${SuperLU_BINARY_DIR}/SRC/superlu_config.h)
configure_file(${SuperLU_SOURCE_DIR}/SRC/superlu_config.h.in ${SuperLU_SOURCE_DIR}/SRC/superlu_config.h)

# Following is to configure a header file for FORTRAN code
configure_file(${SuperLU_SOURCE_DIR}/SRC/superlu_config.h.in ${SuperLU_BINARY_DIR}/FORTRAN/superlu_config.h)

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

add_custom_target(examples)
if (enable_fortran)
  add_subdirectory(FORTRAN)
endif()

add_subdirectory(EXAMPLE)

# To work around CMake bug #18708 that got fixed in 3.28.0, for
# older CMake disable finding Doxygen unless explicitly enabled
if (enable_doc OR (CMAKE_VERSION VERSION_GREATER_EQUAL "3.28.0"))
   add_subdirectory(DOC)
endif()
