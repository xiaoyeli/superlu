##################################################################################
#
#                    Set defaults for XSDK CMake projects
#
##################################################################################

#
# This module implements standard behavior for XSDK CMake projects.  It
# changes the default behavior of CMake on a few fronts.  The main thing it
# does by default in XSDK mode is to ignore compiler vars in the env if they
# are set using CC, CXX, FC and compiler flags CFLAGS, CXXFLAGS, and FFLAGS.
#
# This module must be included after:
#
#   PROJECT(${PROJECT_NAME}  NONE)
#
# is called but before the compilers are defined and processed using:
#
#   ENABLE_LANGAUGE(<LANG>)
#
# The major downside of this specification is that when USE_XSDK_DEFAULTS=TRUE
# but XSDK_USE_COMPILER_ENV_VARS=FALSE, then the default compilers must be
# searched for here instead of inside of ENABLE_LANGAUGE(<LANG>) .  This is
# because if go into ENABLE_LANGAUGE(<LANG>) then CMake will set the default
# compilers by reading the env vars.  If you don't like the default compilers
# picked by this module, set XSDK_USE_COMPILER_ENV_VARS=TRUE (which is the
# CMake default anyway).
# 

FUNCTION(PRINT_VAR  VAR_NAME)
  MESSAGE("${VAR_NAME} = '${${VAR_NAME}}'")
ENDFUNCTION()

SET(USE_XSDK_DEFAULTS  FALSE  CACHE  BOOL
  "Use XSDK defaults and behavior.")
PRINT_VAR(USE_XSDK_DEFAULTS)

IF (USE_XSDK_DEFAULTS)

  SET(XSDK_USE_COMPILER_ENV_VARS  FALSE  CACHE  BOOL
    "When in XSDK mode, read defaults for compilers and flags from env.")
  PRINT_VAR(XSDK_USE_COMPILER_ENV_VARS)

  # Handle env vars CXX, CXXFLAGS, etc. ...
  
  IF (XSDK_USE_COMPILER_ENV_VARS)

    # Announce using env var CXX
    IF (NOT "$ENV{CXX}" STREQUAL "" AND "${CMAKE_CXX_COMPILER}" STREQUAL "")
      MESSAGE("XSDK: Setting CMAKE_CXX_COMPILER from env var CXX='$ENV{CXX}'!")
      SET(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE FILEPATH
        "Set by default from env var CXX")
    ENDIF()

    # Announce using env var CXXFLAGS
    IF (NOT "$ENV{CXXFLAGS}" STREQUAL "" AND "${CMAKE_CXX_FLAGS}" STREQUAL "")
      MESSAGE("XSDK: Setting CMAKE_CXX_FLAGS from env var CXXFLAGS='$ENV{CXXFLAGS}'!")
      SET(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} " CACHE  STRING
        "Set by default from env var CXXFLAGS")
      # NOTE: CMake adds the space after ${CXXFLAGS} so we duplicate that here!
    ENDIF()

  ELSE()

    # Ignore env var CXX
    IF (NOT "$ENV{CXX}" STREQUAL "" AND "${CMAKE_CXX_COMPILER}" STREQUAL "")
      MESSAGE("NOT setting CMAKE_CXX_COMPILER from env var CXX='$ENV{CXX}'!")
      # Got to find the default C++ compiler so ENABLE_LANGAUGE(CXX) does not
      # pick up the CXX set in the env!  No way to avoid this that I can see.
      FIND_PROGRAM(DEFAULT_CXX_COMPILER  NAMES  c++  g++)
      SET(CMAKE_CXX_COMPILER "${DEFAULT_CXX_COMPILER}" CACHE FILEPATH
        "Ignoring default set by env var CXX")
    ENDIF()

    # Ignore env var CXXFLAGS
    IF (NOT "$ENV{CXXFLAGS}" STREQUAL "" AND "${CMAKE_CXX_FLAGS}" STREQUAL "")
      MESSAGE("NOT setting CMAKE_CXX_FLAGS from env var CXXFLAGS='$ENV{CXXFLAGS}'!")
      SET(CMAKE_CXX_FLAGS "" CACHE  STRING
        "Ignoring default set by env var CXXFLAGS")
    ENDIF()

    # ToDo: Handle C and Fortran

  ENDIF()
  
  # Set XSDK defaults for other CMake variables
  
  IF ("${BUILD_SHARED_LIBS}"  STREQUAL  "")
    MESSAGE("XSDK: Setting default for BUILD_SHARED_LIBS to TRUE!")
    SET(BUILD_SHARED_LIBS  TRUE  CACHE  BOOL  "Set by default in XSDK mode")
  ENDIF()
  
  IF ("${CMAKE_BUILD_TYPE}"  STREQUAL  "")
    MESSAGE("XSDK: Setting default for CMAKE_BUILD_TYPE to DEBUG!")
    SET(CMAKE_BUILD_TYPE  DEBUG  CACHE  STRING  "Set by default in XSDK mode")
  ENDIF()

ENDIF()
