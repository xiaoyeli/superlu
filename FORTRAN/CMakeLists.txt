# include the paths for header files
include_directories(${CMAKE_BINARY_DIR}/FORTRAN)
include_directories(${SuperLU_SOURCE_DIR}/SRC)

# build examples for target "all" only if enable_examples is set
if(NOT enable_examples)
  set(_DEPENDENCY_ALL "EXCLUDE_FROM_ALL")
endif()

# targets to build examples
add_dependencies(examples
                 dfexm
                 zfexm
                 test_omp)

set(sources "")
list(APPEND sources c_fortran_dgssv.c c_fortran_zgssv.c)

add_library(superlu_fortran ${sources})
add_library(superlu_fortran-static STATIC ${sources})
# set(targets superlu_fortran)
get_target_property(superlu_version superlu VERSION)
get_target_property(superlu_soversion superlu SOVERSION)
set_target_properties(superlu_fortran PROPERTIES VERSION ${superlu_version})
set_target_properties(superlu_fortran PROPERTIES SOVERSION ${superlu_soversion})
target_link_libraries(superlu_fortran superlu)

set_target_properties(superlu_fortran-static PROPERTIES OUTPUT_NAME superlu_fortran)
set_target_properties(superlu_fortran-static PROPERTIES VERSION ${superlu_version})
set_target_properties(superlu_fortran-static PROPERTIES SOVERSION ${superlu_soversion})
target_link_libraries(superlu_fortran-static superlu)

# depends on FPP defs
add_dependencies(superlu_fortran-static superlu_fortran)

install(TARGETS superlu_fortran
# DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
)
install(TARGETS superlu_fortran-static
# DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN *.mod
  )
install(FILES ${CMAKE_BINARY_DIR}/FORTRAN/superlu_config.fh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Fortran MPI stuff
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

# Libs to be linked with the Fortran codes
set(fortran_link_libs superlu_fortran ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
set(fortran_link_libs superlu_fortran-static ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
set(all_link_libs ${fortran_link_libs} superlu)

#message("!!! in Fortran: MPI_Fortran_LINK_FLAGS='${MPI_Fortran_LINK_FLAGS}'")
#message("!!! in Fortran: all_link_libs='${all_link_libs}'")
#message("!!! in Fortran: cxx_implicit='${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}'")
if (NOT MSVC)
  list(APPEND all_link_libs m)
endif ()

if(enable_examples)
  add_executable(dfexm
               ${_DEPENDENCY_ALL}
               f_main.F90 hbcode1.F90)
  target_link_libraries(dfexm ${all_link_libs})

  add_executable(zfexm
               ${_DEPENDENCY_ALL}
               z_f_main.F90 zhbcode1.F90)
  target_link_libraries(zfexm ${all_link_libs})

  add_executable(test_omp
               ${_DEPENDENCY_ALL}
               test_omp.F hbcode1.F90)
  target_link_libraries(test_omp ${all_link_libs})
endif() # enable_example

# Format superlu_config.fh from superlu_config.h in C
file(READ "${CMAKE_BINARY_DIR}/SRC/superlu_config.h" CONFIG_FILE_H)
# remove C comments /* */
string(REGEX REPLACE "/\\*[^(\\*/)]*\\*/" " " CONFIG_FILE_H_NO_COMMENTS "${CONFIG_FILE_H}")
# remove typedefs
string(REGEX REPLACE "typedef[^;]*;" " " CONFIG_FILE_FH "${CONFIG_FILE_H_NO_COMMENTS}")
file(WRITE "${CMAKE_BINARY_DIR}/FORTRAN/superlu_config.fh" "${CONFIG_FILE_FH}")
