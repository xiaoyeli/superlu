
# include the paths for header files
include_directories(${SuperLU_SOURCE_DIR}/SRC)
include_directories(${SuperLU_BINARY_DIR}/FORTRAN)

set(sources "")

if(enable_double)
  list(APPEND sources c_fortran_dgssv.c hbcode1.F90 f_main.F90)
endif()
if(enable_complex16)
  list(APPEND sources c_fortran_zgssv.c zhbcode1.F90 z_f_main.F90)
endif()

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
add_dependencies(superlu_fortran config_f)
add_dependencies(superlu_fortran-static config_f)
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
install(FILES superlu_config.fh
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

if(enable_double)
## if(FALSE)
  set(F_DEXM f_pddrive.F90)
  add_executable(f_pddrive ${F_DEXM})
  target_link_libraries(f_pddrive ${all_link_libs})
# set_target_properties(f_pddrive PROPERTIES LINKER_LANGUAGE Fortran CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  set_target_properties(f_pddrive PROPERTIES LINKER_LANGUAGE Fortran)
## CXX linker does not work with Intel compiler
# set_target_properties(f_pddrive PROPERTIES LINKER_LANGUAGE CXX LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

  set(F_DEXM3D f_pddrive3d.F90)
  add_executable(f_pddrive3d ${F_DEXM3D})
  target_link_libraries(f_pddrive3d ${all_link_libs})
  set_target_properties(f_pddrive3d PROPERTIES LINKER_LANGUAGE Fortran)

  set(F_5x5 f_5x5.F90 sp_ienv.c)
  add_executable(f_5x5 ${F_5x5})
  target_link_libraries(f_5x5 ${all_link_libs})
  set_target_properties(f_5x5 PROPERTIES LINKER_LANGUAGE Fortran)

endif() ## enable_double

##if(FALSE)
if(enable_complex16)
  set(F_ZEXM f_pzdrive.F90)
  add_executable(f_pzdrive ${F_ZEXM})
  target_link_libraries(f_pzdrive ${all_link_libs})
  set_target_properties(f_pzdrive PROPERTIES LINKER_LANGUAGE Fortran)
#  set_target_properties(f_pzdrive PROPERTIES LINKER_LANGUAGE CXX LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

  set(F_ZEXM3D f_pzdrive3d.F90)
  add_executable(f_pzdrive3d ${F_ZEXM3D})
  target_link_libraries(f_pzdrive3d ${all_link_libs})
  set_target_properties(f_pzdrive3d PROPERTIES LINKER_LANGUAGE Fortran)
#  set_target_properties(f_pzdrive3d PROPERTIES LINKER_LANGUAGE CXX LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif() ## enable_complex16

endif ()  ## enable_example

# Format superlu_config.fh from superlu_config.h in C
add_custom_command(
  OUTPUT superlu_config.fh
  COMMAND sed;'/^\\//;d';<;superlu_config.h;>;temp.fh
  COMMAND sed;'/typedef/;d';<;temp.fh;>;superlu_config.fh
  COMMAND cp;superlu_config.fh;${SuperLU_SOURCE_DIR}/FORTRAN/.
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/FORTRAN
)
add_custom_target(config_f DEPENDS superlu_config.fh)
