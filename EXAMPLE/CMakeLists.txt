include_directories(${SuperLU_SOURCE_DIR}/SRC)

# build examples for target "all" only if enable_examples is set
if(NOT enable_examples)
  set(_DEPENDENCY_ALL "EXCLUDE_FROM_ALL")
endif()

add_executable(superlu_example
               ${_DEPENDENCY_ALL}
               superlu.c)
target_link_libraries(superlu_example superlu)

# examples for complex
add_executable(citersol
               ${_DEPENDENCY_ALL}
               citersol.c cfgmr.c)
target_link_libraries(citersol superlu)

add_executable(citersol1
               ${_DEPENDENCY_ALL}
               citersol1.c cfgmr.c)
target_link_libraries(citersol1 superlu)

add_executable(clinsol
               ${_DEPENDENCY_ALL}
               clinsol.c)
target_link_libraries(clinsol superlu)

add_executable(clinsol1
               ${_DEPENDENCY_ALL}
               clinsol1.c)
target_link_libraries(clinsol1 superlu)

add_executable(clinsolx
               ${_DEPENDENCY_ALL}
               clinsolx.c)
target_link_libraries(clinsolx superlu)

add_executable(clinsolx1
               ${_DEPENDENCY_ALL}
               clinsolx1.c)
target_link_libraries(clinsolx1 superlu)

add_executable(clinsolx2
               ${_DEPENDENCY_ALL}
               clinsolx2.c)
target_link_libraries(clinsolx2 superlu)

add_executable(clinsolx3
               ${_DEPENDENCY_ALL}
               clinsolx3.c)
target_link_libraries(clinsolx3 superlu)

# examples for double
add_executable(ditersol
               ${_DEPENDENCY_ALL}
               ditersol.c dfgmr.c)
target_link_libraries(ditersol superlu)

add_executable(ditersol1
               ${_DEPENDENCY_ALL}
               ditersol1.c dfgmr.c)
target_link_libraries(ditersol1 superlu)

add_executable(dlinsol
               ${_DEPENDENCY_ALL}
               dlinsol.c)
target_link_libraries(dlinsol superlu)

add_executable(dlinsol1
               ${_DEPENDENCY_ALL}
               dlinsol1.c)
target_link_libraries(dlinsol1 superlu)

add_executable(dlinsolx
               ${_DEPENDENCY_ALL}
               dlinsolx.c)
target_link_libraries(dlinsolx superlu)

add_executable(dlinsolx1
               ${_DEPENDENCY_ALL}
               dlinsolx1.c)
target_link_libraries(dlinsolx1 superlu)

add_executable(dlinsolx2
               ${_DEPENDENCY_ALL}
               dlinsolx2.c)
target_link_libraries(dlinsolx2 superlu)

add_executable(dlinsolx3
               ${_DEPENDENCY_ALL}
               dlinsolx3.c)
target_link_libraries(dlinsolx3 superlu)

# examples for float
add_executable(sitersol
               ${_DEPENDENCY_ALL}
               sitersol.c sfgmr.c)
target_link_libraries(sitersol superlu)

add_executable(sitersol1
               ${_DEPENDENCY_ALL}
               sitersol1.c sfgmr.c)
target_link_libraries(sitersol1 superlu)

add_executable(slinsol
               ${_DEPENDENCY_ALL}
               slinsol.c)
target_link_libraries(slinsol superlu)

add_executable(slinsol1
               ${_DEPENDENCY_ALL}
               slinsol1.c)
target_link_libraries(slinsol1 superlu)

add_executable(slinsolx
               ${_DEPENDENCY_ALL}
               slinsolx.c)
target_link_libraries(slinsolx superlu)

add_executable(slinsolx1
               ${_DEPENDENCY_ALL}
               slinsolx1.c)
target_link_libraries(slinsolx1 superlu)

add_executable(slinsolx2
               ${_DEPENDENCY_ALL}
               slinsolx2.c)
target_link_libraries(slinsolx2 superlu)

add_executable(slinsolx3
               ${_DEPENDENCY_ALL}
               slinsolx3.c)
target_link_libraries(slinsolx3 superlu)

# examples for double complex
add_executable(zitersol
               ${_DEPENDENCY_ALL}
               zitersol.c zfgmr.c)
target_link_libraries(zitersol superlu)

add_executable(zitersol1
               ${_DEPENDENCY_ALL}
               zitersol1.c zfgmr.c)
target_link_libraries(zitersol1 superlu)

add_executable(zlinsol
               ${_DEPENDENCY_ALL}
               zlinsol.c)
target_link_libraries(zlinsol superlu)

add_executable(zlinsol1
               ${_DEPENDENCY_ALL}
               zlinsol1.c)
target_link_libraries(zlinsol1 superlu)

add_executable(zlinsolx
               ${_DEPENDENCY_ALL}
               zlinsolx.c)
target_link_libraries(zlinsolx superlu)

add_executable(zlinsolx1
               ${_DEPENDENCY_ALL}
               zlinsolx1.c)
target_link_libraries(zlinsolx1 superlu)

add_executable(zlinsolx2
               ${_DEPENDENCY_ALL}
               zlinsolx2.c)
target_link_libraries(zlinsolx2 superlu)

add_executable(zlinsolx3
               ${_DEPENDENCY_ALL}
               zlinsolx3.c)
target_link_libraries(zlinsolx3 superlu)