include_directories(${SuperLU_SOURCE_DIR}/SRC)

add_executable(superlu_example
               superlu.c)
target_link_libraries(superlu_example superlu)

# examples for complex
add_executable(citersol
               citersol.c cfgmr.c)
target_link_libraries(citersol superlu)

add_executable(citersol1
               citersol1.c cfgmr.c)
target_link_libraries(citersol1 superlu)

add_executable(clinsol
               clinsol.c)
target_link_libraries(clinsol superlu)

add_executable(clinsol1
               clinsol1.c)
target_link_libraries(clinsol1 superlu)

add_executable(clinsolx
               clinsolx.c)
target_link_libraries(clinsolx superlu)

add_executable(clinsolx1
               clinsolx1.c)
target_link_libraries(clinsolx1 superlu)

add_executable(clinsolx2
               clinsolx2.c)
target_link_libraries(clinsolx2 superlu)

add_executable(clinsolx3
               clinsolx3.c)
target_link_libraries(clinsolx3 superlu)

# examples for double
add_executable(ditersol
               ditersol.c dfgmr.c)
target_link_libraries(ditersol superlu)

add_executable(ditersol1
               ditersol1.c dfgmr.c)
target_link_libraries(ditersol1 superlu)

add_executable(dlinsol
               dlinsol.c)
target_link_libraries(dlinsol superlu)
add_executable(dlinsol1
               dlinsol1.c)
target_link_libraries(dlinsol1 superlu)

add_executable(dlinsolx
               dlinsolx.c)
target_link_libraries(dlinsolx superlu)

add_executable(dlinsolx1
               dlinsolx1.c)
target_link_libraries(dlinsolx1 superlu)

add_executable(dlinsolx2
               dlinsolx2.c)
target_link_libraries(dlinsolx2 superlu)

add_executable(dlinsolx3
               dlinsolx3.c)
target_link_libraries(dlinsolx3 superlu)

# examples for float
add_executable(sitersol
               sitersol.c sfgmr.c)
target_link_libraries(sitersol superlu)

add_executable(sitersol1
               sitersol1.c sfgmr.c)
target_link_libraries(sitersol1 superlu)

add_executable(slinsol
               slinsol.c)
target_link_libraries(slinsol superlu)

add_executable(slinsol1
               slinsol1.c)
target_link_libraries(slinsol1 superlu)
add_executable(slinsolx
               slinsolx.c)
target_link_libraries(slinsolx superlu)

add_executable(slinsolx1
               slinsolx1.c)
target_link_libraries(slinsolx1 superlu)

add_executable(slinsolx2
               slinsolx2.c)
target_link_libraries(slinsolx2 superlu)

add_executable(slinsolx3
               slinsolx3.c)
target_link_libraries(slinsolx3 superlu)

# examples for double complex
add_executable(zitersol
               zitersol.c zfgmr.c)
target_link_libraries(zitersol superlu)

add_executable(zitersol1
               zitersol1.c zfgmr.c)
target_link_libraries(zitersol1 superlu)

add_executable(zlinsol
               zlinsol.c)
target_link_libraries(zlinsol superlu)

add_executable(zlinsol1
               zlinsol1.c)
target_link_libraries(zlinsol1 superlu)
add_executable(zlinsolx
               zlinsolx.c)
target_link_libraries(zlinsolx superlu)

add_executable(zlinsolx1
               zlinsolx1.c)
target_link_libraries(zlinsolx1 superlu)

add_executable(zlinsolx2
               zlinsolx2.c)
target_link_libraries(zlinsolx2 superlu)

add_executable(zlinsolx3
               zlinsolx3.c)
target_link_libraries(zlinsolx3 superlu)
