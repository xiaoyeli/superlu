include_directories(${SuperLU_SOURCE_DIR}/SRC)

add_executable(slinsol
               EXCLUDE_FROM_ALL
               slinsol.c)
target_link_libraries(slinsol superlu)
add_test(NAME slinsol
         COMMAND slinsol "${CMAKE_CURRENT_SOURCE_DIR}/g20.rua")
