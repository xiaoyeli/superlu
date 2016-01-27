
if(enable_single)
  set(SLINEXM slinsol.c)

  add_executable(slinsol ${SLINEXM})
  target_link_libraries(slinsol superlu)
endif()
