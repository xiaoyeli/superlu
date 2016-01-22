
if(enable_single)
  set(SLINEXM slinsol.c)

  add_executable(slinsol ${SLINEXM})
  target_link_libraries(slinsol superlu_${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix} ${BLAS_LIB} m)
endif()
