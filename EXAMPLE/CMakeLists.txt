
if(enable_single)
  set(SLINEXM slinsol.c)

  add_executable(slinsol ${SLINEXM})
  target_link_libraries(slinsol superlu)
  # RAB: You might consider at least running this with add_test() to make sure
  # that it does not segfalt.  Even better, you should define some type of
  # pass/fail criteria to make sure this example is not broken.
endif()
