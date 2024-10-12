include_directories(${SuperLU_SOURCE_DIR}/SRC)

# Libs linked to all of the tests
set(test_link_libs matgen)

add_subdirectory(MATGEN)

set(NVAL   9 19)
set(NRHS   2)
set(LWORK  0 10000000)

function(add_superlu_test target input)
  set(TEST_INPUT "${SuperLU_SOURCE_DIR}/EXAMPLE/${input}")

  foreach (s ${NRHS})
    foreach(l ${LWORK})
      # LA tests
      foreach (n ${NVAL})
        set(testName "${target}_${n}_${s}_${l}")
        add_test(
          NAME ${testName}_LA
          COMMAND ${target} -t "LA" -n ${n} -s ${s} -l ${l})
      endforeach()

      # SP tests
      set(testName "${target}_${s}_${l}")
      add_test(
        NAME ${testName}_SP
        COMMAND ${target} -t "SP" -s ${s} -l ${l} -f "${TEST_INPUT}")
    endforeach()
  endforeach()

  # Add getopt dependency to target in case of MSVC
  if(MSVC)
    target_include_directories(${target} PRIVATE ${WinGetOpt_INCLUDE_DIR})
    target_link_libraries(${target} ${WinGetOpt_LIBRARY})
  endif()

endfunction(add_superlu_test)

# prevent creating testing for Windows if WinGetOps is not found
if(MSVC AND NOT WinGetOpt_FOUND)
  message("Disabled tests s_test, d_test, c_test, and z_test because optional dependency WinGetOpt is missing")
  return()
endif()

if(enable_single)
  add_executable(s_test
                 sdrive.c sgst01.c sgst02.c sgst04.c sgst07.c
                 sp_ienv.c sp_sconvert.c)
  target_link_libraries(s_test ${test_link_libs})
  add_superlu_test(s_test g20.rua)
endif()


if(enable_double)
  add_executable(d_test
                 ddrive.c dgst01.c dgst02.c dgst04.c dgst07.c
                 sp_ienv.c sp_dconvert.c)
  target_link_libraries(d_test ${test_link_libs})
  add_superlu_test(d_test g20.rua)
endif()

if(enable_complex)
  add_executable(c_test
                 cdrive.c cgst01.c cgst02.c cgst04.c cgst07.c
                 sp_ienv.c sp_cconvert.c)
  target_link_libraries(c_test ${test_link_libs})
  add_superlu_test(c_test cg20.cua)
endif()


if(enable_complex16)
  add_executable(z_test
                 zdrive.c zgst01.c zgst02.c zgst04.c zgst07.c
                 sp_ienv.c sp_zconvert.c)
  target_link_libraries(z_test ${test_link_libs})
  add_superlu_test(z_test cg20.cua)
endif()
