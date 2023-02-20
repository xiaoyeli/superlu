set(headers 
    ${PROJECT_SOURCE_DIR}/CBLAS/f2c.h
    ${PROJECT_SOURCE_DIR}/CBLAS/slu_Cnames.h
)

# set(sources input_error.c)
set(sources "")

if (enable_single)
    list(APPEND sources
      ${PROJECT_SOURCE_DIR}/CBLAS/isamax.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sasum.c
      ${PROJECT_SOURCE_DIR}/CBLAS/saxpy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/scopy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sdot.c
      ${PROJECT_SOURCE_DIR}/CBLAS/snrm2.c
      ${PROJECT_SOURCE_DIR}/CBLAS/srot.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sscal.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sswap.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sgemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ssymv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/strsv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/sger.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ssyr2.c
    )
endif()

if (enable_double)
    list(APPEND sources
      ${PROJECT_SOURCE_DIR}/CBLAS/idamax.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dasum.c
      ${PROJECT_SOURCE_DIR}/CBLAS/daxpy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dcopy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ddot.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dnrm2.c
      ${PROJECT_SOURCE_DIR}/CBLAS/drot.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dscal.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dswap.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dgemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dsymv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dtrsv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dger.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dsyr2.c
    )
endif()

if (enable_complex)
    list(APPEND sources
      ${PROJECT_SOURCE_DIR}/CBLAS/icamax.c
      ${PROJECT_SOURCE_DIR}/CBLAS/scasum.c
      ${PROJECT_SOURCE_DIR}/CBLAS/caxpy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ccopy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/scnrm2.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cscal.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cswap.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cdotc.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cgemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/chemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ctrsv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cgerc.c
      ${PROJECT_SOURCE_DIR}/CBLAS/cher2.c
    )
endif()

if (enable_complex16)
    list(APPEND sources
      ${PROJECT_SOURCE_DIR}/CBLAS/izamax.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dzasum.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zaxpy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zcopy.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dznrm2.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zscal.c
      ${PROJECT_SOURCE_DIR}/CBLAS/dcabs1.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zswap.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zdotc.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zgemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zhemv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/ztrsv.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zgerc.c
      ${PROJECT_SOURCE_DIR}/CBLAS/zher2.c
    )
endif()

set(CBLAS_SOURCES ${sources} CACHE INTERNAL "")
