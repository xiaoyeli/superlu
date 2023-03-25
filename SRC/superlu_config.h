/* superlu_config.h.in */

#ifndef SUPERLU_CONFIG_H
#define SUPERLU_CONFIG_H

/* Enable parmetis */
/* #undef HAVE_PARMETIS */

/* Enable colamd */
/* #undef HAVE_COLAMD */

/* enable 64bit index mode */
/* #undef XSDK_INDEX_SIZE */

/*
 * Integer type for indexing sparse matrix meta structure
 */
#if (XSDK_INDEX_SIZE == 64)
#include <stdint.h>
#define _LONGINT 1
typedef int64_t int_t;
#define IFMT "%lld"
#else
typedef int int_t; /* default */
#define IFMT "%8d"
#endif

#endif /* SUPERLU_CONFIG_H */

