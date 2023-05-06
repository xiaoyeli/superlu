/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#include "slu_Cnames.h"

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#if 0
typedef long int integer; /* 64 on 64-bit machine */
#endif

typedef int integer;

typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#endif
