
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cscal_(integer *n, singlecomplex *ca, singlecomplex *cx, integer *
	incx)
{


    /* System generated locals */

    singlecomplex q__1;

    /* Local variables */
    integer i, nincx;


/*     scales a vector by a constant.   
       jack dongarra, linpack,  3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	q__1.r = ca->r * CX(i).r - ca->i * CX(i).i, q__1.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	CX(i).r = q__1.r, CX(i).i = q__1.i;
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */

L20:
    for (i = 1; i <= *n; ++i) {
	q__1.r = ca->r * CX(i).r - ca->i * CX(i).i, q__1.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	CX(i).r = q__1.r, CX(i).i = q__1.i;
/* L30: */
    }
    return 0;
} /* cscal_ */

