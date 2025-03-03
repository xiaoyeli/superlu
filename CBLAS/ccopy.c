
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ void ccopy_(integer *n, singlecomplex *cx, integer *incx, singlecomplex *
	cy, integer *incy)
{


    /* System generated locals */


    /* Local variables */
    integer i, ix, iy;


/*     copies a vector, x, to a vector, y.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]


    if (*n <= 0) {
	return;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    for (i = 1; i <= *n; ++i) {
	CY(iy).r = CX(ix).r, CY(iy).i = CX(ix).i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1 */

L20:
    for (i = 1; i <= *n; ++i) {
	CY(i).r = CX(i).r, CY(i).i = CX(i).i;
/* L30: */
    }
    return;
} /* ccopy_ */

