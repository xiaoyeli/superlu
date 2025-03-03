
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ void caxpy_(integer *n, singlecomplex *ca, singlecomplex *cx, integer *
	incx, singlecomplex *cy, integer *incy)
{


    /* System generated locals */

    real r__1, r__2;
    singlecomplex q__1, q__2;

    /* Builtin functions */
    double r_imag(singlecomplex *);

    /* Local variables */
    integer i, ix, iy;


/*     constant times a vector plus a vector.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]


    if (*n <= 0) {
	return;
    }
    if ((r__1 = ca->r, dabs(r__1)) + (r__2 = r_imag(ca), dabs(r__2)) == 0.f) {
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
	q__2.r = ca->r * CX(ix).r - ca->i * CX(ix).i, q__2.i = ca->r * CX(
		ix).i + ca->i * CX(ix).r;
	q__1.r = CY(iy).r + q__2.r, q__1.i = CY(iy).i + q__2.i;
	CY(iy).r = q__1.r, CY(iy).i = q__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1 */

L20:
    for (i = 1; i <= *n; ++i) {
	q__2.r = ca->r * CX(i).r - ca->i * CX(i).i, q__2.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	q__1.r = CY(i).r + q__2.r, q__1.i = CY(i).i + q__2.i;
	CY(i).r = q__1.r, CY(i).i = q__1.i;
/* L30: */
    }
    return;
} /* caxpy_ */

