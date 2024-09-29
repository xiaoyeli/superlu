#include "../../SRC/slu_dcomplex.h"

/* Subroutine */ int zlacgv_slu(int *n, doublecomplex *x, int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLACGV conjugates a complex vector of length N.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The length of the vector X.  N >= 0.   

    X       (input/output) COMPLEX*16 array, dimension   
                           (1+(N-1)*abs(INCX))   
            On entry, the vector of length N to be conjugated.   
            On exit, X is overwritten with conjg(X).   

    INCX    (input) INTEGER   
            The spacing between successive elements of X.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    doublecomplex z__1;

    /* Local variables */
    static int ioff, i;


#define X(I) x[(I)-1]


    if (*incx == 1) {
	for (i = 1; i <= *n; ++i) {
	    d_cnjg(&z__1, &X(i));
	    X(i).r = z__1.r, X(i).i = z__1.i;
/* L10: */
	}
    } else {
	ioff = 1;
	if (*incx < 0) {
	    ioff = 1 - (*n - 1) * *incx;
	}
	for (i = 1; i <= *n; ++i) {
	    d_cnjg(&z__1, &X(ioff));
	    X(ioff).r = z__1.r, X(ioff).i = z__1.i;
	    ioff += *incx;
/* L20: */
	}
    }
    return 0;

/*     End of ZLACGV */

} /* zlacgv_slu */

