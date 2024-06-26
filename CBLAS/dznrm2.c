
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal dznrm2_(integer *n, doublecomplex *x, integer *incx)
{


    /* System generated locals */

    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    doublereal temp, norm, scale;
    integer ix;
    doublereal ssq;


/*  DZNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DZNRM2 := sqrt( conjg( x' )*x )   


    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to ZLASSQ.   
       Sven Hammarling, Nag Ltd.   
    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

	//    if (*n < 1 || *incx < 1) {
    if (*n < 1) {
	norm = 0.;
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL ZLASSQ( N, X, INCX, SCALE, SSQ ) */
	
	int ixinitial = 1;
	if (*incx < 0) {
	    ixinitial = 1 - (*n-1)* (*incx);
	}

	int i;
	for (i = 1,ix = ixinitial; i <= *n; i++, ix += *incx) {
	    //for (ix = 1; *incx < 0 ? ix >= (*n-1)*(*incx)+1 : ix <= (*n-1)*(*incx)+1; ix += *incx) {
	    
	    if (X(ix).r != 0.) {
		temp = (d__1 = X(ix).r, abs(d__1));
		if (scale < temp) {
/* Computing 2nd power */
		    d__1 = scale / temp;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    d__1 = temp / scale;
		    ssq += d__1 * d__1;
		}
	    }
	    if (d_imag(&X(ix)) != 0.) {
		temp = (d__1 = d_imag(&X(ix)), abs(d__1));
		if (scale < temp) {
/* Computing 2nd power */
		    d__1 = scale / temp;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    d__1 = temp / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DZNRM2. */

} /* dznrm2_ */

