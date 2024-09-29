#include "../../SRC/slu_ddefs.h"
#include "powi.h"

#include <stdbool.h>
#include <math.h>

/* Subroutine */ int dlartg_slu(double *f, double *g, double *cs, double *sn, double *r)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine DROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in DBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) DOUBLE PRECISION   
            The first component of vector to be rotated.   

    G       (input) DOUBLE PRECISION   
            The second component of vector to be rotated.   

    CS      (output) DOUBLE PRECISION   
            The cosine of the rotation.   

    SN      (output) DOUBLE PRECISION   
            The sine of the rotation.   

    R       (output) DOUBLE PRECISION   
            The nonzero component of the rotated vector.   

    ===================================================================== 
*/
    /* Initialized data */
    static bool first = true;
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    static int i;
    static double scale;
    static int count;
    static double f1, g1, safmn2, safmx2;
    extern double dmach(char *);
    static double safmin, eps;



    if (first) {
	first = false;
	safmin = dmach("S");
	eps = dmach("E");
	d__1 = dmach("B");
	i__1 = (int) (log(safmin / eps) / log(dmach("B")) / 2.);
	safmn2 = powi(d__1, i__1);
	safmx2 = 1. / safmn2;
    }
    if (*g == 0.) {
	*cs = 1.;
	*sn = 0.;
	*r = *f;
    } else if (*f == 0.) {
	*cs = 0.;
	*sn = 1.;
	*r = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	d__1 = fabs(f1), d__2 = fabs(g1);
	scale = SUPERLU_MAX(d__1,d__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    d__1 = fabs(f1), d__2 = fabs(g1);
	    scale = SUPERLU_MAX(d__1,d__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    d__1 = fabs(f1), d__2 = fabs(g1);
	    scale = SUPERLU_MAX(d__1,d__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	}
	if (fabs(*f) > fabs(*g) && *cs < 0.) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r = -(*r);
	}
    }
    return 0;

/*     End of DLARTG */

} /* dlartg_slu */

