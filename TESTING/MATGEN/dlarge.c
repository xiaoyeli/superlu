/*  -- translated by f2c (version 19940927).
*/

#include "../../SRC/slu_ddefs.h"

/* Table of constant values */

static int c__3 = 3;
static int c__1 = 1;
static double c_b8 = 1.;
static double c_b10 = 0.;

/* Subroutine */ int dlarge_slu(int *n, double *a, int *lda, int
	*iseed, double *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    double d__1;

    /* Local variables */
    extern /* Subroutine */ int dger_(int *, int *, double *,
	    double *, int *, double *, int *, double *,
	    int *);
    extern double dnrm2_(int *, double *, int *);
    static int i;
    extern /* Subroutine */ int dscal_(int *, double *, double *, int *);
    static double wa, wb, wn;
    extern /* Subroutine */ int dlarnv_slu(int *, int *, int *, double *);
    extern int input_error(char *, int *);
    static double tau;


/*  -- LAPACK auxiliary test routine (version 2.0)   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARGE pre- and post-multiplies a real general n by n matrix A   
    with a random orthogonal matrix: A = U*D*U'.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the original n by n matrix A.   
            On exit, A is overwritten by U*A*U' for some random   
            orthogonal matrix U.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= N.   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*lda < SUPERLU_MAX(1,*n)) {
	*info = -3;
    }
    if (*info < 0) {
	i__1 = -(*info);
	input_error("DLARGE", &i__1);
	return 0;
    }

/*     pre- and post-multiply A by random orthogonal matrix */

    for (i = *n; i >= 1; --i) {

/*        generate random reflection */

	i__1 = *n - i + 1;
	dlarnv_slu(&c__3, &iseed[1], &i__1, &work[1]);
	i__1 = *n - i + 1;
	wn = dnrm2_(&i__1, &work[1], &c__1);
	wa = copysign(wn, work[1]);
	if (wn == 0.) {
	    tau = 0.;
	} else {
	    wb = work[1] + wa;
	    i__1 = *n - i;
	    d__1 = 1. / wb;
	    dscal_(&i__1, &d__1, &work[2], &c__1);
	    work[1] = 1.;
	    tau = wb / wa;
	}

/*        multiply A(i:n,1:n) by random reflection from the left */

	i__1 = *n - i + 1;
	dgemv_("Transpose", &i__1, n, &c_b8, &a[i + a_dim1], lda, &work[1], &
		c__1, &c_b10, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	d__1 = -tau;
	dger_(&i__1, n, &d__1, &work[1], &c__1, &work[*n + 1], &c__1, &a[i + 
		a_dim1], lda);

/*        multiply A(1:n,i:n) by random reflection from the right */

	i__1 = *n - i + 1;
	dgemv_("No transpose", n, &i__1, &c_b8, &a[i * a_dim1 + 1], lda, &
		work[1], &c__1, &c_b10, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	d__1 = -tau;
	dger_(n, &i__1, &d__1, &work[*n + 1], &c__1, &work[1], &c__1, &a[i * 
		a_dim1 + 1], lda);
/* L10: */
    }
    return 0;

/*     End of DLARGE */

} /* dlarge_slu */

