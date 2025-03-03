/*  -- translated by f2c (version 19940927).
*/

#include "../../SRC/slu_scomplex.h"
#include "../../SRC/slu_sdefs.h"
#include "../../SRC/slu_cdefs.h"

/* Table of constant values */

static singlecomplex c_b1 = {0.f,0.f};
static singlecomplex c_b2 = {1.f,0.f};
static int c__3 = 3;
static int c__1 = 1;

/* Subroutine */ int clarge_slu(int *n, singlecomplex *a, int *lda, int *
	iseed, singlecomplex *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;
    double d__1;
    singlecomplex q__1;

    /* Local variables */
    static int i;
    extern /* Subroutine */ int cgerc_(int *, int *, singlecomplex *,
	    singlecomplex *, int *, singlecomplex *, int *, singlecomplex *, int *),
	cscal_(int *, singlecomplex *, singlecomplex *, int *);
    extern float scnrm2_(int *, singlecomplex *, int *);
    static singlecomplex wa, wb;
    static float wn;
    extern /* Subroutine */ int clarnv_slu(int *, int *, int *, singlecomplex *);
    extern int input_error(char *, int *);
    static singlecomplex tau;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARGE pre- and post-multiplies a complex general n by n matrix A   
    with a random unitary matrix: A = U*D*U'.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the original n by n matrix A.   
            On exit, A is overwritten by U*A*U' for some random   
            unitary matrix U.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= N.   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

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
	input_error("CLARGE", &i__1);
	return 0;
    }

/*     pre- and post-multiply A by random unitary matrix */

    for (i = *n; i >= 1; --i) {

/*        generate random reflection */

	i__1 = *n - i + 1;
	clarnv_slu(&c__3, &iseed[1], &i__1, &work[1]);
	i__1 = *n - i + 1;
	wn = scnrm2_(&i__1, &work[1], &c__1);
	d__1 = wn / c_abs(&work[1]);
	q__1.r = d__1 * work[1].r, q__1.i = d__1 * work[1].i;
	wa.r = q__1.r, wa.i = q__1.i;
	if (wn == 0.f) {
	    tau.r = 0.f, tau.i = 0.f;
	} else {
	    q__1.r = work[1].r + wa.r, q__1.i = work[1].i + wa.i;
	    wb.r = q__1.r, wb.i = q__1.i;
	    i__1 = *n - i;
	    c_div(&q__1, &c_b2, &wb);
	    cscal_(&i__1, &q__1, &work[2], &c__1);
	    work[1].r = 1.f, work[1].i = 0.f;
	    c_div(&q__1, &wb, &wa);
	    d__1 = q__1.r;
	    tau.r = d__1, tau.i = 0.f;
	}

/*        multiply A(i:n,1:n) by random reflection from the left */

	i__1 = *n - i + 1;
	cgemv_("Conjugate transpose", &i__1, n, &c_b2, &a[i + a_dim1], lda, &
		work[1], &c__1, &c_b1, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	q__1.r = -(double)tau.r, q__1.i = -(double)tau.i;
	cgerc_(&i__1, n, &q__1, &work[1], &c__1, &work[*n + 1], &c__1, &a[i + 
		a_dim1], lda);

/*        multiply A(1:n,i:n) by random reflection from the right */

	i__1 = *n - i + 1;
	cgemv_("No transpose", n, &i__1, &c_b2, &a[i * a_dim1 + 1], lda, &
		work[1], &c__1, &c_b1, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	q__1.r = -(double)tau.r, q__1.i = -(double)tau.i;
	cgerc_(n, &i__1, &q__1, &work[*n + 1], &c__1, &work[1], &c__1, &a[i * 
		a_dim1 + 1], lda);
/* L10: */
    }
    return 0;

/*     End of CLARGE */

} /* clarge_slu */

