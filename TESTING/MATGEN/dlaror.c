/*  -- translated by f2c (version 19940927).
*/
#include <math.h>
#include <string.h>

/* Table of constant values */

static double c_b9 = 0.;
static double c_b10 = 1.;
static int c__3 = 3;
static int c__1 = 1;

/* Subroutine */ int dlaror_slu(char *side, char *init, int *m, int *n,
	double *a, int *lda, int *iseed, double *x, int *
	info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    double d__1;

    /* Local variables */
    static int kbeg;
    extern /* Subroutine */ int dger_(int *, int *, double *,
	    double *, int *, double *, int *, double *,
	    int *);
    static int jcol, irow;
    extern double dnrm2_(int *, double *, int *);
    static int j;
    extern /* Subroutine */ int dscal_(int *, double *, double *,
	    int *);
    extern /* Subroutine */ int dgemv_(char *, int *, int *,
	    double *, double *, int *, double *, int *,
	    double *, double *, int *);
    static int ixfrm, itype, nxfrm;
    static double xnorm;
    extern double dlarnd_slu(int *, int *);
    extern /* Subroutine */ int dlaset_slu(char *, int *, int *,
					double *, double *, double *, int *);
    extern int input_error(char *, int *);
    static double factor, xnorms;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAROR pre- or post-multiplies an M by N matrix A by a random   
    orthogonal matrix U, overwriting A.  A may optionally be initialized 
  
    to the identity matrix before multiplying by U.  U is generated using 
  
    the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409). 
  

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether A is multiplied on the left or right by U. 
  
            = 'L':         Multiply A on the left (premultiply) by U   
            = 'R':         Multiply A on the right (postmultiply) by U'   
            = 'C' or 'T':  Multiply A on the left by U and the right   
                            by U' (Here, U' means U-transpose.)   

    INIT    (input) CHARACTER*1   
            Specifies whether or not A should be initialized to the   
            identity matrix.   
            = 'I':  Initialize A to (a section of) the identity matrix   
                     before applying U.   
            = 'N':  No initialization.  Apply U to the input matrix A.   

            INIT = 'I' may be used to generate square or rectangular   
            orthogonal matrices:   

            For M = N and SIDE = 'L' or 'R', the rows will be orthogonal 
  
            to each other, as will the columns.   

            If M < N, SIDE = 'R' produces a dense matrix whose rows are   
            orthogonal and whose columns are not, while SIDE = 'L'   
            produces a matrix whose rows are orthogonal, and whose first 
  
            M columns are orthogonal, and whose remaining columns are   
            zero.   

            If M > N, SIDE = 'L' produces a dense matrix whose columns   
            are orthogonal and whose rows are not, while SIDE = 'R'   
            produces a matrix whose columns are orthogonal, and whose   
            first M rows are orthogonal, and whose remaining rows are   
            zero.   

    M       (input) INTEGER   
            The number of rows of A.   

    N       (input) INTEGER   
            The number of columns of A.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the array A.   
            On exit, overwritten by U A ( if SIDE = 'L' ),   
             or by A U ( if SIDE = 'R' ),   
             or by U A U' ( if SIDE = 'C' or 'T').   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry ISEED specifies the seed of the random number   
            generator. The array elements should be between 0 and 4095;   
            if not they will be reduced mod 4096.  Also, ISEED(4) must   
            be odd.  The random number generator uses a linear   
            congruential sequence limited to small integers, and so   
            should produce machine independent random numbers. The   
            values of ISEED are changed on exit, and can be used in the   
            next call to DLAROR to continue the same random number   
            sequence.   

    X       (workspace) DOUBLE PRECISION array, dimension (3*MAX( M, N )) 
  
            Workspace of length   
                2*M + N if SIDE = 'L',   
                2*N + M if SIDE = 'R',   
                3*N     if SIDE = 'C' or 'T'.   

    INFO    (output) INTEGER   
            An error flag.  It is set to:   
            = 0:  normal return   
            < 0:  if INFO = -k, the k-th argument had an illegal value   
            = 1:  if the random numbers generated by DLARND are bad.   

    ===================================================================== 
  


       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --x;

    /* Function Body */
    if (*n == 0 || *m == 0) {
	return 0;
    }

    itype = 0;
    if (strncmp(side, "L", 1)==0) {
	itype = 1;
    } else if (strncmp(side, "R", 1)==0) {
	itype = 2;
    } else if (strncmp(side, "C", 1)==0 || strncmp(side, "T", 1)==0) {
	itype = 3;
    }

/*     Check for argument errors. */

    *info = 0;
    if (itype == 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0 || (itype == 3 && *n != *m)) {
	*info = -4;
    } else if (*lda < *m) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	input_error("DLAROR", &i__1);
	return 0;
    }

    if (itype == 1) {
	nxfrm = *m;
    } else {
	nxfrm = *n;
    }

/*     Initialize A to the identity matrix if desired */

    if (strncmp(init, "I", 1)==0) {
	dlaset_slu("Full", m, n, &c_b9, &c_b10, &a[a_offset], lda);
    }

/*     If no rotation possible, multiply by random +/-1   

       Compute rotation by computing Householder transformations   
       H(2), H(3), ..., H(nhouse) */

    i__1 = nxfrm;
    for (j = 1; j <= i__1; ++j) {
	x[j] = 0.;
/* L10: */
    }

    i__1 = nxfrm;
    for (ixfrm = 2; ixfrm <= i__1; ++ixfrm) {
	kbeg = nxfrm - ixfrm + 1;

/*        Generate independent normal( 0, 1 ) random numbers */

	i__2 = nxfrm;
	for (j = kbeg; j <= i__2; ++j) {
	    x[j] = dlarnd_slu(&c__3, &iseed[1]);
/* L20: */
	}

/*        Generate a Householder transformation from the random vector
 X */

	xnorm = dnrm2_(&ixfrm, &x[kbeg], &c__1);
	xnorms = copysign(xnorm, x[kbeg]);
	d__1 = -x[kbeg];
	x[kbeg + nxfrm] = copysign(c_b10, d__1);
	factor = xnorms * (xnorms + x[kbeg]);
	if (fabs(factor) < 1e-20) {
	    *info = 1;
	    input_error("DLAROR", info);
	    return 0;
	} else {
	    factor = 1. / factor;
	}
	x[kbeg] += xnorms;

/*        Apply Householder transformation to A */

	if (itype == 1 || itype == 3) {

/*           Apply H(k) from the left. */

	    dgemv_("T", &ixfrm, n, &c_b10, &a[kbeg + a_dim1], lda, &x[kbeg], &
		    c__1, &c_b9, &x[(nxfrm << 1) + 1], &c__1);
	    d__1 = -factor;
	    dger_(&ixfrm, n, &d__1, &x[kbeg], &c__1, &x[(nxfrm << 1) + 1], &
		    c__1, &a[kbeg + a_dim1], lda);

	}

	if (itype == 2 || itype == 3) {

/*           Apply H(k) from the right. */

	    dgemv_("N", m, &ixfrm, &c_b10, &a[kbeg * a_dim1 + 1], lda, &x[
		    kbeg], &c__1, &c_b9, &x[(nxfrm << 1) + 1], &c__1);
	    d__1 = -factor;
	    dger_(m, &ixfrm, &d__1, &x[(nxfrm << 1) + 1], &c__1, &x[kbeg], &
		    c__1, &a[kbeg * a_dim1 + 1], lda);

	}
/* L30: */
    }

    d__1 = dlarnd_slu(&c__3, &iseed[1]);
    x[nxfrm * 2] = copysign(c_b10, d__1);

/*     Scale the matrix A by D. */

    if (itype == 1 || itype == 3) {
	i__1 = *m;
	for (irow = 1; irow <= i__1; ++irow) {
	    dscal_(n, &x[nxfrm + irow], &a[irow + a_dim1], lda);
/* L40: */
	}
    }

    if (itype == 2 || itype == 3) {
	i__1 = *n;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    dscal_(m, &x[nxfrm + jcol], &a[jcol * a_dim1 + 1], &c__1);
/* L50: */
	}
    }
    return 0;

/*     End of DLAROR */

} /* dlaror_slu */

