#include <string.h>
#include "../../SRC/slu_sdefs.h"

/* Subroutine */ int slaset_slu(char *uplo, int *m, int *n, float *alpha,
	float *beta, float *a, int *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower 
  
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper 
  
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) REAL   
            The constant to which the offdiagonal elements are to be set. 
  

    BETA    (input) REAL   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows: 
  

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    /* Local variables */
    static int i, j;


#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (strncmp(uplo, "U", 1)==0) {

/*        Set the strictly upper triangular or trapezoidal part of the
   
          array to ALPHA. */

	for (j = 2; j <= *n; ++j) {
/* Computing MIN */
	    for (i = 1; i <= SUPERLU_MIN(j-1,*m); ++i) {
		A(i,j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (strncmp(uplo, "L", 1)==0) {

/*        Set the strictly lower triangular or trapezoidal part of the
   
          array to ALPHA. */

	for (j = 1; j <= SUPERLU_MIN(*m,*n); ++j) {
	    for (i = j + 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first SUPERLU_MIN(M,N) diagonal elements to BETA. */

    for (i = 1; i <= SUPERLU_MIN(*m,*n); ++i) {
	A(i,i) = *beta;
/* L70: */
    }

    return 0;

/*     End of SLASET */

} /* slaset_slu */

