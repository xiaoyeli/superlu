/*  -- translated by f2c (version 19940927).
*/
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "../../SRC/slu_sdefs.h"

/* Table of constant values */

static int c__1 = 1;
static float c_b22 = 0.f;
static bool c_true = true;
static bool c_false = false;

/* Subroutine */ int slatms_slu(int *m, int *n, char *dist, int *
	iseed, char *sym, float *d, int *mode, float *cond, float *dmax__,
	int *kl, int *ku, char *pack, float *a, int *lda, float *
	work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    float r__1, r__2, r__3;
    bool L__1;

    /* Local variables */
    static int ilda, icol;
    static float temp;
    static int irow, isym;
    static float c;
    static int i, j, k;
    static float s, alpha, angle;
    static int ipack, ioffg;
    static int iinfo;
    extern /* Subroutine */ int sscal_(int *, float *, float *, int *);
    static int idist, mnmin, iskew;
    static float extra, dummy;
    extern int slatm1_slu(int *, float *, int *, int *,
	    int *, float *, int *, int *);
    static int ic, jc, nc, il, iendch, ir, jr, ipackg, mr;
    extern /* Subroutine */ int slagge_slu(int *, int *, int *,
	    int *, float *, float *, int *, int *, float *, int *
	    );
    static int minlda;
    extern int input_error(char *, int *);
    extern double slarnd_slu(int *, int *);
    static bool iltemp, givens;
    static int ioffst, irsign;
    extern /* Subroutine */ int slartg_slu(float *, float *, float *, float *, float *
	    ), slaset_slu(char *, int *, int *, float *, float *, float *,
	    int *), slagsy_slu(int *, int *, float *, float *,
	    int *, int *, float *, int *), slarot_slu(bool *,
	    bool *, bool *, int *, float *, float *, float *, int *
	    , float *, float *);
    static bool ilextr, topdwn;
    static int ir1, ir2, isympk, jch, llb, jkl, jku, uub;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       SLATMS generates random matrices with specified singular values   
       (or symmetric/hermitian with specified eigenvalues)   
       for testing LAPACK programs.   

       SLATMS operates by applying the following sequence of   
       operations:   

         Set the diagonal to D, where D may be input or   
            computed according to MODE, COND, DMAX, and SYM   
            as described below.   

         Generate a matrix with the appropriate band structure, by one   
            of two methods:   

         Method A:   
             Generate a dense M x N matrix by multiplying D on the left   
                 and the right by random unitary matrices, then:   

             Reduce the bandwidth according to KL and KU, using   
             Householder transformations.   

         Method B:   
             Convert the bandwidth-0 (i.e., diagonal) matrix to a   
                 bandwidth-1 matrix using Givens rotations, "chasing"   
                 out-of-band elements back, much as in QR; then   
                 convert the bandwidth-1 to a bandwidth-2 matrix, etc.   
                 Note that for reasonably small bandwidths (relative to   
                 M and N) this requires less storage, as a dense matrix   
                 is not generated.  Also, for symmetric matrices, only   
                 one triangle is generated.   

         Method A is chosen if the bandwidth is a large fraction of the   
             order of the matrix, and LDA is at least M (so a dense   
             matrix can be stored.)  Method B is chosen if the bandwidth 
  
             is small (< 1/2 N for symmetric, < .3 N+M for   
             non-symmetric), or LDA is less than M and not less than the 
  
             bandwidth.   

         Pack the matrix if desired. Options specified by PACK are:   
            no packing   
            zero out upper half (if symmetric)   
            zero out lower half (if symmetric)   
            store the upper half columnwise (if symmetric or upper   
                  triangular)   
            store the lower half columnwise (if symmetric or lower   
                  triangular)   
            store the lower triangle in banded format (if symmetric   
                  or lower triangular)   
            store the upper triangle in banded format (if symmetric   
                  or upper triangular)   
            store the entire matrix in banded format   
         If Method B is chosen, and band format is specified, then the   
            matrix will be generated in the band format, so no repacking 
  
            will be necessary.   

    Arguments   
    =========   

    M      - INTEGER   
             The number of rows of A. Not modified.   

    N      - INTEGER   
             The number of columns of A. Not modified.   

    DIST   - CHARACTER*1   
             On entry, DIST specifies the type of distribution to be used 
  
             to generate the random eigen-/singular values.   
             'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )   
             'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )   
             'N' => NORMAL( 0, 1 )   ( 'N' for normal )   
             Not modified.   

    ISEED  - INTEGER array, dimension ( 4 )   
             On entry ISEED specifies the seed of the random number   
             generator. They should lie between 0 and 4095 inclusive,   
             and ISEED(4) should be odd. The random number generator   
             uses a linear congruential sequence limited to small   
             integers, and so should produce machine independent   
             random numbers. The values of ISEED are changed on   
             exit, and can be used in the next call to SLATMS   
             to continue the same random number sequence.   
             Changed on exit.   

    SYM    - CHARACTER*1   
             If SYM='S' or 'H', the generated matrix is symmetric, with   
               eigenvalues specified by D, COND, MODE, and DMAX; they   
               may be positive, negative, or zero.   
             If SYM='P', the generated matrix is symmetric, with   
               eigenvalues (= singular values) specified by D, COND,   
               MODE, and DMAX; they will not be negative.   
             If SYM='N', the generated matrix is nonsymmetric, with   
               singular values specified by D, COND, MODE, and DMAX;   
               they will not be negative.   
             Not modified.   

    D      - REAL array, dimension ( MIN( M , N ) )   
             This array is used to specify the singular values or   
             eigenvalues of A (see SYM, above.)  If MODE=0, then D is   
             assumed to contain the singular/eigenvalues, otherwise   
             they will be computed according to MODE, COND, and DMAX,   
             and placed in D.   
             Modified if MODE is nonzero.   

    MODE   - INTEGER   
             On entry this describes how the singular/eigenvalues are to 
  
             be specified:   
             MODE = 0 means use D as input   
             MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND   
             MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND   
             MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))   
             MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)   
             MODE = 5 sets D to random numbers in the range   
                      ( 1/COND , 1 ) such that their logarithms   
                      are uniformly distributed.   
             MODE = 6 set D to random numbers from same distribution   
                      as the rest of the matrix.   
             MODE < 0 has the same meaning as ABS(MODE), except that   
                the order of the elements of D is reversed.   
             Thus if MODE is positive, D has entries ranging from   
                1 to 1/COND, if negative, from 1/COND to 1,   
             If SYM='S' or 'H', and MODE is neither 0, 6, nor -6, then   
                the elements of D will also be multiplied by a random   
                sign (i.e., +1 or -1.)   
             Not modified.   

    COND   - REAL   
             On entry, this is used as described under MODE above.   
             If used, it must be >= 1. Not modified.   

    DMAX   - REAL   
             If MODE is neither -6, 0 nor 6, the contents of D, as   
             computed according to MODE and COND, will be scaled by   
             DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or 
  
             singular value (which is to say the norm) will be abs(DMAX). 
  
             Note that DMAX need not be positive: if DMAX is negative   
             (or zero), D will be scaled by a negative number (or zero). 
  
             Not modified.   

    KL     - INTEGER   
             This specifies the lower bandwidth of the  matrix. For   
             example, KL=0 implies upper triangular, KL=1 implies upper   
             Hessenberg, and KL being at least M-1 means that the matrix 
  
             has full lower bandwidth.  KL must equal KU if the matrix   
             is symmetric.   
             Not modified.   

    KU     - INTEGER   
             This specifies the upper bandwidth of the  matrix. For   
             example, KU=0 implies lower triangular, KU=1 implies lower   
             Hessenberg, and KU being at least N-1 means that the matrix 
  
             has full upper bandwidth.  KL must equal KU if the matrix   
             is symmetric.   
             Not modified.   

    PACK   - CHARACTER*1   
             This specifies packing of matrix as follows:   
             'N' => no packing   
             'U' => zero out all subdiagonal entries (if symmetric)   
             'L' => zero out all superdiagonal entries (if symmetric)   
             'C' => store the upper triangle columnwise   
                    (only if the matrix is symmetric or upper triangular) 
  
             'R' => store the lower triangle columnwise   
                    (only if the matrix is symmetric or lower triangular) 
  
             'B' => store the lower triangle in band storage scheme   
                    (only if matrix symmetric or lower triangular)   
             'Q' => store the upper triangle in band storage scheme   
                    (only if matrix symmetric or upper triangular)   
             'Z' => store the entire matrix in band storage scheme   
                        (pivoting can be provided for by using this   
                        option to store A in the trailing rows of   
                        the allocated storage)   

             Using these options, the various LAPACK packed and banded   
             storage schemes can be obtained:   
             GB               - use 'Z'   
             PB, SB or TB     - use 'B' or 'Q'   
             PP, SP or TP     - use 'C' or 'R'   

             If two calls to SLATMS differ only in the PACK parameter,   
             they will generate mathematically equivalent matrices.   
             Not modified.   

    A      - REAL array, dimension ( LDA, N )   
             On exit A is the desired test matrix.  A is first generated 
  
             in full (unpacked) form, and then packed, if so specified   
             by PACK.  Thus, the first M elements of the first N   
             columns will always be modified.  If PACK specifies a   
             packed or banded storage scheme, all LDA elements of the   
             first N columns will be modified; the elements of the   
             array which do not correspond to elements of the generated   
             matrix are set to zero.   
             Modified.   

    LDA    - INTEGER   
             LDA specifies the first dimension of A as declared in the   
             calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then   
             LDA must be at least M.  If PACK='B' or 'Q', then LDA must   
             be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).   
             If PACK='Z', LDA must be large enough to hold the packed   
             array: MIN( KU, N-1) + MIN( KL, M-1) + 1.   
             Not modified.   

    WORK   - REAL array, dimension ( 3*MAX( N , M ) )   
             Workspace.   
             Modified.   

    INFO   - INTEGER   
             Error code.  On exit, INFO will be set to one of the   
             following values:   
               0 => normal return   
              -1 => M negative or unequal to N and SYM='S', 'H', or 'P'   
              -2 => N negative   
              -3 => DIST illegal string   
              -5 => SYM illegal string   
              -7 => MODE not in range -6 to 6   
              -8 => COND less than 1.0, and MODE neither -6, 0 nor 6   
             -10 => KL negative   
             -11 => KU negative, or SYM='S' or 'H' and KU not equal to KL 
  
             -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N'; 
  
                    or PACK='C' or 'Q' and SYM='N' and KL is not zero;   
                    or PACK='R' or 'B' and SYM='N' and KU is not zero;   
                    or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not 
  
                    N.   
             -14 => LDA is less than M, or PACK='Z' and LDA is less than 
  
                    MIN(KU,N-1) + MIN(KL,M-1) + 1.   
              1  => Error return from SLATM1   
              2  => Cannot scale to DMAX (max. sing. value is 0)   
              3  => Error return from SLAGGE or SLAGSY   

    ===================================================================== 
  


       1)      Decode and Test the input parameters.   
               Initialize flags & seed.   

       Parameter adjustments */
    --iseed;
    --d;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Decode DIST */

    if (strncmp(dist, "U", 1)==0) {
	idist = 1;
    } else if (strncmp(dist, "S", 1)==0) {
	idist = 2;
    } else if (strncmp(dist, "N", 1)==0) {
	idist = 3;
    } else {
	idist = -1;
    }

/*     Decode SYM */

    if (strncmp(sym, "N", 1)==0) {
	isym = 1;
	irsign = 0;
    } else if (strncmp(sym, "P", 1)==0) {
	isym = 2;
	irsign = 0;
    } else if (strncmp(sym, "S", 1)==0) {
	isym = 2;
	irsign = 1;
    } else if (strncmp(sym, "H", 1)==0) {
	isym = 2;
	irsign = 1;
    } else {
	isym = -1;
    }

/*     Decode PACK */

    isympk = 0;
    if (strncmp(pack, "N", 1)==0) {
	ipack = 0;
    } else if (strncmp(pack, "U", 1)==0) {
	ipack = 1;
	isympk = 1;
    } else if (strncmp(pack, "L", 1)==0) {
	ipack = 2;
	isympk = 1;
    } else if (strncmp(pack, "C", 1)==0) {
	ipack = 3;
	isympk = 2;
    } else if (strncmp(pack, "R", 1)==0) {
	ipack = 4;
	isympk = 3;
    } else if (strncmp(pack, "B", 1)==0) {
	ipack = 5;
	isympk = 3;
    } else if (strncmp(pack, "Q", 1)==0) {
	ipack = 6;
	isympk = 2;
    } else if (strncmp(pack, "Z", 1)==0) {
	ipack = 7;
    } else {
	ipack = -1;
    }

/*     Set certain internal parameters */

    mnmin = SUPERLU_MIN(*m,*n);
/* Computing MIN */
    i__1 = *kl, i__2 = *m - 1;
    llb = SUPERLU_MIN(i__1,i__2);
/* Computing MIN */
    i__1 = *ku, i__2 = *n - 1;
    uub = SUPERLU_MIN(i__1,i__2);
/* Computing MIN */
    i__1 = *m, i__2 = *n + llb;
    mr = SUPERLU_MIN(i__1,i__2);
/* Computing MIN */
    i__1 = *n, i__2 = *m + uub;
    nc = SUPERLU_MIN(i__1,i__2);

    if (ipack == 5 || ipack == 6) {
	minlda = uub + 1;
    } else if (ipack == 7) {
	minlda = llb + uub + 1;
    } else {
	minlda = *m;
    }

/*     Use Givens rotation method if bandwidth small enough,   
       or if LDA is too small to store the matrix unpacked. */

    givens = false;
    if (isym == 1) {
/* Computing MAX */
	i__1 = 1, i__2 = mr + nc;
	if ((float) (llb + uub) < (float) SUPERLU_MAX(i__1,i__2) * .3f) {
	    givens = true;
	}
    } else {
	if (llb << 1 < *m) {
	    givens = true;
	}
    }
    if (*lda < *m && *lda >= minlda) {
	givens = true;
    }

/*     Set INFO if an error */

    if (*m < 0) {
	*info = -1;
    } else if (*m != *n && isym != 1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (idist == -1) {
	*info = -3;
    } else if (isym == -1) {
	*info = -5;
    } else if (abs(*mode) > 6) {
	*info = -7;
    } else if (*mode != 0 && abs(*mode) != 6 && *cond < 1.f) {
	*info = -8;
    } else if (*kl < 0) {
	*info = -10;
    } else if (*ku < 0 || (isym != 1 && *kl != *ku)) {
	*info = -11;
    } else if ( ipack == -1 || (isympk == 1 && isym == 1) ||
		(isympk == 2 && isym == 1 && *kl > 0) ||
		(isympk == 3 && isym == 1 && *ku > 0) ||
		(isympk != 0 && *m != *n) ) {
	*info = -12;
    } else if (*lda < SUPERLU_MAX(1,minlda)) {
	*info = -14;
    }

    if (*info != 0) {
	i__1 = -(*info);
	input_error("SLATMS", &i__1);
	return 0;
    }

/*     Initialize random number generator */

    for (i = 1; i <= 4; ++i) {
	iseed[i] = (i__1 = iseed[i], abs(i__1)) % 4096;
/* L10: */
    }

    if (iseed[4] % 2 != 1) {
	++iseed[4];
    }

/*     2)      Set up D  if indicated.   

               Compute D according to COND and MODE */

    slatm1_slu(mode, cond, &irsign, &idist, &iseed[1], &d[1], &mnmin, &iinfo);
    if (iinfo != 0) {
	*info = 1;
	return 0;
    }

/*     Choose Top-Down if D is (apparently) increasing,   
       Bottom-Up if D is (apparently) decreasing. */

    if (fabs(d[1]) <= (r__1 = d[mnmin], fabs(r__1))) {
	topdwn = true;
    } else {
	topdwn = false;
    }

    if (*mode != 0 && abs(*mode) != 6) {

/*        Scale by DMAX */

	temp = fabs(d[1]);
	i__1 = mnmin;
	for (i = 2; i <= i__1; ++i) {
/* Computing MAX */
	    r__2 = temp, r__3 = (r__1 = d[i], fabs(r__1));
	    temp = fmax(r__2,r__3);
/* L20: */
	}

	if (temp > 0.f) {
	    alpha = *dmax__ / temp;
	} else {
	    *info = 2;
	    return 0;
	}

	sscal_(&mnmin, &alpha, &d[1], &c__1);

    }

/*     3)      Generate Banded Matrix using Givens rotations.   
               Also the special case of UUB=LLB=0   

                 Compute Addressing constants to cover all   
                 storage formats.  Whether GE, SY, GB, or SB,   
                 upper or lower triangle or both,   
                 the (i,j)-th element is in   
                 A( i - ISKEW*j + IOFFST, j ) */

    if (ipack > 4) {
	ilda = *lda - 1;
	iskew = 1;
	if (ipack > 5) {
	    ioffst = uub + 1;
	} else {
	    ioffst = 1;
	}
    } else {
	ilda = *lda;
	iskew = 0;
	ioffst = 0;
    }

/*     IPACKG is the format that the matrix is generated in. If this is   
       different from IPACK, then the matrix must be repacked at the   
       end.  It also signals how to compute the norm, for scaling. */

    ipackg = 0;
    slaset_slu("Full", lda, n, &c_b22, &c_b22, &a[a_offset], lda);

/*     Diagonal Matrix -- We are done, unless it   
       is to be stored SP/PP/TP (PACK='R' or 'C') */

    if (llb == 0 && uub == 0) {
	i__1 = ilda + 1;
	scopy_(&mnmin, &d[1], &c__1, &a[1 - iskew + ioffst + a_dim1], &i__1);
	if (ipack <= 2 || ipack >= 5) {
	    ipackg = ipack;
	}

    } else if (givens) {

/*        Check whether to use Givens rotations,   
          Householder transformations, or nothing. */

	if (isym == 1) {

/*           Non-symmetric -- A = U D V */

	    if (ipack > 4) {
		ipackg = ipack;
	    } else {
		ipackg = 0;
	    }

	    i__1 = ilda + 1;
	    scopy_(&mnmin, &d[1], &c__1, &a[1 - iskew + ioffst + a_dim1], &
		    i__1);

	    if (topdwn) {
		jkl = 0;
		i__1 = uub;
		for (jku = 1; jku <= i__1; ++jku) {

/*                 Transform from bandwidth JKL, JKU-1 to 
JKL, JKU   

                   Last row actually rotated is M   
                   Last column actually rotated is MIN( M+
JKU, N )   

   Computing MIN */
		    i__3 = *m + jku;
		    i__2 = SUPERLU_MIN(i__3,*n) + jkl - 1;
		    for (jr = 1; jr <= i__2; ++jr) {
			extra = 0.f;
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = sin(angle);
/* Computing MAX */
			i__3 = 1, i__4 = jr - jkl;
			icol = SUPERLU_MAX(i__3,i__4);
			if (jr < *m) {
/* Computing MIN */
			    i__3 = *n, i__4 = jr + jku;
			    il = SUPERLU_MIN(i__3,i__4) + 1 - icol;
			    L__1 = jr > jkl;
			    slarot_slu(&c_true, &L__1, &c_false, &il, &c, &s, &a[
				    jr - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &extra, &dummy);
			}

/*                    Chase "EXTRA" back up */

			ir = jr;
			ic = icol;
			i__3 = -jkl - jku;
			for (jch = jr - jkl; i__3 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__3) {
			    if (ir < *m) {
				slartg_slu(&a[ir + 1 - iskew * (ic + 1) + ioffst 
					+ (ic + 1) * a_dim1], &extra, &c, &s, 
					&dummy);
			    }
/* Computing MAX */
			    i__4 = 1, i__5 = jch - jku;
			    irow = SUPERLU_MAX(i__4,i__5);
			    il = ir + 2 - irow;
			    temp = 0.f;
			    iltemp = jch > jku;
			    r__1 = -(double)s;
			    slarot_slu(&c_false, &iltemp, &c_true, &il, &c, &
				    r__1, &a[irow - iskew * ic + ioffst + ic *
				     a_dim1], &ilda, &temp, &extra);
			    if (iltemp) {
				slartg_slu(&a[irow + 1 - iskew * (ic + 1) + 
					ioffst + (ic + 1) * a_dim1], &temp, &
					c, &s, &dummy);
/* Computing MAX */
				i__4 = 1, i__5 = jch - jku - jkl;
				icol = SUPERLU_MAX(i__4,i__5);
				il = ic + 2 - icol;
				extra = 0.f;
				L__1 = jch > jku + jkl;
				r__1 = -(double)s;
				slarot_slu(&c_true, &L__1, &c_true, &il, &c, &
					r__1, &a[irow - iskew * icol + ioffst 
					+ icol * a_dim1], &ilda, &extra, &
					temp);
				ic = icol;
				ir = irow;
			    }
/* L30: */
			}
/* L40: */
		    }
/* L50: */
		}

		jku = uub;
		i__1 = llb;
		for (jkl = 1; jkl <= i__1; ++jkl) {

/*                 Transform from bandwidth JKL-1, JKU to 
JKL, JKU   

   Computing MIN */
		    i__3 = *n + jkl;
		    i__2 = SUPERLU_MIN(i__3,*m) + jku - 1;
		    for (jc = 1; jc <= i__2; ++jc) {
			extra = 0.f;
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = sin(angle);
/* Computing MAX */
			i__3 = 1, i__4 = jc - jku;
			irow = SUPERLU_MAX(i__3,i__4);
			if (jc < *n) {
/* Computing MIN */
			    i__3 = *m, i__4 = jc + jkl;
			    il = SUPERLU_MIN(i__3,i__4) + 1 - irow;
			    L__1 = jc > jku;
			    slarot_slu(&c_false, &L__1, &c_false, &il, &c, &s, &
				    a[irow - iskew * jc + ioffst + jc * 
				    a_dim1], &ilda, &extra, &dummy);
			}

/*                    Chase "EXTRA" back up */

			ic = jc;
			ir = irow;
			i__3 = -jkl - jku;
			for (jch = jc - jku; i__3 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__3) {
			    if (ic < *n) {
				slartg_slu(&a[ir + 1 - iskew * (ic + 1) + ioffst 
					+ (ic + 1) * a_dim1], &extra, &c, &s, 
					&dummy);
			    }
/* Computing MAX */
			    i__4 = 1, i__5 = jch - jkl;
			    icol = SUPERLU_MAX(i__4,i__5);
			    il = ic + 2 - icol;
			    temp = 0.f;
			    iltemp = jch > jkl;
			    r__1 = -(double)s;
			    slarot_slu(&c_true, &iltemp, &c_true, &il, &c, &r__1,
				     &a[ir - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &temp, &extra);
			    if (iltemp) {
				slartg_slu(&a[ir + 1 - iskew * (icol + 1) + 
					ioffst + (icol + 1) * a_dim1], &temp, 
					&c, &s, &dummy);
/* Computing MAX */
				i__4 = 1, i__5 = jch - jkl - jku;
				irow = SUPERLU_MAX(i__4,i__5);
				il = ir + 2 - irow;
				extra = 0.f;
				L__1 = jch > jkl + jku;
				r__1 = -(double)s;
				slarot_slu(&c_false, &L__1, &c_true, &il, &c, &
					r__1, &a[irow - iskew * icol + ioffst 
					+ icol * a_dim1], &ilda, &extra, &
					temp);
				ic = icol;
				ir = irow;
			    }
/* L60: */
			}
/* L70: */
		    }
/* L80: */
		}

	    } else {

/*              Bottom-Up -- Start at the bottom right. */

		jkl = 0;
		i__1 = uub;
		for (jku = 1; jku <= i__1; ++jku) {

/*                 Transform from bandwidth JKL, JKU-1 to 
JKL, JKU   

                   First row actually rotated is M   
                   First column actually rotated is MIN( M
+JKU, N )   

   Computing MIN */
		    i__2 = *m, i__3 = *n + jkl;
		    iendch = SUPERLU_MIN(i__2,i__3) - 1;
/* Computing MIN */
		    i__2 = *m + jku;
		    i__3 = 1 - jkl;
		    for (jc = SUPERLU_MIN(i__2,*n) - 1; jc >= i__3; --jc) {
			extra = 0.f;
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = sin(angle);
/* Computing MAX */
			i__2 = 1, i__4 = jc - jku + 1;
			irow = SUPERLU_MAX(i__2,i__4);
			if (jc > 0) {
/* Computing MIN */
			    i__2 = *m, i__4 = jc + jkl + 1;
			    il = SUPERLU_MIN(i__2,i__4) + 1 - irow;
			    L__1 = jc + jkl < *m;
			    slarot_slu(&c_false, &c_false, &L__1, &il, &c, &s, &
				    a[irow - iskew * jc + ioffst + jc * 
				    a_dim1], &ilda, &dummy, &extra);
			}

/*                    Chase "EXTRA" back down */

			ic = jc;
			i__2 = iendch;
			i__4 = jkl + jku;
			for (jch = jc + jkl; i__4 < 0 ? jch >= i__2 : jch <= 
				i__2; jch += i__4) {
			    ilextr = ic > 0;
			    if (ilextr) {
				slartg_slu(&a[jch - iskew * ic + ioffst + ic * 
					a_dim1], &extra, &c, &s, &dummy);
			    }
			    ic = SUPERLU_MAX(1,ic);
/* Computing MIN */
			    i__5 = *n - 1, i__6 = jch + jku;
			    icol = SUPERLU_MIN(i__5,i__6);
			    iltemp = jch + jku < *n;
			    temp = 0.f;
			    i__5 = icol + 2 - ic;
			    slarot_slu(&c_true, &ilextr, &iltemp, &i__5, &c, &s, 
				    &a[jch - iskew * ic + ioffst + ic * 
				    a_dim1], &ilda, &extra, &temp);
			    if (iltemp) {
				slartg_slu(&a[jch - iskew * icol + ioffst + icol 
					* a_dim1], &temp, &c, &s, &dummy);
/* Computing MIN */
				i__5 = iendch, i__6 = jch + jkl + jku;
				il = SUPERLU_MIN(i__5,i__6) + 2 - jch;
				extra = 0.f;
				L__1 = jch + jkl + jku <= iendch;
				slarot_slu(&c_false, &c_true, &L__1, &il, &c, &s,
					 &a[jch - iskew * icol + ioffst + 
					icol * a_dim1], &ilda, &temp, &extra);
				ic = icol;
			    }
/* L90: */
			}
/* L100: */
		    }
/* L110: */
		}

		jku = uub;
		i__1 = llb;
		for (jkl = 1; jkl <= i__1; ++jkl) {

/*                 Transform from bandwidth JKL-1, JKU to 
JKL, JKU   

                   First row actually rotated is MIN( N+JK
L, M )   
                   First column actually rotated is N   

   Computing MIN */
		    i__3 = *n, i__4 = *m + jku;
		    iendch = SUPERLU_MIN(i__3,i__4) - 1;
/* Computing MIN */
		    i__3 = *n + jkl;
		    i__4 = 1 - jku;
		    for (jr = SUPERLU_MIN(i__3,*m) - 1; jr >= i__4; --jr) {
			extra = 0.f;
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = sin(angle);
/* Computing MAX */
			i__3 = 1, i__2 = jr - jkl + 1;
			icol = SUPERLU_MAX(i__3,i__2);
			if (jr > 0) {
/* Computing MIN */
			    i__3 = *n, i__2 = jr + jku + 1;
			    il = SUPERLU_MIN(i__3,i__2) + 1 - icol;
			    L__1 = jr + jku < *n;
			    slarot_slu(&c_true, &c_false, &L__1, &il, &c, &s, &a[
				    jr - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &dummy, &extra);
			}

/*                    Chase "EXTRA" back down */

			ir = jr;
			i__3 = iendch;
			i__2 = jkl + jku;
			for (jch = jr + jku; i__2 < 0 ? jch >= i__3 : jch <= 
				i__3; jch += i__2) {
			    ilextr = ir > 0;
			    if (ilextr) {
				slartg_slu(&a[ir - iskew * jch + ioffst + jch * 
					a_dim1], &extra, &c, &s, &dummy);
			    }
			    ir = SUPERLU_MAX(1,ir);
/* Computing MIN */
			    i__5 = *m - 1, i__6 = jch + jkl;
			    irow = SUPERLU_MIN(i__5,i__6);
			    iltemp = jch + jkl < *m;
			    temp = 0.f;
			    i__5 = irow + 2 - ir;
			    slarot_slu(&c_false, &ilextr, &iltemp, &i__5, &c, &s,
				     &a[ir - iskew * jch + ioffst + jch * 
				    a_dim1], &ilda, &extra, &temp);
			    if (iltemp) {
				slartg_slu(&a[irow - iskew * jch + ioffst + jch *
					 a_dim1], &temp, &c, &s, &dummy);
/* Computing MIN */
				i__5 = iendch, i__6 = jch + jkl + jku;
				il = SUPERLU_MIN(i__5,i__6) + 2 - jch;
				extra = 0.f;
				L__1 = jch + jkl + jku <= iendch;
				slarot_slu(&c_true, &c_true, &L__1, &il, &c, &s, 
					&a[irow - iskew * jch + ioffst + jch *
					 a_dim1], &ilda, &temp, &extra);
				ir = irow;
			    }
/* L120: */
			}
/* L130: */
		    }
/* L140: */
		}
	    }

	} else {

/*           Symmetric -- A = U D U' */

	    ipackg = ipack;
	    ioffg = ioffst;

	    if (topdwn) {

/*              Top-Down -- Generate Upper triangle only */

		if (ipack >= 5) {
		    ipackg = 6;
		    ioffg = uub + 1;
		} else {
		    ipackg = 1;
		}
		i__1 = ilda + 1;
		scopy_(&mnmin, &d[1], &c__1, &a[1 - iskew + ioffg + a_dim1], &
			i__1);

		i__1 = uub;
		for (k = 1; k <= i__1; ++k) {
		    i__4 = *n - 1;
		    for (jc = 1; jc <= i__4; ++jc) {
/* Computing MAX */
			i__2 = 1, i__3 = jc - k;
			irow = SUPERLU_MAX(i__2,i__3);
/* Computing MIN */
			i__2 = jc + 1, i__3 = k + 2;
			il = SUPERLU_MIN(i__2,i__3);
			extra = 0.f;
			temp = a[jc - iskew * (jc + 1) + ioffg + (jc + 1) * 
				a_dim1];
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = sin(angle);
			L__1 = jc > k;
			slarot_slu(&c_false, &L__1, &c_true, &il, &c, &s, &a[
				irow - iskew * jc + ioffg + jc * a_dim1], &
				ilda, &extra, &temp);
/* Computing MIN */
			i__3 = k, i__5 = *n - jc;
			i__2 = SUPERLU_MIN(i__3,i__5) + 1;
			slarot_slu(&c_true, &c_true, &c_false, &i__2, &c, &s, &a[
				(1 - iskew) * jc + ioffg + jc * a_dim1], &
				ilda, &temp, &dummy);

/*                    Chase EXTRA back up the matrix 
*/

			icol = jc;
			i__2 = -k;
			for (jch = jc - k; i__2 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__2) {
			    slartg_slu(&a[jch + 1 - iskew * (icol + 1) + ioffg + 
				    (icol + 1) * a_dim1], &extra, &c, &s, &
				    dummy);
			    temp = a[jch - iskew * (jch + 1) + ioffg + (jch + 
				    1) * a_dim1];
			    i__3 = k + 2;
			    r__1 = -(double)s;
			    slarot_slu(&c_true, &c_true, &c_true, &i__3, &c, &
				    r__1, &a[(1 - iskew) * jch + ioffg + jch *
				     a_dim1], &ilda, &temp, &extra);
/* Computing MAX */
			    i__3 = 1, i__5 = jch - k;
			    irow = SUPERLU_MAX(i__3,i__5);
/* Computing MIN */
			    i__3 = jch + 1, i__5 = k + 2;
			    il = SUPERLU_MIN(i__3,i__5);
			    extra = 0.f;
			    L__1 = jch > k;
			    r__1 = -(double)s;
			    slarot_slu(&c_false, &L__1, &c_true, &il, &c, &r__1, 
				    &a[irow - iskew * jch + ioffg + jch * 
				    a_dim1], &ilda, &extra, &temp);
			    icol = jch;
/* L150: */
			}
/* L160: */
		    }
/* L170: */
		}

/*              If we need lower triangle, copy from upper. No
te that   
                the order of copying is chosen to work for 'q'
 -> 'b' */

		if (ipack != ipackg && ipack != 3) {
		    i__1 = *n;
		    for (jc = 1; jc <= i__1; ++jc) {
			irow = ioffst - iskew * jc;
/* Computing MIN */
			i__2 = *n, i__3 = jc + uub;
			i__4 = SUPERLU_MIN(i__2,i__3);
			for (jr = jc; jr <= i__4; ++jr) {
			    a[jr + irow + jc * a_dim1] = a[jc - iskew * jr + 
				    ioffg + jr * a_dim1];
/* L180: */
			}
/* L190: */
		    }
		    if (ipack == 5) {
			i__1 = *n;
			for (jc = *n - uub + 1; jc <= i__1; ++jc) {
			    i__4 = uub + 1;
			    for (jr = *n + 2 - jc; jr <= i__4; ++jr) {
				a[jr + jc * a_dim1] = 0.f;
/* L200: */
			    }
/* L210: */
			}
		    }
		    if (ipackg == 6) {
			ipackg = ipack;
		    } else {
			ipackg = 0;
		    }
		}
	    } else {

/*              Bottom-Up -- Generate Lower triangle only */

		if (ipack >= 5) {
		    ipackg = 5;
		    if (ipack == 6) {
			ioffg = 1;
		    }
		} else {
		    ipackg = 2;
		}
		i__1 = ilda + 1;
		scopy_(&mnmin, &d[1], &c__1, &a[1 - iskew + ioffg + a_dim1], &
			i__1);

		i__1 = uub;
		for (k = 1; k <= i__1; ++k) {
		    for (jc = *n - 1; jc >= 1; --jc) {
/* Computing MIN */
			i__4 = *n + 1 - jc, i__2 = k + 2;
			il = SUPERLU_MIN(i__4,i__2);
			extra = 0.f;
			temp = a[(1 - iskew) * jc + 1 + ioffg + jc * a_dim1];
			angle = slarnd_slu(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			c = cos(angle);
			s = -(double)sin(angle);
			L__1 = *n - jc > k;
			slarot_slu(&c_false, &c_true, &L__1, &il, &c, &s, &a[(1 
				- iskew) * jc + ioffg + jc * a_dim1], &ilda, &
				temp, &extra);
/* Computing MAX */
			i__4 = 1, i__2 = jc - k + 1;
			icol = SUPERLU_MAX(i__4,i__2);
			i__4 = jc + 2 - icol;
			slarot_slu(&c_true, &c_false, &c_true, &i__4, &c, &s, &a[
				jc - iskew * icol + ioffg + icol * a_dim1], &
				ilda, &dummy, &temp);

/*                    Chase EXTRA back down the matrix
 */

			icol = jc;
			i__4 = *n - 1;
			i__2 = k;
			for (jch = jc + k; i__2 < 0 ? jch >= i__4 : jch <= 
				i__4; jch += i__2) {
			    slartg_slu(&a[jch - iskew * icol + ioffg + icol * 
				    a_dim1], &extra, &c, &s, &dummy);
			    temp = a[(1 - iskew) * jch + 1 + ioffg + jch * 
				    a_dim1];
			    i__3 = k + 2;
			    slarot_slu(&c_true, &c_true, &c_true, &i__3, &c, &s, 
				    &a[jch - iskew * icol + ioffg + icol * 
				    a_dim1], &ilda, &extra, &temp);
/* Computing MIN */
			    i__3 = *n + 1 - jch, i__5 = k + 2;
			    il = SUPERLU_MIN(i__3,i__5);
			    extra = 0.f;
			    L__1 = *n - jch > k;
			    slarot_slu(&c_false, &c_true, &L__1, &il, &c, &s, &a[
				    (1 - iskew) * jch + ioffg + jch * a_dim1],
				     &ilda, &temp, &extra);
			    icol = jch;
/* L220: */
			}
/* L230: */
		    }
/* L240: */
		}

/*              If we need upper triangle, copy from lower. No
te that   
                the order of copying is chosen to work for 'b'
 -> 'q' */

		if (ipack != ipackg && ipack != 4) {
		    for (jc = *n; jc >= 1; --jc) {
			irow = ioffst - iskew * jc;
/* Computing MAX */
			i__2 = 1, i__4 = jc - uub;
			i__1 = SUPERLU_MAX(i__2,i__4);
			for (jr = jc; jr >= i__1; --jr) {
			    a[jr + irow + jc * a_dim1] = a[jc - iskew * jr + 
				    ioffg + jr * a_dim1];
/* L250: */
			}
/* L260: */
		    }
		    if (ipack == 6) {
			i__1 = uub;
			for (jc = 1; jc <= i__1; ++jc) {
			    i__2 = uub + 1 - jc;
			    for (jr = 1; jr <= i__2; ++jr) {
				a[jr + jc * a_dim1] = 0.f;
/* L270: */
			    }
/* L280: */
			}
		    }
		    if (ipackg == 5) {
			ipackg = ipack;
		    } else {
			ipackg = 0;
		    }
		}
	    }
	}

    } else {

/*        4)      Generate Banded Matrix by first   
                  Rotating by random Unitary matrices,   
                  then reducing the bandwidth using Householder   
                  transformations.   

                  Note: we should get here only if LDA .ge. N */

	if (isym == 1) {

/*           Non-symmetric -- A = U D V */

	    slagge_slu(&mr, &nc, &llb, &uub, &d[1], &a[a_offset], lda, &iseed[1],
		     &work[1], &iinfo);
	} else {

/*           Symmetric -- A = U D U' */

	    slagsy_slu(m, &llb, &d[1], &a[a_offset], lda, &iseed[1], &work[1], &
		    iinfo);

	}
	if (iinfo != 0) {
	    *info = 3;
	    return 0;
	}
    }

/*     5)      Pack the matrix */

    if (ipack != ipackg) {
	if (ipack == 1) {

/*           'U' -- Upper triangular, not packed */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = j + 1; i <= i__2; ++i) {
		    a[i + j * a_dim1] = 0.f;
/* L290: */
		}
/* L300: */
	    }

	} else if (ipack == 2) {

/*           'L' -- Lower triangular, not packed */

	    i__1 = *m;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		for (i = 1; i <= i__2; ++i) {
		    a[i + j * a_dim1] = 0.f;
/* L310: */
		}
/* L320: */
	    }

	} else if (ipack == 3) {

/*           'C' -- Upper triangle packed Columnwise. */

	    icol = 1;
	    irow = 0;
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i = 1; i <= i__2; ++i) {
		    ++irow;
		    if (irow > *lda) {
			irow = 1;
			++icol;
		    }
		    a[irow + icol * a_dim1] = a[i + j * a_dim1];
/* L330: */
		}
/* L340: */
	    }

	} else if (ipack == 4) {

/*           'R' -- Lower triangle packed Columnwise. */

	    icol = 1;
	    irow = 0;
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = j; i <= i__2; ++i) {
		    ++irow;
		    if (irow > *lda) {
			irow = 1;
			++icol;
		    }
		    a[irow + icol * a_dim1] = a[i + j * a_dim1];
/* L350: */
		}
/* L360: */
	    }

	} else if (ipack >= 5) {

/*           'B' -- The lower triangle is packed as a band matrix.
   
             'Q' -- The upper triangle is packed as a band matrix.
   
             'Z' -- The whole matrix is packed as a band matrix. 
*/

	    if (ipack == 5) {
		uub = 0;
	    }
	    if (ipack == 6) {
		llb = 0;
	    }

	    i__1 = uub;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__2 = j + llb;
		for (i = SUPERLU_MIN(i__2,*m); i >= 1; --i) {
		    a[i - j + uub + 1 + j * a_dim1] = a[i + j * a_dim1];
/* L370: */
		}
/* L380: */
	    }

	    i__1 = *n;
	    for (j = uub + 2; j <= i__1; ++j) {
/* Computing MIN */
		i__4 = j + llb;
		i__2 = SUPERLU_MIN(i__4,*m);
		for (i = j - uub; i <= i__2; ++i) {
		    a[i - j + uub + 1 + j * a_dim1] = a[i + j * a_dim1];
/* L390: */
		}
/* L400: */
	    }
	}

/*        If packed, zero out extraneous elements.   

          Symmetric/Triangular Packed --   
          zero out everything after A(IROW,ICOL) */

	if (ipack == 3 || ipack == 4) {
	    i__1 = *m;
	    for (jc = icol; jc <= i__1; ++jc) {
		i__2 = *lda;
		for (jr = irow + 1; jr <= i__2; ++jr) {
		    a[jr + jc * a_dim1] = 0.f;
/* L410: */
		}
		irow = 0;
/* L420: */
	    }

	} else if (ipack >= 5) {

/*           Packed Band --   
                1st row is now in A( UUB+2-j, j), zero above it   
                m-th row is now in A( M+UUB-j,j), zero below it   
                last non-zero diagonal is now in A( UUB+LLB+1,j ),
   
                   zero below it, too. */

	    ir1 = uub + llb + 2;
	    ir2 = uub + *m + 2;
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		i__2 = uub + 1 - jc;
		for (jr = 1; jr <= i__2; ++jr) {
		    a[jr + jc * a_dim1] = 0.f;
/* L430: */
		}
/* Computing MAX   
   Computing MIN */
		i__3 = ir1, i__5 = ir2 - jc;
		i__2 = 1, i__4 = SUPERLU_MIN(i__3,i__5);
		i__6 = *lda;
		for (jr = SUPERLU_MAX(i__2,i__4); jr <= i__6; ++jr) {
		    a[jr + jc * a_dim1] = 0.f;
/* L440: */
		}
/* L450: */
	    }
	}
    }

    return 0;

/*     End of SLATMS */

} /* slatms_slu */

