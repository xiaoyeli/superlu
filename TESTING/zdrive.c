/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 7.0.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 * August 2024
 *
 */

/*! \file
 * ZDRIVE is the main test program for the DOUBLE COMPLEX linear 
 * equation driver routines ZGSSV and ZGSSVX.
 * 
 * The program is invoked by a shell script file -- ztest.csh.
 * The output from the tests are written into a file -- ztest.out.
 *
 * \ingroup TestingZ
 */
#include <getopt.h>
#include <string.h>
#include "slu_zdefs.h"
#include "MATGEN/matgen.h"

#define NTESTS    5      /* Number of test types */
#define NTYPES    11     /* Number of matrix types */
#define NTRAN     2    
#define THRESH    20.0
#define FMT1      "%10s:n=%d, test(%d)=%12.5g\n"
#define	FMT2      "%10s:fact=%4d, trans=%4d, equed=%c, n=%d, imat=%d, test(%d)=%12.5g\n"
#define FMT3      "%10s:info=%d, izero=%d, n=%d, nrhs=%d, imat=%d, nfail=%d\n"

static void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *n, int *w, int *relax, int *nrhs, int *maxsuper,
		   int *rowblk, int *colblk, int_t *lwork, double *u, FILE **fp);

/*!
 * Entry point of test program.
 */
int main(int argc, char *argv[])
{
    doublecomplex         *a, *a_save;
    int_t          *asub, *asub_save;
    int_t          *xa, *xa_save;
    SuperMatrix  A, B, X, L, U;
    SuperMatrix  ASAV, AC;
    GlobalLU_t   Glu; /* Not needed on return. */
    mem_usage_t    mem_usage;
    int            *perm_r; /* row permutation from partial pivoting */
    int            *perm_c, *pc_save; /* column permutation */
    int            *etree;
    doublecomplex  zero = {0.0, 0.0};
    double         *R, *C;
    double         *ferr, *berr;
    double         *rwork;
    doublecomplex	   *wwork;
    void           *work = NULL;
    int            nrhs, panel_size, relax;
    int            m, n, info1;
    int_t          nnz, lwork, info;
    doublecomplex         *xact;
    doublecomplex         *rhsb, *solx, *bsav;
    int            ldb, ldx;
    double         rpg, rcond;
    int            i, j, k1;
    double         rowcnd, colcnd, amax;
    int            maxsuper, rowblk, colblk;
    int            prefact, equil, iequed;
    int            nt, nrun, nfail, nerrs, imat, fimat, nimat;
    int            nfact, ifact, itran;
    int            kl, ku, mode, lda, ioff;
    int            zerot; /* indicate whether the matrix is singular */
    int            izero; /* indicate the first column that is entirely zero */
    double         u;
    double         anorm, cndnum;
    doublecomplex         *Afull;
    double         result[NTESTS];
    superlu_options_t options;
    fact_t         fact;
    trans_t        trans;
    SuperLUStat_t  stat;
    static char    matrix_type[8];
    static char    equed[1], path[4], sym[1], dist[1];
    FILE           *fp;

    /* Fixed set of parameters */
    int            iseed[]  = {1988, 1989, 1990, 1991};
    static char    equeds[]  = {'N', 'R', 'C', 'B'};
    static fact_t  facts[] = {FACTORED, DOFACT, SamePattern,
			      SamePattern_SameRowPerm};
    static trans_t transs[]  = {NOTRANS, TRANS, CONJ};

    /* Some function prototypes */ 
    extern int zgst01(int, int, SuperMatrix *, SuperMatrix *, 
		      SuperMatrix *, int *, int *, double *);
    extern int zgst02(trans_t, int, int, int, SuperMatrix *, doublecomplex *,
                      int, doublecomplex *, int, double *resid);
    extern int zgst04(int, int, doublecomplex *, int, 
                      doublecomplex *, int, double rcond, double *resid);
    extern int zgst07(trans_t, int, int, SuperMatrix *, doublecomplex *, int,
                         doublecomplex *, int, doublecomplex *, int, 
                         double *, double *, double *);
    extern int zlatb4_slu(char *, int *, int *, int *, char *, int *, int *, 
	               double *, int *, double *, char *);
    extern int zlatms_slu(int *, int *, char *, int *, char *, double *d,
                       int *, double *, double *, int *, int *,
                       char *, doublecomplex *, int *, doublecomplex *, int *);
    extern int sp_zconvert(int, int, doublecomplex *, int, int, int,
	                   doublecomplex *a, int_t *, int_t *, int_t *);


    /* Executable statements */

    strcpy(path, "ZGE");
    nrun  = 0;
    nfail = 0;
    nerrs = 0;

    /* Defaults */
    lwork      = 0;
    n          = 1;
    nrhs       = 1;
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    u          = 1.0;
    strcpy(matrix_type, "LA");
    parse_command_line(argc, argv, matrix_type, &n,
		       &panel_size, &relax, &nrhs, &maxsuper,
		       &rowblk, &colblk, &lwork, &u, &fp);
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    fprintf(stderr, "expert: cannot allocate %lld bytes\n", (long long) lwork);
	    exit (-1);
	}
    }

    /* Set the default input options. */
    set_default_options(&options);
    options.DiagPivotThresh = u;
    options.PrintStat = NO;
    options.PivotGrowth = YES;
    options.ConditionNumber = YES;
    options.IterRefine = SLU_DOUBLE;
    
    if ( strcmp(matrix_type, "LA") == 0 ) {
	/* Test LAPACK matrix suite. */
	m = n;
	lda = SUPERLU_MAX(n, 1);
	nnz = n * n;        /* upper bound */
	fimat = 1;
	nimat = NTYPES;
	Afull = doublecomplexCalloc(lda * n);
	zallocateA(n, nnz, &a, &asub, &xa);
    } else {
	/* Read a sparse matrix */
	fimat = nimat = 0;
	zreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
    }

    zallocateA(n, nnz, &a_save, &asub_save, &xa_save);
    rhsb = doublecomplexMalloc(m * nrhs);
    bsav = doublecomplexMalloc(m * nrhs);
    solx = doublecomplexMalloc(n * nrhs);
    xact = doublecomplexMalloc(n * nrhs);
    wwork = doublecomplexCalloc( SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs) );

    ldb  = m;
    ldx  = n;
    zCreate_Dense_Matrix(&B, m, nrhs, rhsb, ldb, SLU_DN, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix(&X, n, nrhs, solx, ldx, SLU_DN, SLU_Z, SLU_GE);
    etree   = int32Malloc(n);
    perm_r  = int32Malloc(n);
    perm_c  = int32Malloc(n);
    pc_save = int32Malloc(n);
    R       = (double *) SUPERLU_MALLOC(m*sizeof(double));
    C       = (double *) SUPERLU_MALLOC(n*sizeof(double));
    ferr    = (double *) SUPERLU_MALLOC(nrhs*sizeof(double));
    berr    = (double *) SUPERLU_MALLOC(nrhs*sizeof(double));
    j = SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs);    
    rwork   = (double *) SUPERLU_MALLOC(j*sizeof(double));
    for (i = 0; i < j; ++i) rwork[i] = 0.;
    if ( !R ) ABORT("SUPERLU_MALLOC fails for R");
    if ( !C ) ABORT("SUPERLU_MALLOC fails for C");
    if ( !ferr ) ABORT("SUPERLU_MALLOC fails for ferr");
    if ( !berr ) ABORT("SUPERLU_MALLOC fails for berr");
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");

    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i] = i;
    options.ColPerm = MY_PERMC;

    for (imat = fimat; imat <= nimat; ++imat) { /* All matrix types */
	
	if ( imat ) {

	    /* Skip types 5, 6, or 7 if the matrix size is too small. */
	    zerot = (imat >= 5 && imat <= 7);
	    if ( zerot && n < imat-4 )
		continue;
	    
	    /* Set up parameters with ZLATB4 and generate a test matrix
	       with ZLATMS.  */
	    zlatb4_slu(path, &imat, &n, &n, sym, &kl, &ku, &anorm, &mode,
		    &cndnum, dist);

	    zlatms_slu(&n, &n, dist, iseed, sym, &rwork[0], &mode, &cndnum,
		    &anorm, &kl, &ku, "No packing", Afull, &lda,
		    &wwork[0], &info1);

	    if ( info1 ) {
		printf(FMT3, "ZLATMS", info1, izero, n, nrhs, imat, nfail);
		continue;
	    }

	    /* For types 5-7, zero one or more columns of the matrix
	       to test that INFO is returned correctly.   */
	    if ( zerot ) {
		if ( imat == 5 ) izero = 1;
		else if ( imat == 6 ) izero = n;
		else izero = n / 2 + 1;
		ioff = (izero - 1) * lda;
		if ( imat < 7 ) {
		    for (i = 0; i < n; ++i) Afull[ioff + i] = zero;
		} else {
		    for (j = 0; j < n - izero + 1; ++j)
			for (i = 0; i < n; ++i)
			    Afull[ioff + i + j*lda] = zero;
		}
	    } else {
		izero = n+1; /* none of the column is zero */
	    }

	    /* Convert to sparse representation. */
	    sp_zconvert(n, n, Afull, lda, kl, ku, a, asub, xa, &nnz);

	} else {
	    izero = n+1; /* none of the column is zero */
	    zerot = 0;
	}
	
	zCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_Z, SLU_GE);

	/* Save a copy of matrix A in ASAV */
	zCreate_CompCol_Matrix(&ASAV, m, n, nnz, a_save, asub_save, xa_save,
			      SLU_NC, SLU_Z, SLU_GE);
	zCopy_CompCol_Matrix(&A, &ASAV);
	
	/* Form exact solution. */
	zGenXtrue(n, nrhs, xact, ldx);
	
	StatInit(&stat);

	for (iequed = 0; iequed < 4; ++iequed) {
	    *equed = equeds[iequed];
	    if (iequed == 0) nfact = 4;
	    else nfact = 1; /* Only test factored, pre-equilibrated matrix */

	    for (ifact = 0; ifact < nfact; ++ifact) {
		fact = facts[ifact];
		options.Fact = fact;

		for (equil = 0; equil < 2; ++equil) {
		    options.Equil = equil;
		    prefact   = ( options.Fact == FACTORED ||
				  options.Fact == SamePattern_SameRowPerm );
                                /* Need a first factor */

		    /* Restore the matrix A. */
		    zCopy_CompCol_Matrix(&ASAV, &A);
			
		    if ( zerot ) {
                        if ( prefact ) continue;
		    } else if ( options.Fact == FACTORED ) {
                        if ( equil || iequed ) {
			    /* Compute row and column scale factors to
			       equilibrate matrix A.    */
			    zgsequ(&A, R, C, &rowcnd, &colcnd, &amax, &info1);

			    /* Force equilibration. */
			    if ( !info && n > 0 ) {
				if ( strncmp(equed, "R", 1)==0 ) {
				    rowcnd = 0.;
				    colcnd = 1.;
				} else if ( strncmp(equed, "C", 1)==0 ) {
				    rowcnd = 1.;
				    colcnd = 0.;
				} else if ( strncmp(equed, "B", 1)==0 ) {
				    rowcnd = 0.;
				    colcnd = 0.;
				}
			    }
			
			    /* Equilibrate the matrix. */
			    zlaqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
		    }
		    
		    if ( prefact ) { /* Need a factor for the first time */
			
		        /* Save Fact option. */
		        fact = options.Fact;
			options.Fact = DOFACT;

			/* Preorder the matrix, obtain the column etree. */
			sp_preorder(&options, &A, perm_c, etree, &AC);

			/* Factor the matrix AC. */
			zgstrf(&options, &AC, relax, panel_size,
                               etree, work, lwork, perm_c, perm_r, &L, &U,
                               &Glu, &stat, &info);

			if ( info ) { 
                            printf("** First factor: info %lld, equed %c\n",
				   (long long) info, *equed);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %lld bytes\n",
                                        (long long) info - n);
                                exit(0);
                            }
                        }
	
                        Destroy_CompCol_Permuted(&AC);
			
		        /* Restore Fact option. */
			options.Fact = fact;
		    } /* if .. first time factor */
		    
		    for (itran = 0; itran < NTRAN; ++itran) {
			trans = transs[itran];
                        options.Trans = trans;

			/* Restore the matrix A. */
			zCopy_CompCol_Matrix(&ASAV, &A);
			
 			/* Set the right hand side. */
			zFillRHS(trans, nrhs, xact, ldx, &A, &B);
			zCopy_Dense_Matrix(m, nrhs, rhsb, ldb, bsav, ldb);

			/*----------------
			 * Test zgssv
			 *----------------*/
			if ( options.Fact == DOFACT && itran == 0) {
                            /* Not yet factored, and untransposed */
	
			    zCopy_Dense_Matrix(m, nrhs, rhsb, ldb, solx, ldx);
			    zgssv(&options, &A, perm_c, perm_r, &L, &U, &X,
                                  &stat, &info);
			    
			    if ( info && info != izero ) {
                                printf(FMT3, "zgssv",
				       (int) info, izero, n, nrhs, imat, nfail);
			    } else {
                                /* Reconstruct matrix from factors and compute residual.
				 * Only compute the leading 'izero' nonzero columns.
				 */
                                zgst01(m, izero-1, &A, &L, &U, perm_c, perm_r,
                                         &result[0]);
				nt = 1;
				if ( izero == (n+1) ) {
				    /* Compute residual of the computed
				       solution. */
				    zCopy_Dense_Matrix(m, nrhs, rhsb, ldb,
						       wwork, ldb);
				    zgst02(trans, m, n, nrhs, &A, solx,
                                              ldx, wwork,ldb, &result[1]);
				    nt = 2;
				}
				
				/* Print information about the tests that
				   did not pass the threshold.      */
				for (i = 0; i < nt; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT1, "zgssv", n, i,
					       result[i]);
					++nfail;
				    }
				}
				nrun += nt;
			    } /* else .. info == 0 */

			    /* Restore perm_c. */
			    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i];

		            if (lwork == 0) {
			        Destroy_SuperNode_Matrix(&L);
			        Destroy_CompCol_Matrix(&U);
			    }
			} /* if .. end of testing zgssv */
    
			/*----------------
			 * Test zgssvx
			 *----------------*/
    
			/* Equilibrate the matrix if fact = FACTORED and
			   equed = 'R', 'C', or 'B'.   */
			if ( options.Fact == FACTORED &&
			     (equil || iequed) && n > 0 ) {
			    zlaqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
			
			/* Solve the system and compute the condition number
			   and error bounds using zgssvx.      */
			zgssvx(&options, &A, perm_c, perm_r, etree,
                               equed, R, C, &L, &U, work, lwork, &B, &X, &rpg,
                               &rcond, ferr, berr, &Glu,
			       &mem_usage, &stat, &info);

			if ( info && info != izero ) {
			    printf(FMT3, "zgssvx",
				   (int) info, izero, n, nrhs, imat, nfail);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %.0f bytes\n",
                                        mem_usage.total_needed);
                                exit(0);
                            }
			} else {
			    if ( !prefact ) {
			    	/* Reconstruct matrix from factors and compute residual.
				 * Only compute the leading 'izero' nonzero columns.
				 */
                                zgst01(m, izero-1, &A, &L, &U, perm_c, perm_r,
                                         &result[0]);
				k1 = 0;
			    } else {
			   	k1 = 1;
			    }

			    if ( !info ) {
				/* Compute residual of the computed solution.*/
				zCopy_Dense_Matrix(m, nrhs, bsav, ldb,
						  wwork, ldb);
				zgst02(trans, m, n, nrhs, &ASAV, solx, ldx,
					  wwork, ldb, &result[1]);

				/* Check solution from generated exact
				   solution. */
				zgst04(n, nrhs, solx, ldx, xact, ldx, rcond,
					  &result[2]);

				/* Check the error bounds from iterative
				   refinement. */
				zgst07(trans, n, nrhs, &ASAV, bsav, ldb,
					  solx, ldx, xact, ldx, ferr, berr,
					  &result[3]);

				/* Print information about the tests that did
				   not pass the threshold.    */
				for (i = k1; i < NTESTS; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT2, "zgssvx",
					       options.Fact, trans, *equed,
					       n, imat, i, result[i]);
					++nfail;
				    }
				}
				nrun += NTESTS;
			    } /* if .. info == 0 */
			} /* else .. end of testing zgssvx */

		    } /* for itran ... */

		    if ( lwork == 0 ) {
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
		    }

		} /* for equil ... */
	    } /* for ifact ... */
	} /* for iequed ... */
#if 0    
    if ( !info ) {
	PrintPerf(&L, &U, &mem_usage, rpg, rcond, ferr, berr, equed);
    }
#endif
        Destroy_SuperMatrix_Store(&A);
        Destroy_SuperMatrix_Store(&ASAV);
        StatFree(&stat);

    } /* for imat ... */

    /* Print a summary of the results. */
    PrintSumm("ZGE", nfail, nrun, nerrs);

    if ( strcmp(matrix_type, "LA") == 0 ) SUPERLU_FREE (Afull);
    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (bsav);
    SUPERLU_FREE (solx);    
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (pc_save);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    SUPERLU_FREE (rwork);
    SUPERLU_FREE (wwork);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
#if 0
    Destroy_CompCol_Matrix(&A);
    Destroy_CompCol_Matrix(&ASAV);
#else
    SUPERLU_FREE(a); SUPERLU_FREE(asub); SUPERLU_FREE(xa);
    SUPERLU_FREE(a_save); SUPERLU_FREE(asub_save); SUPERLU_FREE(xa_save);
#endif
    if ( lwork > 0 ) {
	SUPERLU_FREE (work);
	Destroy_SuperMatrix_Store(&L);
	Destroy_SuperMatrix_Store(&U);
    }

    return (nfail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
} /* end main */

/*!
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
static void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *n, int *w, int *relax, int *nrhs, int *maxsuper,
		   int *rowblk, int *colblk, int_t *lwork, double *u, FILE **fp)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "ht:n:w:r:s:m:b:c:l:u:f:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-r <int> - granularity of relaxed supernodes\n");
	    exit(1);
	    break;
	  case 't': strcpy(matrix_type, optarg);
	            break;
	  case 'n': *n = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'r': *relax = atoi(optarg); 
	            break;
	  case 's': *nrhs = atoi(optarg); 
	            break;
	  case 'm': *maxsuper = atoi(optarg); 
	            break;
	  case 'b': *rowblk = atoi(optarg); 
	            break;
	  case 'c': *colblk = atoi(optarg); 
	            break;
	  case 'l': *lwork = atoi(optarg); 
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
          case 'f':
                    if ( !(*fp = fopen(optarg, "r")) ) {
                        ABORT("File does not exist");
                    }
                    printf(".. test sparse matrix in file: %s\n", optarg);
                    break;
  	}
    }
}
