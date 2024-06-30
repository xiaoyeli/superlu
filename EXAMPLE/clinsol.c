/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

/*! \file
 * \brief LU factorization from CGSTRF (CGSSV)
 *
 * \ingroup Example
 */

#include "slu_cdefs.h"

int main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    singlecomplex   *a;
    int_t    *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, m, n;
    int_t    nnz, info;
    singlecomplex   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE      *fp = stdin;
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

    /* Read the matrix in Harwell-Boeing format. */
    creadhb(fp, &m, &n, &nnz, &a, &asub, &xa);

    cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_C, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", (int)A.nrow, (int)A.ncol, (int)Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = singlecomplexMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    cCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_C, SLU_GE);
    xact = singlecomplexMalloc(n * nrhs);
    ldx = n;
    cGenXtrue(n, nrhs, xact, ldx);
    cFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = int32Malloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = int32Malloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    cgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        singlecomplex *sol = (singlecomplex*) ((DNformat*) B.Store)->nzval; 
        (void)sol;  // suppress unused variable warning

	 /* Compute the infinity norm of the error. */
	cinf_norm_error(nrhs, &B, xact);

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %lld\n", (long long) Lstore->nnz);
    	printf("No of nonzeros in factor U = %lld\n", (long long) Ustore->nnz);
    	printf("No of nonzeros in L+U = %lld\n", (long long) Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
	
	cQuerySpace(&L, &U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	
    } else {
	printf("cgssv() error returns INFO= %lld\n", (long long)info);
	if ( info <= n ) { /* factorization completes */
	    cQuerySpace(&L, &U, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	}
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif

    return EXIT_SUCCESS;
}

