/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 5.0) --
 * Lawrence Berkeley National Laboratory
 * November, 2010
 * August, 2011
 */

/*! \file
 * \brief Example #1 showing how to use ILU to precondition GMRES
 *
 * This example shows that ILU is computed from the equilibrated matrix,
 * and the preconditioned GMRES is applied to the equilibrated system.
 * The driver routine DGSISX is called twice to perform factorization
 * and apply preconditioner separately.
 * 
 * Note that DGSISX performs the following factorization:
 *     Pr*Dr*A*Dc*Pc^T ~= LU
 * with Pr being obtained from MC64 statically then partial pivoting
 * dynamically. On return, A is overwritten as A1 = Dr*A*Dc.
 *
 * We can solve the transformed system, A1*y = Dr*B, using FGMRES.
 * B is first overwritten as Dr*B.
 * Then GMRES step requires requires 2 procedures:
 *   1) Apply preconditioner M^{-1} = Pc^T*U^{-1}*L^{-1}*Pr
 *   2) Matrix-vector multiplication: w = A1*v
 * 
 * \ingroup Example
 */

#include "slu_ddefs.h"

superlu_options_t *GLOBAL_OPTIONS;
double *GLOBAL_R, *GLOBAL_C;
int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
SuperLUStat_t *GLOBAL_STAT;
mem_usage_t   *GLOBAL_MEM_USAGE;

/*!
 * \brief Performs dgsisx with original matrix A.
 *
 * See documentation of dgsisx for more details.
 *
 * \param [in] n     Dimension of matrices
 * \param [out] x    Solution
 * \param [in,out] y Right-hand side
 */
void dpsolve(int n, double x[], double y[])
{
    SuperMatrix *A = GLOBAL_A, *L = GLOBAL_L, *U = GLOBAL_U;
    SuperLUStat_t *stat = GLOBAL_STAT;
    int *perm_c = GLOBAL_PERM_C, *perm_r = GLOBAL_PERM_R;
    char equed[1] = {'N'};
    double *R = GLOBAL_R, *C = GLOBAL_C;
    superlu_options_t *options = GLOBAL_OPTIONS;
    mem_usage_t  *mem_usage = GLOBAL_MEM_USAGE;
    int_t info;
    static DNformat X, Y;
    static SuperMatrix XX = {SLU_DN, SLU_D, SLU_GE, 1, 1, &X};
    static SuperMatrix YY = {SLU_DN, SLU_D, SLU_GE, 1, 1, &Y};
    double rpg, rcond;

    XX.nrow = YY.nrow = n;
    X.lda = Y.lda = n;
    X.nzval = x;
    Y.nzval = y;

#if 0
    dcopy_(&n, y, &i_1, x, &i_1);
    dgstrs(NOTRANS, L, U, perm_c, perm_r, &XX, stat, &info);
#else
    dgsisx(options, A, perm_c, perm_r, NULL, equed, R, C,
	   L, U, NULL, 0, &YY, &XX, &rpg, &rcond, NULL,
	   mem_usage, stat, &info);
#endif
}


/*!
 * \brief Performs matrix-vector multiplication sp_dgemv with original matrix A.
 *
 * The operations is y := alpha*A*x + beta*y. See documentation of sp_dgemv
 * for further details.
 *
 * \param [in] alpha Scalar factor for A*x
 * \param [in] x Vector to multiply with A
 * \param [in] beta Scalar factor for y
 * \param [in,out] y Vector to add to to matrix-vector multiplication and
 *                   storage for result.
 */
void dmatvec_mult(double alpha, double x[], double beta, double y[])
{
    SuperMatrix *A = GLOBAL_A;

    sp_dgemv("N", alpha, A, x, 1, beta, y, 1);
}

int main(int argc, char *argv[])
{
    void dmatvec_mult(double alpha, double x[], double beta, double y[]);
    void dpsolve(int n, double x[], double y[]);
    extern int dfgmr( int n,
	void (*matvec_mult)(double, double [], double, double []),
	void (*psolve)(int n, double [], double[]),
	double *rhs, double *sol, double tol, int restrt, int *itmax,
	FILE *fits);
    extern int dfill_diag(int n, NCformat *Astore);

    char     equed[1] = {'B'};
    trans_t  trans;
    SuperMatrix A, L, U;
    SuperMatrix B, X;
    NCformat *Astore;
    NCformat *Ustore;
    SCformat *Lstore;
    GlobalLU_t	   Glu; /* facilitate multiple factorizations with 
                           SamePattern_SameRowPerm                  */
    double   *a;
    int_t    *asub, *xa;
    int      *etree;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    int      nrhs, ldx, m, n;
    int_t    info, nnz, lwork;
    double   *rhsb, *rhsx, *xact;
    double   *work = NULL;
    double   *R, *C;
    double   rpg, rcond;
    double zero = 0.0;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE 	  *fp = stdin;

    double *x, *b;

#ifdef DEBUG
    extern int num_drop_L, num_drop_U;
#endif

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    trans = NOTRANS;

    /* Set the default input options:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 0.1; //different from complete LU
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
	options.RowPerm = LargeDiag_MC64;
	options.ILU_DropTol = 1e-4;
	options.ILU_FillTol = 1e-2;
	options.ILU_FillFactor = 10.0;
	options.ILU_DropRule = DROP_BASIC | DROP_AREA;
	options.ILU_Norm = INF_NORM;
	options.ILU_MILU = SILU;
     */
    ilu_set_default_options(&options);

    /* Modify the defaults. */
    options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
    options.ConditionNumber = YES;/* Compute reciprocal condition number */

    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) ABORT("Malloc fails for work[].");
    }

    /* Read matrix A from a file in Harwell-Boeing format.*/
    if (argc < 2)
    {
	printf("Usage:\n%s [OPTION] < [INPUT] > [OUTPUT]\nOPTION:\n"
		"-h -hb:\n\t[INPUT] is a Harwell-Boeing format matrix.\n"
		"-r -rb:\n\t[INPUT] is a Rutherford-Boeing format matrix.\n"
		"-t -triplet:\n\t[INPUT] is a triplet format matrix.\n",
		argv[0]);
        return EXIT_FAILURE;
    }
    else
    {
	switch (argv[1][1])
	{
	    case 'H':
	    case 'h':
		printf("Input a Harwell-Boeing format matrix:\n");
		dreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
		break;
	    case 'R':
	    case 'r':
		printf("Input a Rutherford-Boeing format matrix:\n");
		dreadrb(&m, &n, &nnz, &a, &asub, &xa);
		break;
	    case 'T':
	    case 't':
		printf("Input a triplet format matrix:\n");
		dreadtriple(&m, &n, &nnz, &a, &asub, &xa);
		break;
	    default:
		printf("Unrecognized format.\n");
		return EXIT_FAILURE;
	}
    }

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa,
                                SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    dfill_diag(n, Astore);
    printf("Dimension %dx%d; # nonzeros %d\n", (int)A.nrow, (int)A.ncol, (int)Astore->nnz);
    fflush(stdout);

    /* Generate the right-hand side */
    if ( !(rhsb = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if ( !(etree = int32Malloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = int32Malloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = int32Malloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
	ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
	ABORT("SUPERLU_MALLOC fails for C[].");

    info = 0;
#ifdef DEBUG
    num_drop_L = 0;
    num_drop_U = 0;
#endif

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Compute the incomplete factorization and compute the condition number
       and pivot growth using dgsisx. */
    B.ncol = 0;  /* not to perform triangular solution */
    dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work,
	   lwork, &B, &X, &rpg, &rcond, &Glu, &mem_usage, &stat, &info);

    /* Set RHS for GMRES. */
    if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
    if (*equed == 'R' || *equed == 'B') {
	for (int i = 0; i < n; ++i) b[i] = rhsb[i] * R[i];
    } else {
	for (int i = 0; i < m; i++) b[i] = rhsb[i];
    }

    printf("dgsisx(): info %lld, equed %c\n", (long long)info, equed[0]);
    if (info > 0 || rcond < 1e-8 || rpg > 1e8)
	printf("WARNING: This preconditioner might be unstable.\n");

    if ( info == 0 || info == n+1 ) {
	if ( options.PivotGrowth == YES )
	    printf("Recip. pivot growth = %e\n", rpg);
	if ( options.ConditionNumber == YES )
	    printf("Recip. condition number = %e\n", rcond);
    } else if ( info > 0 && lwork == -1 ) {
	printf("** Estimated memory: %lld bytes\n", (long long)info - n);
    }

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;
    printf("n(A) = %d, nnz(A) = %lld\n", n, (long long) Astore->nnz);
    printf("No of nonzeros in factor L = %lld\n", (long long) Lstore->nnz);
    printf("No of nonzeros in factor U = %lld\n", (long long) Ustore->nnz);
    printf("No of nonzeros in L+U = %lld\n", (long long) Lstore->nnz + Ustore->nnz - n);
    printf("Fill ratio: nnz(F)/nnz(A) = %.1f\n",
	    ((double)(Lstore->nnz) + (double)(Ustore->nnz) - (double)n)
	    / (double)Astore->nnz);
    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
    fflush(stdout);

    /* Set the global variables. */
    GLOBAL_A = &A;
    GLOBAL_L = &L;
    GLOBAL_U = &U;
    GLOBAL_STAT = &stat;
    GLOBAL_PERM_C = perm_c;
    GLOBAL_PERM_R = perm_r;
    GLOBAL_OPTIONS = &options;
    GLOBAL_R = R;
    GLOBAL_C = C;
    GLOBAL_MEM_USAGE = &mem_usage;

    /* Set the options to do solve-only. */
    options.Fact = FACTORED;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;

    /* Set the variables used by GMRES. */
    int restrt = SUPERLU_MIN(n / 3 + 1, 50);
    int maxit = 1000;
    int iter = maxit;
    double resid = 1e-8;
    if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");

    if (info <= n + 1)
    {
	int i_1 = 1, nnz32;
	double maxferr = 0.0, nrmA, nrmB, res, t;
	extern double dnrm2_(int *, double [], int *);
	extern void daxpy_(int *, double *, double [], int *, double [], int *);

	/* Initial guess */
	for (int i = 0; i < n; i++) x[i] = zero;

	t = SuperLU_timer_();

	/* Call GMRES */
	dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, stdout);

	t = SuperLU_timer_() - t;

	/* Output the result. */
	nnz32 = Astore->nnz;
	nrmA = dnrm2_(&nnz32, (double *)((NCformat *)A.Store)->nzval,
		&i_1);
	nrmB = dnrm2_(&m, b, &i_1);
	sp_dgemv("N", -1.0, &A, x, 1, 1.0, b, 1);
	res = dnrm2_(&m, b, &i_1);
	resid = res / nrmB;
	printf("||A||_F = %.1e, ||B||_2 = %.1e, ||B-A*X||_2 = %.1e, "
		"relres = %.1e\n", nrmA, nrmB, res, resid);

	if (iter >= maxit)
	{
	    if (resid >= 1.0) iter = -180;
	    else if (resid > 1e-8) iter = -111;
	}
	printf("iteration: %d\nresidual: %.1e\nGMRES time: %.2f seconds.\n",
		iter, resid, t);

	/* Scale the solution back if equilibration was performed. */
	if (*equed == 'C' || *equed == 'B') 
	    for (int i = 0; i < n; i++) x[i] *= C[i];

	for (int i = 0; i < m; i++) {
	    maxferr = SUPERLU_MAX(maxferr, fabs(x[i] - xact[i]));
        }
	printf("||X-X_true||_oo = %.1e\n", maxferr);
    }
#ifdef DEBUG
    printf("%d entries in L and %d entries in U dropped.\n",
	    num_drop_L, num_drop_U);
#endif
    fflush(stdout);

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork >= 0 ) {
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
    }
    SUPERLU_FREE(b);
    SUPERLU_FREE(x);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif

    return EXIT_SUCCESS;
}
