/* ITSOL COPYRIGHT

Copyright (C) 2006, the University of Minnesota 

ITSOL is free software; you can redistribute it and/or modify it under
the terms of  the GNU General Public License as  published by the Free
Software Foundation [version 2 of the License, or any later version]
For details, see 

http://www.gnu.org/licenses/gpl-2.0.txt

A copy of the GNU licencing agreement is attached to the ITSOL package
in the file GNU.  For additional information contact the Free Software
Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 

DISCLAIMER
----------

This program  is distributed in the  hope that it will  be useful, but
WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
General Public License for more details. 

For information on ITSOL contact saad@cs.umn.edu
*/


/*! \file
 * \brief Flexible GMRES from ITSOL developed by Yousef Saad.
 *
 * \ingroup Example
 */

#include "slu_sdefs.h"

#define  epsmac  1.0e-16

extern float sdot_(int *, float [], int *, float [], int *);
extern float snrm2_(int *, float [], int *);

/*!
 * \brief Simple version of the ARMS preconditioned FGMRES algorithm.
 *
 *  Y. S. Dec. 2000. -- Apr. 2008
 *
 *  internal work arrays:
 *  vv      = work array of length [im+1][n] (used to store the Arnoldi
 *            basis)
 *  hh      = work array of length [im][im+1] (Householder matrix)
 *  z       = work array of length [im][n] to store preconditioned vectors
 *
 * \param [in] n         Dimension of vectors and matrices.
 * \param [in] smatvec   Operation for matrix-vector multiplication.
 * \param [in] spsolve   (right) preconditioning operation. Can be a NULL pointer (GMRES without preconditioner)
 * \param [in] rhs       Real vector of length n containing the right hand side.
 * \param [in,out] sol   In: Real vector of length n containing an initial guess to the solution on input.
 *                       Out: Contains an approximate solution (upon successful return).
 * \param [in] tol       Tolerance for stopping iteration
 * \param [in] im        Krylov subspace dimension
 * \param [in,out] itmax In: max number of iterations allowed.
 *                       Out: number of steps required to converge.
 * \param [in] fits      If NULL, no output. If not NULL, file handle to output "resid vs time and its".
 * \return Whether the algorithm finished successfully.
 */
int sfgmr(int n,
     void (*smatvec) (float, float[], float, float[]),
     void (*spsolve) (int, float[], float[]),
     float *rhs, float *sol, double tol, int im, int *itmax, FILE * fits)
{
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***
+-----------------------------------------------------------------------
| This is a simple version of the ARMS preconditioned FGMRES algorithm.
+-----------------------------------------------------------------------
| Y. S. Dec. 2000. -- Apr. 2008
+-----------------------------------------------------------------------
| on entry:
|----------
|
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension
| (itmax) = max number of iterations allowed.
| fits    = NULL: no output
|        != NULL: file handle to output " resid vs time and its"
|
| on return:
|----------
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| itmax   = has changed. It now contains the number of steps required
|           to converge --
+-----------------------------------------------------------------------
| internal work arrays:
|----------
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Householder matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| matvec - matrix-vector multiplication operation
| psolve - (right) preconditioning operation
|	   psolve can be a NULL pointer (GMRES without preconditioner)
+---------------------------------------------------------------------*/

    int maxits = *itmax;
    int its, i_1 = 1, i_2 = 2;
    float eps1 = 0.0;
    float **hh, *c, *s, *rs;
    float **vv, **z;
    float zero = 0.0;
    float one = 1.0;

    its = 0;
    vv = (float **)SUPERLU_MALLOC((im + 1) * sizeof(float *));
    for (int i = 0; i <= im; i++) vv[i] = floatMalloc(n);
    z = (float **)SUPERLU_MALLOC(im * sizeof(float *));
    hh = (float **)SUPERLU_MALLOC(im * sizeof(float *));
    for (int i = 0; i < im; i++)
    {
	hh[i] = floatMalloc(i + 2);
	z[i] = floatMalloc(n);
    }
    c = floatMalloc(im);
    s = floatMalloc(im);
    rs = floatMalloc(im + 1);

    /*---- outer loop starts here ----*/
    do
    {
	/*---- compute initial residual vector ----*/
	smatvec(one, sol, zero, vv[0]);
	for (int j = 0; j < n; j++)
	    vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
	float beta = snrm2_(&n, vv[0], &i_1);

	/*---- print info if fits != null ----*/
	if (fits != NULL && its == 0)
	    fprintf(fits, "%8d   %10.2e\n", its, beta);
	/*if ( beta <= tol * dnrm2_(&n, rhs, &i_1) )*/
	if ( !(beta > tol * snrm2_(&n, rhs, &i_1)) )
	    break;
	float t = 1.0 / beta;

	/*---- normalize: vv[0] = vv[0] / beta ----*/
	for (int j = 0; j < n; j++)
	    vv[0][j] = vv[0][j] * t;
	if (its == 0)
	    eps1 = tol * beta;

	/*---- initialize 1-st term of rhs of hessenberg system ----*/
	rs[0] = beta;
	int i = 0;
	for (i = 0; i < im; i++)
	{
	    its++;
	    int i1 = i + 1;

	    /*------------------------------------------------------------
	    |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
	    +-----------------------------------------------------------*/
	    if (spsolve)
		spsolve(n, z[i], vv[i]);
	    else
		scopy_(&n, vv[i], &i_1, z[i], &i_1);

	    /*---- matvec operation w = A z_{j} = A M^{-1} v_{j} ----*/
	    smatvec(one, z[i], zero, vv[i1]);

	    /*------------------------------------------------------------
	    |     modified gram - schmidt...
	    |     h_{i,j} = (w,v_{i})
	    |     w  = w - h_{i,j} v_{i}
	    +------------------------------------------------------------*/
	    float t0 = snrm2_(&n, vv[i1], &i_1);
	    for (int j = 0; j <= i; j++)
	    {
		float negt;
		float tt = sdot_(&n, vv[j], &i_1, vv[i1], &i_1);
		hh[i][j] = tt;
		negt = -tt;
		saxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
	    }

	    /*---- h_{j+1,j} = ||w||_{2} ----*/
	    t = snrm2_(&n, vv[i1], &i_1);
	    while (t < 0.5 * t0)
	    {
		t0 = t;
		for (int j = 0; j <= i; j++)
		{
		    float negt;
		    float tt = sdot_(&n, vv[j], &i_1, vv[i1], &i_1);
		    hh[i][j] += tt;
		    negt = -tt;
		    saxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
		}
		t = snrm2_(&n, vv[i1], &i_1);
	    }

	    hh[i][i1] = t;

	    if (t != 0.0)
	    {
		/*---- v_{j+1} = w / h_{j+1,j} ----*/
		t = 1.0 / t;
		for (int k = 0; k < n; k++)
		    vv[i1][k] = vv[i1][k] * t;
	    }
	    /*---------------------------------------------------
	    |     done with modified gram schmidt and arnoldi step
	    |     now  update factorization of hh
	    +--------------------------------------------------*/

	    /*--------------------------------------------------------
	    |   perform previous transformations  on i-th column of h
	    +-------------------------------------------------------*/
	    for (int k = 1; k <= i; k++)
	    {
		int k1 = k - 1;
		float tt = hh[i][k1];
		hh[i][k1] = c[k1] * tt + s[k1] * hh[i][k];
		hh[i][k] = -s[k1] * tt + c[k1] * hh[i][k];
	    }

	    float gam = sqrt(pow(hh[i][i], 2) + pow(hh[i][i1], 2));

	    /*---------------------------------------------------
	    |     if gamma is zero then any small value will do
	    |     affect only residual estimate
	    +--------------------------------------------------*/
	    /* if (gam == 0.0) gam = epsmac; */

	    /*---- get next plane rotation ---*/
	    if (gam == 0.0)
	    {
		c[i] = one;
		s[i] = zero;
	    }
            else
	    {
		c[i] = hh[i][i] / gam;
		s[i] = hh[i][i1] / gam;
	    }

	    rs[i1] = -s[i] * rs[i];
	    rs[i] = c[i] * rs[i];

	    /*----------------------------------------------------
	    |   determine residual norm and test for convergence
	    +---------------------------------------------------*/
	    hh[i][i] = c[i] * hh[i][i] + s[i] * hh[i][i1];
	    beta = fabs(rs[i1]);
	    if (fits != NULL)
		fprintf(fits, "%8d   %10.2e\n", its, beta);
	    if (beta <= eps1 || its >= maxits)
		break;
	}

	if (i == im) i--;

	/*---- now compute solution. 1st, solve upper triangular system ----*/
	rs[i] = rs[i] / hh[i][i];

	for (int ii = 1; ii <= i; ii++)
	{
	    int k = i - ii;
	    float tt = rs[k];
	    for (int j = k + 1; j <= i; j++)
		tt = tt - hh[j][k] * rs[j];
	    rs[k] = tt / hh[k][k];
	}

	/*---- linear combination of v[i]'s to get sol. ----*/
	for (int j = 0; j <= i; j++)
	{
	    float tt = rs[j];
	    for (int k = 0; k < n; k++)
		sol[k] += tt * z[j][k];
	}

	/* calculate the residual and output */
	smatvec(one, sol, zero, vv[0]);
	for (int j = 0; j < n; j++)
	    vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */

	/*---- print info if fits != null ----*/
	beta = snrm2_(&n, vv[0], &i_1);

	/*---- restart outer loop if needed ----*/
	/*if (beta >= eps1 / tol)*/
	if ( !(beta < eps1 / tol) )
	{
	    its = maxits + 10;
	    break;
	}
	if (beta <= eps1)
	    break;
    } while(its < maxits);

    int retval = (its >= maxits);
    for (int i = 0; i <= im; i++)
	SUPERLU_FREE(vv[i]);
    SUPERLU_FREE(vv);
    for (int i = 0; i < im; i++)
    {
	SUPERLU_FREE(hh[i]);
	SUPERLU_FREE(z[i]);
    }
    SUPERLU_FREE(hh);
    SUPERLU_FREE(z);
    SUPERLU_FREE(c);
    SUPERLU_FREE(s);
    SUPERLU_FREE(rs);

    *itmax = its;

    return retval;
} /*----end of fgmr ----*/
