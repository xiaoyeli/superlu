/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */

/*! \file
 * Convert a full matrix into a sparse matrix format.
 *
 * \ingroup TestingC
 */

#include "slu_cdefs.h"

/*!
 * Convert a full matrix into a sparse matrix format.
 *
 * For complex float.
 */
int
sp_cconvert(int m, int n, singlecomplex *A, int lda, int kl, int ku,
	   singlecomplex *a, int_t *asub, int_t *xa, int_t *nnz)
{
    int_t     lasta = 0;
    int_t     i, j, ilow, ihigh;
    int_t     *row;
    singlecomplex  *val;

    for (j = 0; j < n; ++j) {
	xa[j] = lasta;
	val = &a[xa[j]];
	row = &asub[xa[j]];

	ilow = SUPERLU_MAX(0, j - ku);
	ihigh = SUPERLU_MIN(n-1, j + kl);
	for (i = ilow; i <= ihigh; ++i) {
	    val[i-ilow] = A[i + j*lda];
	    row[i-ilow] = i;
	}
	lasta += ihigh - ilow + 1;
    }

    xa[n] = *nnz = lasta;
    return 0;
}


