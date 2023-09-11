/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! \file
 * \brief Read matrix in triplet format.
 *
 * \ingroup Example
 */

#include <stdio.h>
#include "slu_ddefs.h"

#undef EXPAND_SYM

/*!
 * \brief Read matrix in triplet format from stdin.
 *
 * File format: Triplet in a line for each nonzero entry:
 * \code row col value \endcode or
 * \code row col real_part imaginary_part \endcode
 *
 * \param [out] m      Number of rows in the matrix
 * \param [out] n      Number of columns in the matrix
 * \param [out] nonz   Number of non-zero entries in the matrix
 * \param [out] rowind Contains the row subscripts of nonzeros in columns of matrix A
 * \param [out] nzval  The numerical values
 * \param [out] colptr Column i of A is given by (*nzval)[k], k = (*rowind)[i],...,(*rowind)[i+1]-1.
 */
void
dreadtriple_noheader(int *m, int *n, int_t *nonz,
                     double **nzval, int_t **rowind, int_t **colptr)
{
    int    i, j, jsize, minn = 100;
    double *a, *val, vali;
    int_t  *asub, *xa, k, nnz, nz, new_nonz;
    int    *row, *col;
    int    zero_base = 0, ret_val = 0;
    FILE *fp = stdin;

    /* First pass: determine N and NNZ */
    nz = *n = 0;

    ret_val = fscanf(fp, "%d%d%lf\n", &i, &j, &vali);
    printf("%d\t%d\t %f\n", i, j, vali);

    while (ret_val != EOF) {
	*n = SUPERLU_MAX(*n, i);
	*n = SUPERLU_MAX(*n, j);
	minn = SUPERLU_MIN(minn, i);
	minn = SUPERLU_MIN(minn, j);
	++nz;

        ret_val = fscanf(fp, "%d%d%lf\n", &i, &j, &vali);
    }
    
    if ( minn == 0 ) { /* zero-based indexing */
	zero_base = 1;
	++(*n);
	printf("triplet file: row/col indices are zero-based.\n");
    } else {
	printf("triplet file: row/col indices are one-based.\n");
    }

    *m = *n;
    *nonz = nz;
    rewind(fp);  // Move to the start of the input file 

#ifdef EXPAND_SYM
    new_nonz = 2 * *nonz - *n;
#else
    new_nonz = *nonz;
#endif

    /* Second pass: read the actual matrix values */
    printf("m %d, n %d, nonz %lld, new_nonz %lld\n", *m, *n,
	   (long long)*nonz, (long long)new_nonz);
    dallocateA(*n, new_nonz, nzval, rowind, colptr); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;

    if ( !(val = (double *) SUPERLU_MALLOC(new_nonz * sizeof(double))) )
        ABORT("Malloc fails for val[]");
    if ( !(row = (int *) SUPERLU_MALLOC(new_nonz * sizeof(int))) )
        ABORT("Malloc fails for row[]");
    if ( !(col = (int *) SUPERLU_MALLOC(new_nonz * sizeof(int))) )
        ABORT("Malloc fails for col[]");

    for (j = 0; j < *n; ++j) xa[j] = 0;

    /* Read into the triplet array from a file */
    for (nnz = 0, nz = 0; nnz < *nonz; ++nnz) {
	
	fscanf(fp, "%d%d%lf\n", &row[nz], &col[nz], &val[nz]);

	if ( !zero_base ) {
	    /* Change to 0-based indexing. */
	    --row[nz];
	    --col[nz];
	}

	if (row[nz] < 0 || row[nz] >= *m || col[nz] < 0 || col[nz] >= *n
	    /*|| val[nz] == 0.*/) {
	    fprintf(stderr, "nz %lld, (%d, %d) = %e out of bound, removed\n",
		    (long long)nz, row[nz], col[nz], val[nz]);
	    exit(-1);
	} else {
	    ++xa[col[nz]];
#ifdef EXPAND_SYM
	    if ( row[nz] != col[nz] ) { /* Excluding diagonal */
	      ++nz;
	      row[nz] = col[nz-1];
	      col[nz] = row[nz-1];
	      val[nz] = val[nz-1];
	      ++xa[col[nz]];
	    }
#endif
	    ++nz;
	}
    }

    *nonz = nz;
    
#ifdef EXPAND_SYM
    printf("new_nonz after symmetric expansion:\t%lld\n", (long long)*nonz);
#endif


    /* Initialize the array of column pointers */
    k = 0;
    jsize = xa[0];
    xa[0] = 0;
    for (j = 1; j < *n; ++j) {
	k += jsize;
	jsize = xa[j];
	xa[j] = k;
    }

    /* Copy the triplets into the column oriented storage */
    for (nz = 0; nz < *nonz; ++nz) {
	j = col[nz];
	k = xa[j];
	asub[k] = row[nz];
	a[k] = val[nz];
	++xa[j];
    }

    /* Reset the column pointers to the beginning of each column */
    for (j = *n; j > 0; --j)
	xa[j] = xa[j-1];
    xa[0] = 0;

    SUPERLU_FREE(val);
    SUPERLU_FREE(row);
    SUPERLU_FREE(col);

    
#ifdef CHK_INPUT
    for (i = 0; i < *n; i++) {
	printf("Col %d, xa[i] %d, xa[i+1] %d\n", i, xa[i], xa[i+1]);
	for (k = xa[i]; k < xa[i+1]; k++)
	    printf("%d\t%16.10f\n", asub[k], a[k]);
    }
    exit(-1);
#endif

}

#if 0
void dreadrhs(int m, double *b)
{
    FILE *fp, *fopen();
    int i, j;

    if ( !(fp = fopen("b.dat", "r")) ) {
        fprintf(stderr, "zreadrhs: file does not exist\n");
	exit(-1);
    }
    for (i = 0; i < m; ++i)
      fscanf(fp, "%lf\n", &b[i]);

    fclose(fp);
}
#endif
