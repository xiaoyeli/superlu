/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_ddefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 7.0.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * August 2024
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: \#supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numerical LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_dSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_dSP_DEFS

/*
 * File name:		dsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "superlu_config.h"
#include "supermatrix.h"
#include "slu_util.h"


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
dgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int_t *info);
extern void
dgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int_t lwork, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *info);
    /* ILU */
extern void
dgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
dgsisx(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
       int *etree, char *equed, double *R, double *C,
       SuperMatrix *L, SuperMatrix *U, void *work, int_t lwork,
       SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth, double *rcond,
       GlobalLU_t *Glu, mem_usage_t *mem_usage, SuperLUStat_t *stat, int_t *info);


/*! \brief Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int_t, double *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompRow_Matrix(SuperMatrix *, int, int, int_t, double *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void dCompRow_to_CompCol(int, int, int_t, double*, int_t*, int_t*,
		                   double **, int_t **, int_t **);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, int, int, int_t, double *, 
		         int_t *, int_t *, int_t *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix(int, int, double *, int, double *, int);

extern void    dallocateA (int, int_t, double **, int_t **, int_t **);
extern void    dgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int_t, int *, int *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int_t *info);
extern int_t   dsnode_dfs (const int, const int, const int_t *, const int_t *,
			     const int_t *, int_t *, int *, GlobalLU_t *);
extern int     dsnode_bmod (const int, const int, const int, double *,
                              double *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, double *, int *, int *, int *,
			   int_t *, int *, int *, int_t *, GlobalLU_t *);
extern void    dpanel_bmod (const int, const int, const int, const int,
                           double *, double *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int     dcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int_t *, int *, int *, int_t *, GlobalLU_t *);
extern int     dcolumn_bmod (const int, const int, double *,
			   double *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
extern int     dcopy_to_ucol (int, int, int *, int *, int *,
                              double *, GlobalLU_t *);         
extern int     dpivotL (const int, const double, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpruneL (const int, const int *, const int, const int,
			  const int *, const int *, int_t *, GlobalLU_t *);
extern void    dreadmt (int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dGenXtrue (int, int, double *, int);
extern void    dFillRHS (trans_t, int, double *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, const int *, const int *,
                        SuperMatrix *, SuperLUStat_t*, int *);
/* ILU */
extern void    dgsitrf (superlu_options_t*, SuperMatrix*, int, int, int*,
		        void *, int_t, int *, int *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int_t *info);
extern int     dldperm(int, int, int_t, int_t [], int_t [], double [],
                        int [],	double [], double []);
extern int     ilu_dsnode_dfs (const int, const int, const int_t *, const int_t *,
			       const int_t *, int *, GlobalLU_t *);
extern void    ilu_dpanel_dfs (const int, const int, const int, SuperMatrix *,
			       int *, int *, double *, double *, int *, int *,
			       int *, int *, int *, int_t *, GlobalLU_t *);
extern int     ilu_dcolumn_dfs (const int, const int, int *, int *, int *,
				int *, int *, int *, int *, int_t *, GlobalLU_t *);
extern int     ilu_dcopy_to_ucol (int, int, int *, int *, int *,
                                  double *, int, milu_t, double, int,
                                  double *, int *, GlobalLU_t *, double *);
extern int     ilu_dpivotL (const int, const double, int *, int *, int, int *,
			    int *, int *, int *, double, milu_t,
                            double, GlobalLU_t *, SuperLUStat_t*);
extern int     ilu_ddrop_row (superlu_options_t *, int, int, double,
                              int, int *, double *, GlobalLU_t *, 
                              double *, double *, int);


/*! \brief Driver related */

extern void    dgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, int *);
extern void    dlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
extern void    dgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, int *);
extern double   dPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    dgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, int *);

extern int     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, SuperLUStat_t*, int *);
extern int     sp_dgemv (char *, double, SuperMatrix *, double *,
			int, double, double *, int);

extern int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);
extern         double dmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int_t   dLUMemInit (fact_t, void *, int_t, int, int, int_t, int,
                            double, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int **, double **);
extern void    dSetRWork (int, int, double *, double **, double **);
extern void    dLUWorkFree (int *, double *, GlobalLU_t *);
extern int_t   dLUMemXpand (int, int_t, MemType, int_t *, GlobalLU_t *);

extern double  *doubleMalloc(size_t);
extern double  *doubleCalloc(size_t);
extern int_t   dmemory_usage(const int_t, const int_t, const int_t, const int);
extern int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern int     ilu_dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    dreadhb(FILE *, int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dreadrb(int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dreadtriple(int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dreadtriple_noheader(int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dreadMM(FILE *, int *, int *, int_t *, double **, int_t **, int_t **);
extern void    dfill (double *, int, double);
extern void    dinf_norm_error (int, SuperMatrix *, double *);
extern double  dqselect(int, double *, int);


/*! \brief Routines for debugging */
extern void    dPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    dPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    dPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    dprint_lu_col(char *, int, int, int_t *, GlobalLU_t *);
extern int     print_double_vec(const char *, int, const double *);
extern void    dcheck_tempv(int, double *);

/*! \brief BLAS */

extern void dcopy_(int *, double *, int *, double *, int *);
extern void daxpy_(int *, double *, double *, int *, double *, int *);
extern void dgemm_(const char*, const char*, const int*, const int*, const int*,
                  const double*, const double*, const int*, const double*,
		  const int*, const double*, double*, const int*);
extern void dtrsv_(char*, char*, char*, int*, double*, int*,
                  double*, int*);
extern void dtrsm_(char*, char*, char*, char*, int*, int*,
                  double*, double*, int*, double*, int*);
extern void dgemv_(char *, int *, int *, double *, double *a, int *,
                  double *, int *, double *, double *, int *);

extern void dusolve(int, int, double*, double*);
extern void dlsolve(int, int, double*, double*);
extern void dmatvec(int, int, int, double*, double*, double*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_dSP_DEFS */

