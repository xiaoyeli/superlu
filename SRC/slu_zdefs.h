/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_zdefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
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
 *	after the numercial LU factorization phase.
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
#ifndef __SUPERLU_zSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_zSP_DEFS

/*
 * File name:		zsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

/* Define my integer type int_t */
#ifdef _LONGINT
typedef long long int int_t;
#define IFMT "%lld"
#else
typedef int int_t; /* default */
#define IFMT "%8d"
#endif

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"
#include "slu_dcomplex.h"


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
zgssv(superlu_options_t *, SuperMatrix *, int_t *, int_t *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int_t *);
extern void
zgssvx(superlu_options_t *, SuperMatrix *, int_t *, int_t *, int_t *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int_t, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *);
    /* ILU */
extern void
zgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
zgsisx(superlu_options_t *, SuperMatrix *, int_t *, int_t *, int_t *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int_t, SuperMatrix *, SuperMatrix *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *);


/*! \brief Supernodal LU factor related */
extern void
zCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, doublecomplex *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_CompRow_Matrix(SuperMatrix *, int_t, int_t, int_t, doublecomplex *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, doublecomplex *, int_t,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, int_t, int_t, int_t, doublecomplex *, 
		         int_t *, int_t *, int_t *, int_t *, int_t *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_Dense_Matrix(int_t, int_t, doublecomplex *, int_t, doublecomplex *, int_t);

extern void    countnz (const int_t, int_t *, int_t *, int_t *, GlobalLU_t *);
extern void    ilu_countnz (const int_t, int_t *, int_t *, GlobalLU_t *);
extern void    fixupL (const int_t, const int_t *, GlobalLU_t *);

extern void    zallocateA (int_t, int_t, doublecomplex **, int_t **, int_t **);
extern void    zgstrf (superlu_options_t*, SuperMatrix*,
                       int_t, int_t, int_t*, void *, int_t, int_t *, int_t *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int_t *);
extern int_t     zsnode_dfs (const int_t, const int_t, const int_t *, const int_t *,
			     const int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     zsnode_bmod (const int_t, const int_t, const int_t, doublecomplex *,
                              doublecomplex *, GlobalLU_t *, SuperLUStat_t*);
extern void    zpanel_dfs (const int_t, const int_t, const int_t, SuperMatrix *,
			   int_t *, int_t *, doublecomplex *, int_t *, int_t *, int_t *,
			   int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern void    zpanel_bmod (const int_t, const int_t, const int_t, const int_t,
                           doublecomplex *, doublecomplex *, int_t *, int_t *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int_t     zcolumn_dfs (const int_t, const int_t, int_t *, int_t *, int_t *, int_t *,
			   int_t *, int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     zcolumn_bmod (const int_t, const int_t, doublecomplex *,
			   doublecomplex *, int_t *, int_t *, int_t,
                           GlobalLU_t *, SuperLUStat_t*);
extern int_t     zcopy_to_ucol (int_t, int_t, int_t *, int_t *, int_t *,
                              doublecomplex *, GlobalLU_t *);         
extern int_t     zpivotL (const int_t, const double, int_t *, int_t *, 
                         int_t *, int_t *, int_t *, GlobalLU_t *, SuperLUStat_t*);
extern void    zpruneL (const int_t, const int_t *, const int_t, const int_t,
			  const int_t *, const int_t *, int_t *, GlobalLU_t *);
extern void    zreadmt (int *, int *, int *, doublecomplex **, int **, int **);
extern void    zGenXtrue (int, int, doublecomplex *, int);
extern void    zFillRHS (trans_t, int, doublecomplex *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    zgstrs (trans_t, SuperMatrix *, SuperMatrix *, int_t *, int_t *,
                        SuperMatrix *, SuperLUStat_t*, int_t *);
/* ILU */
extern void    zgsitrf (superlu_options_t*, SuperMatrix*, int_t, int_t, int_t*,
		        void *, int_t, int_t *, int_t *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int_t *);
extern int     zldperm(int_t, int_t, int_t, int_t [], int_t [], doublecomplex [],
                        int_t [],	double [], double []);
extern int_t     ilu_zsnode_dfs (const int_t, const int_t, const int_t *, const int_t *,
			       const int_t *, int_t *, GlobalLU_t *);
extern void    ilu_zpanel_dfs (const int_t, const int_t, const int_t, SuperMatrix *,
			       int_t *, int_t *, doublecomplex *, double *, int_t *, int_t *,
			       int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     ilu_zcolumn_dfs (const int_t, const int_t, int_t *, int_t *, int_t *,
				int_t *, int_t *, int_t *, int_t *, int_t *,
				GlobalLU_t *);
extern int_t     ilu_zcopy_to_ucol (int_t, int_t, int_t *, int_t *, int_t *,
                                  doublecomplex *, int_t, milu_t, double, int_t,
                                  doublecomplex *, int_t *, GlobalLU_t *, double *);
extern int_t     ilu_zpivotL (const int_t, const double, int_t *, int_t *, int_t, int_t *,
			    int_t *, int_t *, int_t *, double, milu_t,
                            doublecomplex, GlobalLU_t *, SuperLUStat_t*);
extern int_t     ilu_zdrop_row (superlu_options_t *, int_t, int_t, double,
                              int_t, int_t *, double *, GlobalLU_t *, 
                              double *, double *, int_t);


/*! \brief Driver related */

extern void    zgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, int_t *);
extern void    zlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
extern void    zgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, int_t *);
extern double   zPivotGrowth(int_t, SuperMatrix *, int_t *, 
                            SuperMatrix *, SuperMatrix *);
extern void    zgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int_t *, int_t *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, int_t *);

extern int     sp_ztrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, doublecomplex *, SuperLUStat_t*, int *);
extern int     sp_zgemv (char *, doublecomplex, SuperMatrix *, doublecomplex *,
			int, doublecomplex, doublecomplex *, int);

extern int     sp_zgemm (char *, char *, int, int, int, doublecomplex,
			SuperMatrix *, doublecomplex *, int, doublecomplex, 
			doublecomplex *, int);
extern         double dmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int_t     zLUMemInit (fact_t, void *, int_t, int_t, int_t, int_t, int_t,
                            double, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int_t **, doublecomplex **);
extern void    zSetRWork (int_t, int_t, doublecomplex *, doublecomplex **, doublecomplex **);
extern void    zLUWorkFree (int_t *, doublecomplex *, GlobalLU_t *);
extern int_t     zLUMemXpand (int_t, int_t, MemType, int_t *, GlobalLU_t *);

extern doublecomplex  *doublecomplexMalloc(int_t);
extern doublecomplex  *doublecomplexCalloc(int_t);
extern double  *doubleMalloc(int_t);
extern double  *doubleCalloc(int_t);
extern int_t     zmemory_usage(const int_t, const int_t, const int_t, const int_t);
extern int_t     zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern int_t     ilu_zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    zreadhb(FILE *, int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void    zreadrb(int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void    zreadtriple(int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void    zreadMM(FILE *, int_t *, int_t *, int_t *, doublecomplex **, int_t **, int_t **);
extern void    zCompRow_to_CompCol(int_t, int_t, int_t, doublecomplex*, int_t*, int_t*,
		                   doublecomplex **, int_t **, int_t **);
extern void    zfill (doublecomplex *, int, doublecomplex);
extern void    zinf_norm_error (int, SuperMatrix *, doublecomplex *);
extern double  dqselect(int, double *, int);


/*! \brief Routines for debugging */
extern void    zPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    zPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    zPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    zprint_lu_col(char *, int, int, int *, GlobalLU_t *);
extern int     print_double_vec(char *, int, double *);
extern void    zcheck_tempv(int, doublecomplex *);

/*! \brief BLAS */

extern int zgemm_(const char*, const char*, const int*, const int*, const int*,
                  const doublecomplex*, const doublecomplex*, const int*, const doublecomplex*,
		  const int*, const doublecomplex*, doublecomplex*, const int*);
extern int ztrsv_(char*, char*, char*, int*, doublecomplex*, int*,
                  doublecomplex*, int*);
extern int ztrsm_(char*, char*, char*, char*, int*, int*,
                  doublecomplex*, doublecomplex*, int*, doublecomplex*, int*);
extern int zgemv_(char *, int *, int *, doublecomplex *, doublecomplex *a, int *,
                  doublecomplex *, int *, doublecomplex *, doublecomplex *, int *);

extern void zusolve(int, int, doublecomplex*, doublecomplex*);
extern void zlsolve(int, int, doublecomplex*, doublecomplex*);
extern void zmatvec(int, int, int, doublecomplex*, doublecomplex*, doublecomplex*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_zSP_DEFS */

