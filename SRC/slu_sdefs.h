/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_sdefs.h
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
#ifndef __SUPERLU_sSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_sSP_DEFS

/*
 * File name:		ssp_defs.h
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


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
sgssv(superlu_options_t *, SuperMatrix *, int_t *, int_t *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int_t *);
extern void
sgssvx(superlu_options_t *, SuperMatrix *, int_t *, int_t *, int_t *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int_t, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *);
    /* ILU */
extern void
sgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
sgsisx(superlu_options_t *, SuperMatrix *, int_t *, int_t *, int_t *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int_t, SuperMatrix *, SuperMatrix *, float *, float *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *);


/*! \brief Supernodal LU factor related */
extern void
sCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, float *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_CompRow_Matrix(SuperMatrix *, int_t, int_t, int_t, float *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
sCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, float *, int_t,
		     Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_SuperNode_Matrix(SuperMatrix *, int_t, int_t, int_t, float *, 
		         int_t *, int_t *, int_t *, int_t *, int_t *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_Dense_Matrix(int_t, int_t, float *, int_t, float *, int_t);

extern void    countnz (const int_t, int_t *, int_t *, int_t *, GlobalLU_t *);
extern void    ilu_countnz (const int_t, int_t *, int_t *, GlobalLU_t *);
extern void    fixupL (const int_t, const int_t *, GlobalLU_t *);

extern void    sallocateA (int_t, int_t, float **, int_t **, int_t **);
extern void    sgstrf (superlu_options_t*, SuperMatrix*,
                       int_t, int_t, int_t*, void *, int_t, int_t *, int_t *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int_t *);
extern int_t     ssnode_dfs (const int_t, const int_t, const int_t *, const int_t *,
			     const int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     ssnode_bmod (const int_t, const int_t, const int_t, float *,
                              float *, GlobalLU_t *, SuperLUStat_t*);
extern void    spanel_dfs (const int_t, const int_t, const int_t, SuperMatrix *,
			   int_t *, int_t *, float *, int_t *, int_t *, int_t *,
			   int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern void    spanel_bmod (const int_t, const int_t, const int_t, const int_t,
                           float *, float *, int_t *, int_t *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int_t     scolumn_dfs (const int_t, const int_t, int_t *, int_t *, int_t *, int_t *,
			   int_t *, int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     scolumn_bmod (const int_t, const int_t, float *,
			   float *, int_t *, int_t *, int_t,
                           GlobalLU_t *, SuperLUStat_t*);
extern int_t     scopy_to_ucol (int_t, int_t, int_t *, int_t *, int_t *,
                              float *, GlobalLU_t *);         
extern int_t     spivotL (const int_t, const double, int_t *, int_t *, 
                         int_t *, int_t *, int_t *, GlobalLU_t *, SuperLUStat_t*);
extern void    spruneL (const int_t, const int_t *, const int_t, const int_t,
			  const int_t *, const int_t *, int_t *, GlobalLU_t *);
extern void    sreadmt (int *, int *, int *, float **, int **, int **);
extern void    sGenXtrue (int_t, int_t, float *, int_t);
extern void    sFillRHS (trans_t, int_t, float *, int_t, SuperMatrix *,
			  SuperMatrix *);
extern void    sgstrs (trans_t, SuperMatrix *, SuperMatrix *, int_t *, int_t *,
                        SuperMatrix *, SuperLUStat_t*, int_t *);
/* ILU */
extern void    sgsitrf (superlu_options_t*, SuperMatrix*, int_t, int_t, int_t*,
		        void *, int_t, int_t *, int_t *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int_t *);
extern int     sldperm(int_t, int_t, int_t, int_t [], int_t [], float [],
                        int_t [],	float [], float []);
extern int_t     ilu_ssnode_dfs (const int_t, const int_t, const int_t *, const int_t *,
			       const int_t *, int_t *, GlobalLU_t *);
extern void    ilu_spanel_dfs (const int_t, const int_t, const int_t, SuperMatrix *,
			       int_t *, int_t *, float *, float *, int_t *, int_t *,
			       int_t *, int_t *, int_t *, int_t *, GlobalLU_t *);
extern int_t     ilu_scolumn_dfs (const int_t, const int_t, int_t *, int_t *, int_t *,
				int_t *, int_t *, int_t *, int_t *, int_t *,
				GlobalLU_t *);
extern int_t     ilu_scopy_to_ucol (int_t, int_t, int_t *, int_t *, int_t *,
                                  float *, int_t, milu_t, double, int_t,
                                  float *, int_t *, GlobalLU_t *, float *);
extern int_t     ilu_spivotL (const int_t, const double, int_t *, int_t *, int_t, int_t *,
			    int_t *, int_t *, int_t *, double, milu_t,
                            float, GlobalLU_t *, SuperLUStat_t*);
extern int_t     ilu_sdrop_row (superlu_options_t *, int_t, int_t, double,
                              int_t, int_t *, double *, GlobalLU_t *, 
                              float *, float *, int_t);


/*! \brief Driver related */

extern void    sgsequ (SuperMatrix *, float *, float *, float *,
			float *, float *, int_t *);
extern void    slaqgs (SuperMatrix *, float *, float *, float,
                        float, float, char *);
extern void    sgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, int_t *);
extern float   sPivotGrowth(int_t, SuperMatrix *, int_t *, 
                            SuperMatrix *, SuperMatrix *);
extern void    sgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int_t *, int_t *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, int_t *);

extern int_t     sp_strsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, float *, SuperLUStat_t*, int_t *);
extern int     sp_sgemv (char *, float, SuperMatrix *, float *,
			int, float, float *, int);

extern int     sp_sgemm (char *, char *, int, int, int, float,
			SuperMatrix *, float *, int, float, 
			float *, int);
extern         float smach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int_t     sLUMemInit (fact_t, void *, int_t, int_t, int_t, int_t, int_t,
                            float, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int_t **, float **);
extern void    sSetRWork (int_t, int_t, float *, float **, float **);
extern void    sLUWorkFree (int_t *, float *, GlobalLU_t *);
extern int_t     sLUMemXpand (int_t, int_t, MemType, int_t *, GlobalLU_t *);

extern float  *floatMalloc(int_t);
extern float  *floatCalloc(int_t);
extern int_t     smemory_usage(const int_t, const int_t, const int_t, const int_t);
extern int_t     sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern int_t     ilu_sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    sreadhb(FILE *, int_t *, int_t *, int_t *, float **, int_t **, int_t **);
extern void    sreadrb(int_t *, int_t *, int_t *, float **, int_t **, int_t **);
extern void    sreadtriple(int_t *, int_t *, int_t *, float **, int_t **, int_t **);
extern void    sreadMM(FILE *, int_t *, int_t *, int_t *, float **, int_t **, int_t **);
extern void    sCompRow_to_CompCol(int_t, int_t, int_t, float*, int_t*, int_t*,
		                   float **, int_t **, int_t **);
extern void    sfill (float *, int_t, float);
extern void    sinf_norm_error (int_t, SuperMatrix *, float *);
extern float  sqselect(int, float *, int);


/*! \brief Routines for debugging */
extern void    sPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    sPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    sPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    sprint_t_lu_col(char *, int_t, int_t, int_t *, GlobalLU_t *);
extern int     print_double_vec(char *, int, double *);
extern void    scheck_tempv(int_t, float *);

/*! \brief BLAS */

extern int sgemm_(const char*, const char*, const int*, const int*, const int*,
                  const float*, const float*, const int*, const float*,
		  const int*, const float*, float*, const int*);
extern int strsv_(char*, char*, char*, int*, float*, int*,
                  float*, int*);
extern int strsm_(char*, char*, char*, char*, int*, int*,
                  float*, float*, int*, float*, int*);
extern int sgemv_(char *, int *, int *, float *, float *a, int *,
                  float *, int *, float *, float *, int *);

extern void susolve(int, int, float*, float*);
extern void slsolve(int, int, float*, float*);
extern void smatvec(int, int, int, float*, float*, float*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_sSP_DEFS */

