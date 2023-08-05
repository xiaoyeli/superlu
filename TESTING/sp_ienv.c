/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! \file
 * Choose machine-dependent parameters for the
 * local environment.
 *
 * History: Modified from lapack routine ILAENV
 *
 * \ingroup Testing
 */

#include "../SRC/slu_Cnames.h"

/*!
 * sp_ienv() is inquired to choose machine-dependent parameters for the
 * local environment.
 *
 * See ISPEC for a description of the parameters.
 *
 * This version provides a set of parameters which should give good,
 * but not optimal, performance on many of the currently available
 * computers.  Users are encouraged to modify this subroutine to set
 * the tuning parameters for their particular machine using the option
 * and problem size information in the arguments.
 *
 * \param[in] ispec   Specifies the parameter to be returned as the value of SP_IENV.<br/>
 *         = 1: the panel size w; a panel consists of w consecutive
 *         columns of matrix A in the process of Gaussian elimination.
 *         The best value depends on machine's cache characters.<br/>
 *         = 2: the relaxation parameter relax; if the number of
 *         nodes (columns) in a subtree of the elimination tree is less
 *         than relax, this subtree is considered as one supernode,
 *         regardless of the their row structures.<br/>
 *         = 3: the maximum size for a supernode;<br/>
 *         = 4: the minimum row dimension for 2-D blocking to be used;<br/>
 *         = 5: the minimum column dimension for 2-D blocking to be used;<br/>
 *         = 6: the estimated fills factor for L and U, compared with A;
 *
 * \result >= 0: the value of the parameter specified by ISPEC<br/>
 *         < 0:  if SP_IENV = -k, the k-th argument had an illegal value.
 */
int
sp_ienv(int ispec)
{
    switch (ispec) {
	case 1: return (3);
	case 2: return (2);
	case 3: return (10);
	case 4: return (20);
	case 5: return (10);
        case 6: return (2);
    }

    /* Invalid value for ISPEC */
    return (-1);

} /* sp_ienv_ */
