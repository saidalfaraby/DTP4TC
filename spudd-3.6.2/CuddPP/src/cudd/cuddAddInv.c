/**CFile***********************************************************************

  FileName    [cuddAddInv.c]

  PackageName [cudd]

  Synopsis    [Function to compute the scalar inverse of an ADD.]

  Description [External procedures included in this module:
		<ul>
		<li> Cudd_addScalarInverse()
		</ul>
	    Internal procedures included in this module:
		<ul>
		<li> cuddAddScalarInverseRecur()
		</ul>]

  Author      [Fabio Somenzi]

  Copyright   [This file was created at the University of Colorado at
  Boulder.  The University of Colorado at Boulder makes no warranty
  about the suitability of this software for any purpose.  It is
  presented on an AS IS basis.]

******************************************************************************/

#include "util.h"
#include "cuddInt.h"


/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Stucture declarations                                                     */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] DD_UNUSED = "$Id: cuddAddInv.c,v 2.0 2003/02/07 00:12:30 staubin Exp $";
#endif


/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/


/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/


/**AutomaticEnd***************************************************************/


/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/


/**Function********************************************************************

  Synopsis    [Computes the scalar inverse of an ADD.]
  
  Description [Computes an n ADD where the discriminants are the
  multiplicative inverses of the corresponding discriminants of the
  argument ADD.  Returns a pointer to the resulting ADD in case of
  success. Returns NULL if any discriminants smaller than epsilon is
  encountered.]

  SideEffects [None]

******************************************************************************/
DdNode *
Cudd_addScalarInverse(
  DdManager * dd,
  DdNode * f,
  DdNode * epsilon)
{
    DdNode *res;

    if (!cuddIsConstant(epsilon)) {
	(void) fprintf(dd->err,"Invalid epsilon\n");
	return(NULL);
    }
    do {
	dd->reordered = 0;
	res  = cuddAddScalarInverseRecur(dd,f,epsilon);
    } while (dd->reordered == 1);
    return(res);

} /* end of Cudd_addScalarInverse */

/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/


/**Function********************************************************************

  Synopsis    [Performs the recursive step of addScalarInverse.]

  Description [Returns a pointer to the resulting ADD in case of
  success. Returns NULL if any discriminants smaller than epsilon is
  encountered.]

  SideEffects [None]

******************************************************************************/
DdNode *
cuddAddScalarInverseRecur(
  DdManager * dd,
  DdNode * f,
  DdNode * epsilon)
{
    DdNode *t, *e, *res;
    ///xxxx
    //CUDD_VALUE_TYPE value;
    Terminal terminal;
    //(*terminal).set(0.0);
    statLine(dd);
    if (cuddIsConstant(f)) {
	if (ddAbs((*cuddV(f))) < (*cuddV(epsilon))) return(NULL);
	//value = 1.0 / cuddV(f);
	//value = cuddV(f) / 1.0;
	terminal = (*cuddV(f)) / 1.0;
	//res = cuddUniqueConst(dd,value);
	res = cuddUniqueConst(dd,&terminal);
	return(res);
    }

    res = cuddCacheLookup2(dd,Cudd_addScalarInverse,f,epsilon);
    if (res != NULL) return(res);

    t = cuddAddScalarInverseRecur(dd,cuddT(f),epsilon);
    if (t == NULL) return(NULL);
    cuddRef(t);

    e = cuddAddScalarInverseRecur(dd,cuddE(f),epsilon);
    if (e == NULL) {
	Cudd_RecursiveDeref(dd, t);
	return(NULL);
    }
    cuddRef(e);

    res = (t == e) ? t : cuddUniqueInter(dd,(int)f->index,t,e);
    if (res == NULL) {
	Cudd_RecursiveDeref(dd, t);
	Cudd_RecursiveDeref(dd, e);
	return(NULL);
    }

    cuddCacheInsert2(dd,Cudd_addScalarInverse,f,epsilon,res);

    return(res);

} /* end of cuddAddScalarInverseRecur */


/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/

