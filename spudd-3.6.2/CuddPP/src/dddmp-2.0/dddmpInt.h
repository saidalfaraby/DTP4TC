/**CHeaderFile*****************************************************************

  FileName     [dddmpInt.h]

  PackageName  [dddmp]

  Synopsis     [Low level functions to read in and write out bdds to file]

  Description  [A set of internal low-level routines of the dddmp package
    doing:
    <ul>
      <li> read and write of node codes in binary mode,
      <li> read and write of integers in binary mode,
      <li> marking/unmarking nodes as visited,
      <li> numbering nodes.
    </ul>
    ]

  Author       [Gianpiero Cabodi and Stefano Quer]

  Copyright    [This file was created at the Politecnico di Torino,
    Torino, Italy. 
    The  Politecnico di Torino makes no warranty about the suitability of 
    this software for any purpose.  
    It is presented on an AS IS basis.
    ]

******************************************************************************/

#ifndef _DDDMPINT
#define _DDDMPINT

#include "dddmp.h"
#include "cuddInt.h"

/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/

/* constants for code fields */
#define DDDMP_TERMINAL      0
#define DDDMP_ABSOLUTE_ID   1
#define DDDMP_RELATIVE_ID   2
#define DDDMP_RELATIVE_1    3

#define DDDMP_MAXSTRLEN 500

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Structure declarations                                                     */
/*---------------------------------------------------------------------------*/

/**Struct**********************************************************************
 Synopsis    [used in binary mode to store code info of a dd node]
 Description [V , T , E store the mode used to represent variable, Then 
              and Else indexes. An index is either an absolute
	      ( DDDMP_ABSOLUTE_ID ),
              a relative numbers ( DDDMP_RELATIVE_ID , DDDMP_RELATIVE_1 ) or 
              a terminal node ( DDDMP_TERMINAL ) .
	      Ecomp is used for the complemented edge attribute.
             ]
 SideEffect  [none]
 SeeAlso     [DddmpWriteCode DddmpReadCode] 
******************************************************************************/

struct binary_dd_code {
  unsigned  Unused : 1;
  unsigned  V      : 2;
  unsigned  T      : 2;
  unsigned  Ecompl : 1;
  unsigned  E      : 2;
};

/**Struct*********************************************************************

 Synopsis    [BDD file header]

 Description [Structure containing the BDD header file infos]

******************************************************************************/

struct Dddmp_Hdr_s {
  char *ver;
  char mode;
  Dddmp_DecompType dd_type;
  Dddmp_VarInfoType varinfo;
  char *dd;
  int nnodes;
  int nVars;
  int nsuppvars;
  char **orderedVarNames;
  char **suppVarNames;
  int *ids;
  int *permids;
  int *auxids;
  int *cnfids;
  int nRoots;
  int *rootids;
  char **rootnames;
  int nAddedCnfVar;
  int nVarsCnf;
  int nClausesCnf;  
};	

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

/**Macro***********************************************************************

  Synopsis     [Memory Allocation Macro for DDDMP]

  Description  []

  SideEffects  [None]

  SeeAlso      []

******************************************************************************/

#ifdef ALLOC
#  define DDDMP_ALLOC(type, num)	ALLOC(type,num)
#else
#  define DDDMP_ALLOC(type, num)	\
     ((type *) malloc(sizeof(type) * (num)))
#endif

/**Macro***********************************************************************

  Synopsis     [Memory Free Macro for DDDMP]

  Description  []

  SideEffects  [None]

  SeeAlso      []

******************************************************************************/

#ifdef FREE
#define DDDMP_FREE(p)  (FREE(p))
#else
#define DDDMP_FREE(p)	\
    ((p)!=NULL)?(free(p)):0)
#endif


/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */
/*---------------------------------------------------------------------------*/

EXTERN int DddmpWriteCode(FILE *fp, struct binary_dd_code code);
EXTERN int DddmpReadCode(FILE *fp, struct binary_dd_code *pcode);
EXTERN int DddmpWriteInt(FILE *fp, int id);
EXTERN int DddmpReadInt(FILE *fp, int *pid);
EXTERN int DddmpNumberAddNodes(DdManager *ddMgr, DdNode **f, int n);
EXTERN void DddmpUnnumberAddNodes(DdManager *ddMgr, DdNode **f, int n);
EXTERN void DddmpWriteNodeIndexAdd(DdNode *f, int id);
EXTERN int DddmpReadNodeIndexAdd(DdNode *f);
EXTERN int DddmpVisitedAdd(DdNode *f);
EXTERN void DddmpSetVisitedAdd(DdNode *f);
EXTERN void DddmpClearVisitedAdd(DdNode *f);
EXTERN int DddmpNumberBddNodes(DdManager *ddMgr, DdNode **f, int n);
EXTERN void DddmpUnnumberBddNodes(DdManager *ddMgr, DdNode **f, int n);
EXTERN void DddmpWriteNodeIndexBdd(DdNode *f, int id);
EXTERN int DddmpReadNodeIndexBdd(DdNode *f);
EXTERN int DddmpVisitedBdd(DdNode *f);
EXTERN void DddmpSetVisitedBdd(DdNode *f);
EXTERN void DddmpClearVisitedBdd(DdNode *f);
EXTERN int DddmpNumberDdNodesCnf(DdManager *ddMgr, DdNode **f, int rootN, int *cnfIds, int id);
EXTERN int DddmpDdNodesCountEdgesAndNumber(DdManager *ddMgr, DdNode **f, int rootN, int edgeInTh, int pathLengthTh, int *cnfIds, int id);
EXTERN void DddmpUnnumberDdNodesCnf(DdManager *ddMgr, DdNode **f, int rootN);
EXTERN int DddmpPrintBddAndNext(DdManager *ddMgr, DdNode **f, int rootN);
EXTERN int DddmpWriteNodeIndexCnf(DdNode *f, int id);
EXTERN int DddmpVisitedCnf(DdNode *f);
EXTERN void DddmpSetVisitedCnf(DdNode *f);
EXTERN int DddmpReadNodeIndexCnf(DdNode *f);
EXTERN int DddmpCuddDdArrayStoreBdd(Dddmp_DecompType dd_type, DdManager *ddMgr, char *ddname, int nRoots, DdNode **f, char **rootnames, char **varnames, int *auxids, int mode, Dddmp_VarInfoType varinfo, char *fname, FILE *fp);
EXTERN int DddmpCuddBddArrayStore(Dddmp_DecompType dd_type, DdManager *ddMgr, char *ddname, int nRoots, DdNode **f, char **rootnames, char **varnames, int *auxids, int mode, Dddmp_VarInfoType varinfo, char *fname, FILE *fp);
EXTERN int QsortStrcmp(const void *ps1, const void *ps2);
EXTERN int FindVarname(char *name, char **array, int n);
EXTERN char * DddmpStrDup(char *str);
EXTERN char ** DddmpStrArrayDup(char **array, int n);
EXTERN char ** DddmpStrArrayRead(FILE *fp, int n);
EXTERN int DddmpStrArrayWrite(FILE *fp, char **array, int n);
EXTERN void DddmpStrArrayFree(char **array, int n);
EXTERN int * DddmpIntArrayDup(int *array, int n);
EXTERN int * DddmpIntArrayRead(FILE *fp, int n);
EXTERN int DddmpIntArrayWrite(FILE *fp, int *array, int n);

/**AutomaticEnd***************************************************************/

#endif
