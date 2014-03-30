/**CHeaderFile*****************************************************************

  FileName    [dddmpInt.h]

  PackageName [dddmp]

  Synopsis    [Low level functions to read in and write out bdds to file]

  Description [A set of internal low-level routines of the dddmp package doing:
               <ul>
                 <li> read and write of node codes in binary mode,
                 <li> read and write of integers in binary mode,
		 <li> marking/unmarking nodes as visited,
                 <li> numbering nodes.
	       </ul>
              ]

  Author      [Gianpiero Cabodi & Stefano Quer]

  Copyright   [Politecnico di Torino(Italy) ]

******************************************************************************/

#ifndef _DDDMPINT
#define _DDDMPINT

#include	"dddmp.h"
#include	"cuddInt.h"


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

/**Enum************************************************************************

  Synopsis    [Type for supported decomposition types.]

  Description [Type for supported decomposition types. Used internally
  to select the proper type (bdd, add, ...). ]

******************************************************************************/
typedef enum {
    DDDMP_BDD,
    DDDMP_ADD
} Dddmp_DecompType;

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
  int nvars;
  int nsuppvars;
  char **varnames;
  int *ids;
  int *permids;
  int *auxids;
  int nroots;
  int *rootids;
  char **rootnames;
};	

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

#ifdef ALLOC
#define DDDMP_ALLOC(type, num)	ALLOC(type,num)
#else
#define DDDMP_ALLOC(type, num)	\
    ((type *) malloc(sizeof(type) * (num)))
#endif

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
EXTERN int DddmpNumberDdNodes(DdManager *dd, DdNode **f, int n);
EXTERN void DddmpUnnumberDdNodes(DdManager *dd, DdNode **f, int n);
EXTERN void DddmpWriteNodeIndex(DdNode *f, int id);
EXTERN int DddmpReadNodeIndex(DdNode *f);
EXTERN int DddmpVisited(DdNode *f);
EXTERN void DddmpSetVisited(DdNode *f);
EXTERN void DddmpClearVisited(DdNode *f);
EXTERN int NumberNodeRecur(DdNode *f, int id);
EXTERN int DddmpCuddDdArrayStore(Dddmp_DecompType dd_type, DdManager *dd, char *ddname, int nroots, DdNode **f, char **rootnames, char **varnames, int *auxids, int mode, Dddmp_VarInfoType varinfo, char *fname, FILE *fp);
EXTERN int DddmpCuddDdArrayLoad(Dddmp_DecompType dd_type, DdManager *dd, Dddmp_RootMatchType rootmatchmode, char **rootmatchnames, Dddmp_VarMatchType varmatchmode, char **varmatchnames, int *varmatchauxids, int *varcomposeids, int mode, char *file, FILE *fp, DdNode ***pproots);
EXTERN Dddmp_Hdr_t * DddmpBddReadHeader(char *file, FILE *fp);
EXTERN void DddmpFreeHeader(Dddmp_Hdr_t *Hdr);

/**AutomaticEnd***************************************************************/

#endif
