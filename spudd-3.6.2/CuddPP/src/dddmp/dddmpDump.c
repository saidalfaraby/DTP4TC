/**CFile**********************************************************************

  FileName    [dddmpDump.c]

  PackageName [dddmp]

  Synopsis    [Functions to read in and write out bdds to file]

  Description [Functions to read in and write out bdds to file.  BDDs
  are represended on file either in text or binary format under the
  following rules.  A file contains a forest of BDDs (a vector of
  Boolean functions).  BDD nodes are numbered with contiguous numbers,
  from 1 to NNodes (total number of nodes on a file). 0 is not used to
  allow negative node indexes for complemented edges.  A file contains
  a header, including information about variables and roots to BDD
  functions, followed by the list of nodes.  BDD nodes are listed
  according to their numbering, and in the present implementation
  numbering follows a post-order strategy, in such a way that a node
  is never listed before its Then/Else children.  ]

  Author      [Gianpiero Cabodi & Stefano Quer]

  Copyright   [Politecnico di Torino (Italy)]

******************************************************************************/
#include "dddmpInt.h"

/*---------------------------------------------------------------------------*/
/* Stucture declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

#define matchkeywd(str,key) (strncmp(str,key,strlen(key))==0)

/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

static int StoreNodeRecur(DdManager *dd, DdNode *f, int mode, int *supportids, char **varnames, int *outids, FILE *fp);
static int QsortStrcmp(const void *ps1, const void *ps2);
static int FindVarname(char *name, char **array, int n);
static char * DddmpStrDup(char *str);
static char ** DddmpStrArrayDup(char **array, int n);
static char ** DddmpStrArrayRead(FILE *fp, int n);
static int DddmpStrArrayWrite(FILE *fp, char **array, int n);
static void DddmpStrArrayFree(char **array, int n);
static int * DddmpIntArrayDup(int *array, int n);
static int * DddmpIntArrayRead(FILE *fp, int n);
static int DddmpIntArrayWrite(FILE *fp, int *array, int n);

/**AutomaticEnd***************************************************************/

/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis    [Writes a dump file representing the argument BDD.]

  Description [Dumps the argument BDD to file. Dumping is done through
               Dddmp_cuddBddArrayStore, And a dummy array of 1 BDD root is
               used for this purpose.
               ]

  SideEffects [Nodes are temporarily removed from unique hash. They are 
re-linked after the store operation in a modified order.]

  SeeAlso     [Dddmp_cuddBddLoad Dddmp_cuddBddArrayLoad]

******************************************************************************/
int
Dddmp_cuddBddStore (
  DdManager *dd          /* manager */,
  char *ddname           /* dd name (or NULL) */,
  DdNode *f              /* BDD root to be stored */,
  char **varnames        /* array of variable names (or NULL) */,
  int *auxids            /* array of converted var ids */,
  int mode               /* storing mode selector */,
  Dddmp_VarInfoType varinfo /* extra info for variables in text mode */,
  char *fname            /* file name */,
  FILE *fp               /* pointer to the store file */ 
)
{
  DdNode *FArray[1];

  FArray[0] = f;
  return (Dddmp_cuddBddArrayStore (dd,ddname,1,FArray,NULL,
          varnames,auxids,mode,varinfo,fname,fp));

}

/**Function********************************************************************

  Synopsis    [Writes a dump file representing the argument Array of BDDs.]

  Description [Dumps the argument array of BDDs to file. Dumping is
  either in text or binary form.  BDDs are stored to the fp (already
  open) file if not NULL. Otherwise the file whose name is fname is opened
  in write mode.  The header has the same format for both textual and
  binary dump.  Names are allowed for input variables (vnames) and for
  represented functions (rnames).  For sake of generality and because
  of dynamic variable ordering both variable IDs and permuted IDs are
  included. New IDs are also supported (auxids).  Variables are identified with incremental       
  numbers. according with their positiom in the support set.
  In text mode, an extra info may be added, chosen among the following options:
  name, ID, PermID, or an auxiliary id.  Since conversion from DD pointers to integers is
  required, DD nodes are temporarily removed from the unique
  hash table. This allows the use of the next field to store node IDs.  ]

  SideEffects [Nodes are temporarily removed from the unique hash
  table. They are re-linked after the store operation in a modified
  order.]

  SeeAlso     [Dddmp_cuddBddStore, Dddmp_cuddBddLoad, Dddmp_cuddBddArrayLoad]

******************************************************************************/
int
Dddmp_cuddBddArrayStore (
  DdManager *dd             /* manager */,
  char *ddname              /* dd name (or NULL) */,
  int nroots                /* number of output BDD roots to be stored */,
  DdNode **f                /* array of BDD roots to be stored */,
  char **rootnames          /* array of root names (or NULL) */,
  char **varnames           /* array of variable names (or NULL) */,
  int *auxids               /* array of converted var IDs */,
  int mode                  /* storing mode selector */,
  Dddmp_VarInfoType varinfo /* extra info for variables in text mode */,
  char *fname               /* file name */,
  FILE *fp                  /* pointer to the store file */ 
)
{

  return(DddmpCuddDdArrayStore(DDDMP_BDD,dd,ddname,nroots,f,
           rootnames,varnames,auxids,mode,varinfo,fname,fp));

}

/**Function********************************************************************

  Synopsis    [Writes a dump file representing the argument ADD.]

  Description [Dumps the argument ADD to file. Dumping is done through
               Dddmp_cuddAddArrayStore, And a dummy array of 1 ADD root is
               used for this purpose.
               ]

  SideEffects [Nodes are temporarily removed from unique hash. They are 
re-linked after the store operation in a modified order.]

  SeeAlso     [Dddmp_cuddAddLoad Dddmp_cuddAddArrayLoad]

******************************************************************************/
int
Dddmp_cuddAddStore (
  DdManager *dd          /* manager */,
  char *ddname           /* dd name (or NULL) */,
  DdNode *f              /* BDD root to be stored */,
  char **varnames        /* array of variable names (or NULL) */,
  int *auxids            /* array of converted var ids */,
  int mode               /* storing mode selector */,
  Dddmp_VarInfoType varinfo /* extra info for variables in text mode */,
  char *fname            /* file name */,
  FILE *fp               /* pointer to the store file */ 
)
{
  DdNode *FArray[1];

  FArray[0] = f;
  return (Dddmp_cuddAddArrayStore (dd,ddname,1,FArray,NULL,
          varnames,auxids,mode,varinfo,fname,fp));

}

/**Function********************************************************************

  Synopsis    [Writes a dump file representing the argument Array of ADDs.]

  Description [Dumps the argument array of ADDs to file. Dumping is
  either in text or binary form. see the corresponding BDD dump function
  for further details.

  SideEffects [Nodes are temporarily removed from the unique hash
  table. They are re-linked after the store operation in a modified
  order.]

  SeeAlso     [Dddmp_cuddAddStore, Dddmp_cuddAddLoad, Dddmp_cuddAddArrayLoad]

******************************************************************************/
int
Dddmp_cuddAddArrayStore (
  DdManager *dd             /* manager */,
  char *ddname              /* dd name (or NULL) */,
  int nroots                /* number of output BDD roots to be stored */,
  DdNode **f                /* array of BDD roots to be stored */,
  char **rootnames          /* array of root names (or NULL) */,
  char **varnames           /* array of variable names (or NULL) */,
  int *auxids               /* array of converted var IDs */,
  int mode                  /* storing mode selector */,
  Dddmp_VarInfoType varinfo /* extra info for variables in text mode */,
  char *fname               /* file name */,
  FILE *fp                  /* pointer to the store file */ 
)
{

  return(DddmpCuddDdArrayStore(DDDMP_ADD,dd,ddname,nroots,f,
           rootnames,varnames,auxids,mode,varinfo,fname,fp));

}


/**Function********************************************************************

  Synopsis    [Reads a dump file representing the argument BDD.]

  Description [Reads a dump file representing the argument BDD.
  Dddmp_cuddBddArrayLoad is used through a dummy array.  ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddBddStore, Dddmp_cuddBddArrayLoad]

******************************************************************************/
DdNode *
Dddmp_cuddBddLoad (
  DdManager *dd           /* manager */,
  Dddmp_VarMatchType varmatchmode /* storing mode selector */,
  char **varmatchnames    /* array of variable names (accessed by IDs) */,
  int  *varmatchauxids    /* array of variable auxids (accessed by IDs) */,
  int  *varcomposeids     /* array of new ids (accessed by ids) */,
  int mode      /* requested input file format (checked against file format)*/,
  char *file		  /* file name */,
  FILE *fp                /* file pointer */
)
{
  DdNode *f , **FArray;
  int i, nroots;

  nroots = Dddmp_cuddBddArrayLoad(dd,DDDMP_ROOT_MATCHLIST,NULL,
             varmatchmode,varmatchnames,varmatchauxids,varcomposeids,
             mode,file,fp,&FArray);

  if (nroots == 0)
    return (NULL);
  else {
    f = FArray[0];
    if (nroots > 1) {
      printf ("Warning: %d BDD roots found in file. Only first retrieved.\n",
              nroots);
      for (i=1; i<nroots; i++)
        Cudd_RecursiveDeref(dd,FArray[i]);
    } 
    DDDMP_FREE(FArray);
    return f;
  }

}

/**Function********************************************************************

  Synopsis    [Reads a dump file representing the argument BDDs.]

  Description [Reads a dump file representing the argument BDDs. The header is
  common to both text and binary mode. The node list is either 
  in text or binary format. A dynamic vector of DD pointers 
  is allocated to support conversion from DD indexes to pointers.
  Several criteria are supported for variable match between file
  and dd manager. Several changes/permutations/compositions are allowed
  for variables while loading DDs. Variable of the dd manager are allowed 
  to match with variables on file on ids, permids, varnames, 
  varauxids; also direct composition between ids and 
  composeids is supported. More in detail:
  <ol>
  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  allows the loading of a DD keeping variable IDs unchanged
  (regardless of the variable ordering of the reading manager); this
  is useful, for example, when swapping DDs to file and restoring them
  later from file, after possible variable reordering activations.
  
  <li> varmatchmode=DDDMP_VAR_MATCHPERMIDS <p>
  is used to allow variable match according to the position in the ordering.
  
  <li> varmatchmode=DDDMP_VAR_MATCHNAMES <p>
  requires a non NULL varmatchnames parameter; this is a vector of
  strings in one-to-one correspondence with variable IDs of the
  reading manager. Variables in the DD file read are matched with
  manager variables according to their name (a non NULL varnames
  parameter was required while storing the DD file).

  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  has a meaning similar to DDDMP_VAR_MATCHNAMES, but integer auxiliary
  IDs are used instead of strings; the additional non NULL
  varmatchauxids parameter is needed.

  <li> varmatchmode=DDDMP_VAR_COMPOSEIDS <p>
  uses the additional varcomposeids parameter is used as array of
  variable ids to be composed with ids stored in file.
  </ol>

  In the present implementation, the array varnames (3), varauxids (4)
  and composeids (5) need to have one entry for each variable in the 
  DD manager (NULL pointers are allowed for unused variables
  in varnames). Hence variables need to be already present in the 
  manager. All arrays are sorted according to IDs.
  ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddBddArrayStore]

******************************************************************************/
int
Dddmp_cuddBddArrayLoad (
  DdManager *dd           /* manager */,
  Dddmp_RootMatchType rootmatchmode /* storing mode selector */,
  char **rootmatchnames   /* sorted names for loaded roots */,
  Dddmp_VarMatchType varmatchmode /* storing mode selector */,
  char **varmatchnames    /* array of variable names (accessed by ids) */,
  int  *varmatchauxids    /* array of variable auxids (accessed by ids) */,
  int  *varcomposeids     /* array of new ids (accessed by ids) */,
  int mode                /* requested input file format (checked against file format)*/,
  char *file		  /* file name */,
  FILE *fp                /* file pointer */,
  DdNode ***pproots       /* array of returned BDD roots (by reference) */
)
{
   return (DddmpCuddDdArrayLoad(DDDMP_BDD,dd,rootmatchmode,rootmatchnames,
             varmatchmode,varmatchnames,varmatchauxids,varcomposeids,
             mode,file,fp,pproots));
}

/**Function********************************************************************

  Synopsis    [Reads a dump file representing the argument ADD.]

  Description [Reads a dump file representing the argument ADD.
  Dddmp_cuddAddArrayLoad is used through a dummy array.  ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddAddStore, Dddmp_cuddAddArrayLoad]

******************************************************************************/
DdNode *
Dddmp_cuddAddLoad (
  DdManager *dd           /* manager */,
  Dddmp_VarMatchType varmatchmode /* storing mode selector */,
  char **varmatchnames    /* array of variable names (accessed by IDs) */,
  int  *varmatchauxids    /* array of variable auxids (accessed by IDs) */,
  int  *varcomposeids     /* array of new ids (accessed by ids) */,
  int mode      /* requested input file format (checked against file format)*/,
  char *file		  /* file name */,
  FILE *fp                /* file pointer */
)
{
  DdNode *f , **FArray;
  int i, nroots;

  nroots = Dddmp_cuddAddArrayLoad(dd,DDDMP_ROOT_MATCHLIST,NULL,
             varmatchmode,varmatchnames,varmatchauxids,varcomposeids,
             mode,file,fp,&FArray);

  if (nroots == 0)
    return (NULL);
  else {
    f = FArray[0];
    if (nroots > 1) {
      printf ("Warning: %d BDD roots found in file. Only first retrieved.\n",
              nroots);
      for (i=1; i<nroots; i++)
        Cudd_RecursiveDeref(dd,FArray[i]);
    } 
    DDDMP_FREE(FArray);
    return f;
  }

}

/**Function********************************************************************

  Synopsis    [Reads a dump file representing the argument ADDs.]

  Description [Reads a dump file representing the argument ADDs. See 
  BDD load functions for detailed explanation. ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddBddArrayStore]

******************************************************************************/
int
Dddmp_cuddAddArrayLoad (
  DdManager *dd           /* manager */,
  Dddmp_RootMatchType rootmatchmode /* storing mode selector */,
  char **rootmatchnames   /* sorted names for loaded roots */,
  Dddmp_VarMatchType varmatchmode /* storing mode selector */,
  char **varmatchnames    /* array of variable names (accessed by ids) */,
  int  *varmatchauxids    /* array of variable auxids (accessed by ids) */,
  int  *varcomposeids     /* array of new ids (accessed by ids) */,
  int mode                /* requested input file format (checked against file format)*/,
  char *file		  /* file name */,
  FILE *fp                /* file pointer */,
  DdNode ***pproots       /* array of returned BDD roots (by reference) */
)
{
   return (DddmpCuddDdArrayLoad(DDDMP_ADD,dd,rootmatchmode,rootmatchnames,
             varmatchmode,varmatchnames,varmatchauxids,varcomposeids,
             mode,file,fp,pproots));
}

/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis    [Writes a dump file representing the argument Array of BDDs/ADDs.]

  Description [Dumps the argument array of BDDs/ADDs to file. Internal function
  doing inner steps of store for BDDs and ADDs. ADD store is presently supported 
  only with the text format. ]

  SideEffects [Nodes are temporarily removed from the unique hash
  table. They are re-linked after the store operation in a modified
  order.]

  SeeAlso     [Dddmp_cuddBddStore, Dddmp_cuddBddLoad, Dddmp_cuddBddArrayLoad]

******************************************************************************/
int
DddmpCuddDdArrayStore (
  Dddmp_DecompType dd_type  /* selects the proper decomp type: BDD or ADD */,
  DdManager *dd             /* manager */,
  char *ddname              /* dd name (or NULL) */,
  int nroots                /* number of output BDD roots to be stored */,
  DdNode **f                /* array of BDD roots to be stored */,
  char **rootnames          /* array of root names (or NULL) */,
  char **varnames           /* array of variable names (or NULL) */,
  int *auxids               /* array of converted var IDs */,
  int mode                  /* storing mode selector */,
  Dddmp_VarInfoType varinfo /* extra info for variables in text mode */,
  char *fname               /* file name */,
  FILE *fp                  /* pointer to the store file */ 
)
{
  DdNode      *support = NULL;
  DdNode      *scan;
  int         *ids = NULL;
  int         *permids = NULL;
  int         *invpermids = NULL;
  int         *supportids = NULL;
  int         *outids = NULL;
  char        **outvarnames = NULL;
  int         nvars = dd->size;
  int         nnodes;
  int         retval;
  int         i, var;
  int         close_fp = 0;


  if (fp == NULL) {
    /* 
     * File needs to be opened in the proper mode.
     */
    fp = fopen (fname, "w");
    if (fp == NULL) {
      (void) fprintf (stdout,"DdStore: Error opening %s\n",fname);
      goto failure;
    }
    close_fp = 1;
  }
  /* 
   * Force binary mode if automatic.
   */
  switch (mode) {
    case DDDMP_MODE_TEXT:
    case DDDMP_MODE_BINARY:
         break;
    case DDDMP_MODE_DEFAULT:
         mode = DDDMP_MODE_BINARY;
         break;
  }

  /* 
   * Alloc vectors for variable IDs, perm IDs and support IDs.
   *  +1 to include a slot for terminals.
   */
  ids = DDDMP_ALLOC(int,nvars);
  permids = DDDMP_ALLOC(int,nvars);
  invpermids = DDDMP_ALLOC(int,nvars);
  supportids = DDDMP_ALLOC(int,nvars+1);
  if ((ids == NULL)||(permids == NULL)||
      (invpermids == NULL)||(supportids == NULL)) {
    (void) fprintf (stdout,"DdStore: Error allocating memory\n");
    goto failure;
  }
    
  for (i = 0; i < nvars; i++) {
    ids[i] = permids[i] = invpermids[i] = supportids[i] = -1;
  }

  /* 
   * Take the union of the supports of each output function.
   * skip NULL functions.
   * Set permids and invpermids of support variables to the proper values.
   */
  for (i=0; i < nroots; i++) {
    if (f[i] == NULL) continue;
    support = Cudd_Support(dd,f[i]);
    if (support == NULL) {
      (void) fprintf (stdout,"DdStore Error: NULL support returned\n");
      goto failure;
    }
    cuddRef(support);
    scan = support;
    while (!cuddIsConstant(scan)) {
	ids[scan->index] = scan->index;
	permids[scan->index] = dd->perm[scan->index];
      invpermids[dd->perm[scan->index]] = scan->index;
	scan = cuddT(scan);
    }
    Cudd_RecursiveDeref(dd,support);
  }
  support = NULL; /* so that we do not try to free it in case of failure */

  /*
   * Set supportids to incremental (shrinked) values following the ordering.
   */
  for (i=0, var=0; i < nvars; i++) {
    if (invpermids[i] >= 0)
      supportids[invpermids[i]] = var++;
  }
  supportids[nvars] = var; /* set a dummy id for terminal nodes */

  /*
   * select conversion array for extra var info
   */
  switch (mode) {
    case DDDMP_MODE_TEXT:
      switch (varinfo) {
        case DDDMP_VARIDS:
          outids = ids;
          break;
        case DDDMP_VARPERMIDS:
          outids = permids;
          break;
        case DDDMP_VARAUXIDS:
          outids = auxids;
          break;
        case DDDMP_VARNAMES:
          outvarnames = varnames;
          break;
        case DDDMP_VARDEFAULT:
          break;
      }
      break;
    case DDDMP_MODE_BINARY:
      outids = NULL;
      break;
  }

  /* 
   * number dd nodes and count them (numbering is from 1 to nnodes)
   */
  nnodes = DddmpNumberDdNodes(dd,f,nroots);

  /* 
   * START HEADER 
   */

#ifdef DDDMP_VERSION
  retval = fprintf(fp,".ver %s\n", DDDMP_VERSION);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }
#endif


  if (dd_type == DDDMP_ADD) {
    retval = fprintf(fp,".add\n");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
  }

  retval = fprintf(fp,".mode %c\n", mode);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  if (mode == DDDMP_MODE_TEXT) {
    retval = fprintf(fp,".varinfo %d\n", varinfo);
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
  }

  if (ddname != NULL) {
    retval = fprintf(fp,".dd %s\n",ddname);
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
  }

  retval = fprintf(fp,".nnodes %d\n", nnodes);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  retval = fprintf(fp,".nvars %d\n", nvars);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  retval = fprintf(fp,".nsuppvars %d\n", var);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  if (varnames != NULL) {

    /* 
     * Write the var names by scanning the ids array.
     */

    retval = fprintf(fp,".varnames");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }

    for (i = 0; i < nvars; i++) {
	if (ids[i] >= 0) {
        if (varnames[ids[i]] == NULL) {
          (void) fprintf (stdout,"DdStore Warning: null variable name. DUMMY%d generated\n",i);
          varnames[ids[i]] = DDDMP_ALLOC(char,10);
          if (varnames[ids[i]] == NULL) {
            (void) fprintf (stdout,"DdStore: Error allocating memory\n");
            goto failure;
          }
          sprintf(varnames[ids[i]],"DUMMY%d",i);
        }
        retval = fprintf(fp," %s", varnames[ids[i]]);
        if (retval == EOF) {
          (void) fprintf (stdout,"DdStore: Error writing to file\n");
          goto failure;
        }
	}
    }

    retval = fprintf(fp,"\n");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }

  }

  /*
   * Write the var ids by scanning the ids array. 
   */
  retval = fprintf(fp,".ids");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }
  for (i = 0; i < nvars; i++) {
    if (ids[i] >= 0) {
	retval = fprintf(fp," %d", i);
      if (retval == EOF) {
        (void) fprintf (stdout,"DdStore: Error writing to file\n");
        goto failure;
      }
    }
  }
  retval = fprintf(fp,"\n");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  /*
   * Write the var permids by scanning the permids array. 
   */
  retval = fprintf(fp,".permids");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }
  for (i = 0; i < nvars; i++) {
    if (permids[i] >= 0) {
      retval = fprintf(fp," %d", permids[i]);
      if (retval == EOF) {
        (void) fprintf (stdout,"DdStore: Error writing to file\n");
        goto failure;
      }
    }
  }
  retval = fprintf(fp,"\n");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  if (auxids != NULL) {
  
    /*
     * Write the var auxids by scanning the ids array. 
     */
    retval = fprintf(fp,".auxids");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
    for (i = 0; i < nvars; i++) {
      if (ids[i] >= 0) {
        retval = fprintf(fp," %d", auxids[i]);
        if (retval == EOF) {
          (void) fprintf (stdout,"DdStore: Error writing to file\n");
          goto failure;
        }
      }
    }
    retval = fprintf(fp,"\n");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
  }

  /* 
   * Write the roots info. 
   */
  retval = fprintf(fp,".nroots %d\n", nroots);
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  if (rootnames != NULL) {

    /* 
     * Write the root names. 
     */

    retval = fprintf(fp,".rootnames");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
    for (i = 0; i < nroots; i++) {
      if (rootnames[i] == NULL) {
        (void) fprintf (stdout,"DdStore Warning: null variable name. ROOT%d generated\n",i);
        rootnames[i] = DDDMP_ALLOC(char,10);
        if (rootnames[i] == NULL) {
          (void) fprintf (stdout,"DdStore: Error allocating memory\n");
          goto failure;
        }
        sprintf(rootnames[ids[i]],"ROOT%d",i);
      }
      retval = fprintf(fp," %s", rootnames[i]);
      if (retval == EOF) {
        (void) fprintf (stdout,"DdStore: Error writing to file\n");
        goto failure;
      }
    }

    retval = fprintf(fp,"\n");
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }

  }

  retval = fprintf(fp,".rootids");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  /* 
   * Write BDD indexes of function roots.
   * Use negative integers for complemented edges. 
   */

  for (i = 0; i < nroots; i++) {
    if (f[i] == NULL) {
      (void) fprintf (stdout,"DdStore Warning: %d-th root is NULL\n",i);
      retval = fprintf(fp," 0");
    }
    if (Cudd_IsComplement(f[i])) {
      retval = fprintf(fp," -%d", DddmpReadNodeIndex(Cudd_Regular(f[i])));
    } 
    else {
      retval = fprintf(fp," %d", DddmpReadNodeIndex(Cudd_Regular(f[i])));
    }
    if (retval == EOF) {
      (void) fprintf (stdout,"DdStore: Error writing to file\n");
      goto failure;
    }
  }

  retval = fprintf(fp,"\n");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  retval = fprintf(fp,".nodes\n");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  /* 
   * END HEADER
   */

  /* 
   * Call the function that really gets the job done.
   */

  for (i = 0; i < nroots; i++) {
    if (f[i] != NULL) {
      retval = StoreNodeRecur(dd,Cudd_Regular(f[i]),
                              mode,supportids,outvarnames,outids,fp);
      if (retval == 0) {
        (void) fprintf (stdout,"DdStore: Error in recursive node store\n");
        goto failure;
      }
    }
  }

  /* 
   * Write trailer and return.
   */

  retval = fprintf(fp,".end\n");
  if (retval == EOF) {
    (void) fprintf (stdout,"DdStore: Error writing to file\n");
    goto failure;
  }

  if (close_fp)
    fclose (fp);

  DddmpUnnumberDdNodes(dd,f,nroots);
  DDDMP_FREE(ids);
  DDDMP_FREE(permids);
  DDDMP_FREE(invpermids);
  DDDMP_FREE(supportids);

  return(1);

  failure:

  if (ids != NULL) DDDMP_FREE(ids);
  if (permids != NULL) DDDMP_FREE(permids);
  if (invpermids != NULL) DDDMP_FREE(invpermids);
  if (supportids != NULL) DDDMP_FREE(supportids);
  if (support != NULL) Cudd_RecursiveDeref(dd,support);

  return(0);

}

/**Function********************************************************************

  Synopsis    [Reads a dump file representing the argument BDDs.]

  Description [Reads a dump file representing the argument BDDs. The header is
  common to both text and binary mode. The node list is either 
  in text or binary format. A dynamic vector of DD pointers 
  is allocated to support conversion from DD indexes to pointers.
  Several criteria are supported for variable match between file
  and dd manager. Several changes/permutations/compositions are allowed
  for variables while loading DDs. Variable of the dd manager are allowed 
  to match with variables on file on ids, permids, varnames, 
  varauxids; also direct composition between ids and 
  composeids is supported. More in detail:
  <ol>
  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  allows the loading of a DD keeping variable IDs unchanged
  (regardless of the variable ordering of the reading manager); this
  is useful, for example, when swapping DDs to file and restoring them
  later from file, after possible variable reordering activations.
  
  <li> varmatchmode=DDDMP_VAR_MATCHPERMIDS <p>
  is used to allow variable match according to the position in the ordering.
  
  <li> varmatchmode=DDDMP_VAR_MATCHNAMES <p>
  requires a non NULL varmatchnames parameter; this is a vector of
  strings in one-to-one correspondence with variable IDs of the
  reading manager. Variables in the DD file read are matched with
  manager variables according to their name (a non NULL varnames
  parameter was required while storing the DD file).

  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  has a meaning similar to DDDMP_VAR_MATCHNAMES, but integer auxiliary
  IDs are used instead of strings; the additional non NULL
  varmatchauxids parameter is needed.

  <li> varmatchmode=DDDMP_VAR_COMPOSEIDS <p>
  uses the additional varcomposeids parameter is used as array of
  variable ids to be composed with ids stored in file.
  </ol>

  In the present implementation, the array varnames (3), varauxids (4)
  and composeids (5) need to have one entry for each variable in the 
  DD manager (NULL pointers are allowed for unused variables
  in varnames). Hence variables need to be already present in the 
  manager. All arrays are sorted according to IDs.
  ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddBddArrayStore]

******************************************************************************/
int
DddmpCuddDdArrayLoad (
  Dddmp_DecompType dd_type /* selects the proper decomp type: BDD or ADD */,
  DdManager *dd           /* manager */,
  Dddmp_RootMatchType rootmatchmode /* storing mode selector */,
  char **rootmatchnames   /* sorted names for loaded roots */,
  Dddmp_VarMatchType varmatchmode /* storing mode selector */,
  char **varmatchnames    /* array of variable names (accessed by ids) */,
  int  *varmatchauxids    /* array of variable auxids (accessed by ids) */,
  int  *varcomposeids     /* array of new ids (accessed by ids) */,
  int mode                /* requested input file format (checked against file format)*/,
  char *file		  /* file name */,
  FILE *fp                /* file pointer */,
  DdNode ***pproots       /* array of returned BDD roots (by reference) */
)
{
Dddmp_Hdr_t *Hdr;
DdNode *f, *T, *E;
struct binary_dd_code code;
char buf[DDDMP_MAXSTRLEN];
int id, size, maxv;
int i, j, k, maxaux, 
    var, vT, vE, idT, idE;
double addConstant;
int  *permsupport = NULL;
int  *convertids = NULL;
int  *invconvertids = NULL;
int  *invauxids = NULL;
char **sortedvarnames = NULL;
int  nddvars, nroots;
DdNode **pnodes = NULL;
unsigned char *pvars1byte = NULL;
unsigned short *pvars2byte = NULL;
DdNode **proots = NULL;       /* array of BDD roots to be loaded */
int close_fp = 0;

  *pproots = NULL;

  if (fp == NULL) {
    fp = fopen (file, "r");
    if (fp == NULL) {
      (void) fprintf (stdout,"DdLoad: Error opening %s\n",file);
      goto failure;
    }
    close_fp = 1;
  }

  nddvars = dd->size;

  Hdr = DddmpBddReadHeader(NULL,fp);

  nroots = Hdr->nroots;

  if (Hdr->dd_type != dd_type) {
    (void) fprintf (stdout,"DdLoad Error: dd_type mismatch\n");
    if (Hdr->dd_type == DDDMP_BDD)
      (void) fprintf (stdout,"BDD found\n");
    if (Hdr->dd_type == DDDMP_ADD)
      (void) fprintf (stdout,"ADD found\n");
    if (dd_type == DDDMP_BDD)
      (void) fprintf (stdout,"when loading a BDD\n");
    if (dd_type == DDDMP_ADD)
      (void) fprintf (stdout,"when loading an ADD\n");
    goto failure;
  }

  if (Hdr->mode != mode) {
    if (mode == DDDMP_MODE_DEFAULT) {
      mode = Hdr->mode;
    }
    else {
      (void) fprintf (stdout,"DdLoad Error: mode mismatch\n");
      goto failure;
    }
  }

  /*
   * for each variable in the support, the relative position in the ordering
   * (within the support only) is computed
   */

  permsupport = DDDMP_ALLOC(int,Hdr->nsuppvars);
  if (permsupport == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }
  for (i=0,k=0; i < Hdr->nvars; i++) { 
    for (j=0; j < Hdr->nsuppvars; j++) { 
      if (Hdr->permids[j] == i) {
        permsupport[j] = k++;
      }
    }
  }
  assert (k==Hdr->nsuppvars);

  if (Hdr->varnames != NULL) {
    /*
     *  Varnames are sorted for binary search
     */
    sortedvarnames = DDDMP_ALLOC(char *,Hdr->nsuppvars);
    if (sortedvarnames == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
    for (i=0; i < Hdr->nsuppvars; i++) {
      if (Hdr->varnames[i] == NULL) {
        (void) fprintf (stdout,"DdLoad Error: support variable name missing in file\n");
        goto failure;
      } 
      sortedvarnames[i] = Hdr->varnames[i];
    }    
    
    qsort((void *)sortedvarnames,Hdr->nsuppvars,sizeof(char *),QsortStrcmp);
    
  }

  /*
   * convertids is the array used to convert variable ids from positional (shrinked)
   * ids used within the DD file. Positions in the file are from 0 to nsuppvars-1.
   */ 

  convertids = DDDMP_ALLOC(int,Hdr->nsuppvars);
  if (convertids == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }

  again_matchmode:
  switch (varmatchmode) {
    case DDDMP_VAR_MATCHIDS:
      for (i=0; i<Hdr->nsuppvars; i++)
        convertids[permsupport[i]] = Hdr->ids[i];
      break;
    case DDDMP_VAR_MATCHPERMIDS:
      for (i=0; i<Hdr->nsuppvars; i++)
        convertids[permsupport[i]] = Cudd_ReadInvPerm(dd,Hdr->permids[i]);
      break;
    case DDDMP_VAR_MATCHAUXIDS:
      if (Hdr->auxids == NULL) {
        (void) fprintf (stdout,"DdLoad Error: variable auxids matching requested\n");
        (void) fprintf (stdout,"but .auxids not found in BDD file\n");
        (void) fprintf (stdout,"Matching IDs forced.\n");
        varmatchmode = DDDMP_VAR_MATCHIDS;
        goto again_matchmode;
      }
      /* find max auxid value to alloc invaux array */
      for (i=0,maxaux= -1; i<nddvars; i++)
        if (varmatchauxids[i]>maxaux)
          maxaux = varmatchauxids[i];
      /* generate invaux array */
      invauxids = DDDMP_ALLOC(int,maxaux+1);
      if (invauxids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i<=maxaux; i++)
        invauxids[i] = -1;
      for (i=0; i<Hdr->nsuppvars; i++)
        invauxids[varmatchauxids[Hdr->ids[i]]] = Hdr->ids[i];
      /* generate convertids array */
      for (i=0; i<Hdr->nsuppvars; i++) {
        if ((Hdr->auxids[i]>maxaux) || (invauxids[Hdr->auxids[i]]<0)) {
          (void) fprintf (stdout,
                   "DdLoad Error: auxid %d not found in DD manager. ID matching forced (%d)\n", 
                   Hdr->auxids[i], i);
          (void) fprintf (stdout,"Beware of possible overlappings with other variables\n"); 
          convertids[permsupport[i]]=i;
        }
        else
          convertids[permsupport[i]] = invauxids[Hdr->auxids[i]];
      }
      break;
    case DDDMP_VAR_MATCHNAMES:
      if (Hdr->varnames == NULL) {
        (void) fprintf (stdout,"DdLoad Error: variable names matching requested\n");
        (void) fprintf (stdout,"but .varnames not found in BDD file\n");
        (void) fprintf (stdout,"Matching IDs forced.\n");
        varmatchmode = DDDMP_VAR_MATCHIDS;
        goto again_matchmode;
      }
      /* generate invaux array */
      invauxids = DDDMP_ALLOC(int,Hdr->nsuppvars);
      if (invauxids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i<Hdr->nsuppvars; i++)
        invauxids[i] = -1;
      for (i=0; i<nddvars; i++) {
        if (varmatchnames[i]==NULL) {
          (void) fprintf (stdout,"DdLoad Warning: NULL match variable name (id: %d). Ignored.\n",
                                 i);
        }
        else
          if ((j=FindVarname(varmatchnames[i],sortedvarnames,Hdr->nsuppvars))>=0) {
            assert(j<Hdr->nsuppvars);
            invauxids[j] = i;
          }
      }
      /* generate convertids array */
      for (i=0; i<Hdr->nsuppvars; i++) {
        assert (Hdr->varnames[i] != NULL);
        j=FindVarname(Hdr->varnames[i],sortedvarnames,Hdr->nsuppvars);
        assert((j>=0)&&(j<Hdr->nsuppvars));
        if (invauxids[j]<0) {
          (void) fprintf (stdout,
              "DdLoad Error: varname %s not found in DD manager. ID matching forced (%d)\n", 
               Hdr->varnames[i],i);
          convertids[permsupport[i]]=i;
        }
        else
          convertids[permsupport[i]] = invauxids[j];
      }
      break;
    case DDDMP_VAR_COMPOSEIDS:
      for (i=0; i<Hdr->nsuppvars; i++)
        convertids[permsupport[i]] = varcomposeids[Hdr->ids[i]];
      break;
  }

  maxv= -1;
  for (i=0; i<Hdr->nsuppvars; i++)
    if (convertids[i] > maxv)
      maxv = convertids[i];
 
  invconvertids = DDDMP_ALLOC(int,maxv+1);
  if (invconvertids == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }
  for (i=0; i<=maxv; i++)
    invconvertids[i]= -1;
  for (i=0; i<Hdr->nsuppvars; i++)
    invconvertids[convertids[i]] = i;

  pnodes = DDDMP_ALLOC(DdNode *,(Hdr->nnodes+1));
  if (pnodes == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }

  if (Hdr->nsuppvars < 256) {
    pvars1byte = DDDMP_ALLOC(unsigned char,(Hdr->nnodes+1));
    if (pvars1byte == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
  }
  else if (Hdr->nsuppvars < 0xffff) {
    pvars2byte = DDDMP_ALLOC(unsigned short,(Hdr->nnodes+1));
    if (pvars2byte == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
  }
  else {
    printf ("DdLoad Error: more than %d variables. Not supported.\n", 0xffff);
    goto failure;
  }

  for (i=1; i<=Hdr->nnodes; i++) {

    if (feof(fp)) {
      (void) fprintf (stdout,"DdLoad Error: unexpected EOF while reading DD nodes\n");
      goto failure;
    }

    switch (mode) {

      case DDDMP_MODE_TEXT:

        switch (Hdr->varinfo) {
          case DDDMP_VARIDS:
          case DDDMP_VARPERMIDS:
          case DDDMP_VARAUXIDS:
          case DDDMP_VARNAMES:
            if (fscanf(fp,"%d %*s %s %d %d\n", 
                           &id, buf, &idT, &idE) < 4) {
              (void) fprintf (stdout,"DdLoad: Error reading nodes in text mode\n");
              goto failure;
            }
            break;
          case DDDMP_VARDEFAULT:
            if (fscanf(fp,"%d %s %d %d\n", 
                           &id, buf, &idT, &idE) < 4) {
              (void) fprintf (stdout,"DdLoad: Error reading nodes in text mode\n");
              goto failure;
            }
            break;
        }
#ifdef DDDMP_DEBUG
        assert (id == i);
#endif
        if (idT==0 && idE==0)
	  {
	    /* leaf node: a constant */
	    if (strcmp(buf,"1")==0) {
	      pnodes[i] = Cudd_ReadOne(dd);       
	    }
	    else {
	      /* this is an ADD constant ! */
	      if (strcmp(buf,"0")==0) {
		pnodes[i] = Cudd_ReadZero(dd);       
	      }
	      else {
		addConstant = atof(buf);
		Terminal ADDconst;
		ADDconst.set(addConstant);
		//pnodes[i] = Cudd_addConst(dd, (CUDD_VALUE_TYPE) addConstant);       
		pnodes[i] = Cudd_addConst(dd,&ADDconst);       
	      }
	    }
	    if (pnodes[i]==NULL)
	      goto failure;
	    continue;
	  }else{
#ifdef DDDMP_DEBUG
          assert (idT > 0);
#endif
          var = atoi(buf);
          T = pnodes[idT];
          if(idE<0) {
            idE = -idE;
            E = pnodes[idE];
            E = Cudd_Not(E);
          }
          else
            E = pnodes[idE];
        }

      break;

      case DDDMP_MODE_BINARY:

        if (DddmpReadCode(fp,&code) == 0)
          goto failure;

        switch (code.V) {
        case DDDMP_TERMINAL:     
          /* only 1 terminal presently supported */    
          pnodes[i] = Cudd_ReadOne(dd);       
          continue; 
          break;
        case DDDMP_RELATIVE_1:
          break;
        case DDDMP_RELATIVE_ID:
        case DDDMP_ABSOLUTE_ID:
          size = DddmpReadInt(fp,&var);
          if (size == 0)
            goto failure;
          break;
        }
        switch (code.T) {
        case DDDMP_TERMINAL:     
          idT = 1;
          break;
        case DDDMP_RELATIVE_1:
          idT = i-1;
          break;
        case DDDMP_RELATIVE_ID:
          size = DddmpReadInt(fp,&id);
          if (size == 0)  goto failure;
          idT = i-id;
          break;
        case DDDMP_ABSOLUTE_ID:
          size = DddmpReadInt(fp,&idT);
          if (size == 0)  goto failure;
          break;
        }
        switch (code.E) {
        case DDDMP_TERMINAL:     
          idE = 1;
          break;
        case DDDMP_RELATIVE_1:
          idE = i-1;
          break;
        case DDDMP_RELATIVE_ID:
          size = DddmpReadInt(fp,&id);
          if (size == 0)  goto failure;
          idE = i-id;
          break;
        case DDDMP_ABSOLUTE_ID:
          size = DddmpReadInt(fp,&idE);
          if (size == 0)  goto failure;
          break;
        }

#ifdef DDDMP_DEBUG
      assert(idT<i);
#endif
      T = pnodes[idT];
      if (cuddIsConstant(T))
        vT = Hdr->nsuppvars;
      else {
        if (pvars1byte != NULL)
          vT = pvars1byte[idT];
        else if (pvars2byte != NULL)
          vT = pvars2byte[idT];
        else
          vT = invconvertids[T->index];
      }
#ifdef DDDMP_DEBUG
      assert (vT>0);
      assert (vT<=Hdr->nsuppvars);
#endif

#ifdef DDDMP_DEBUG
      assert(idE<i);
#endif
      E = pnodes[idE];
      if (cuddIsConstant(E))
        vE = Hdr->nsuppvars;
      else {
        if (pvars1byte != NULL)
          vE = pvars1byte[idE];
        else if (pvars2byte != NULL)
          vE = pvars2byte[idE];
        else
          vE = invconvertids[E->index];
      }
#ifdef DDDMP_DEBUG
      assert (vE>0);
      assert (vE<=Hdr->nsuppvars);
#endif
  
      switch (code.V) {
        case DDDMP_TERMINAL:     
        case DDDMP_ABSOLUTE_ID:
          break;
        case DDDMP_RELATIVE_1:
          var = (vT<vE) ? vT-1 : vE-1;
          break;
        case DDDMP_RELATIVE_ID:
          var = (vT<vE) ? vT-var : vE-var;
          break;
      }

      if (code.Ecompl)
        E = Cudd_Not(E);

#ifdef DDDMP_DEBUG
      assert (var<Hdr->nsuppvars);
#endif

      break;

    }

    if (pvars1byte != NULL)
      pvars1byte[i] = (unsigned char) var;
    else if (pvars2byte != NULL)
      pvars2byte[i] = (unsigned short) var;

    var = convertids[var];
    switch (dd_type) {
      case DDDMP_BDD: 
        pnodes[i] = Cudd_bddIte(dd, Cudd_bddIthVar(dd,var), T, E);
        break;
      case DDDMP_ADD: 
      { 
        DdNode *tmp = Cudd_addIthVar(dd,var);
        Cudd_Ref(tmp);
        pnodes[i] = Cudd_addIte(dd, tmp, T, E);
        Cudd_RecursiveDeref(dd,tmp);
        break;
      }
    }
    cuddRef(pnodes[i]);

  }

  fgets(buf, DDDMP_MAXSTRLEN-1,fp);
  if (!matchkeywd(buf, ".end")) {
    (void) fprintf (stdout,"DdLoad Error: .end not found\n");
    goto failure;
  }

  if (close_fp)
    fclose (fp);

  /* BDD Roots */
  proots = DDDMP_ALLOC(DdNode *,nroots);
  if (proots == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }

  for(i=0; i<nroots; ++i) {
    switch (rootmatchmode) {
      case DDDMP_ROOT_MATCHNAMES:
        for (j=0; j<nroots; j++) {
          if (strcmp(rootmatchnames[i],Hdr->rootnames[j])==0)
            break;
        }
        if (j>=nroots) { /* rootname not found */
          printf ("Warning: unable to match root name <%s>\n",
                  rootmatchnames[i]);
        }
        break; 
      case DDDMP_ROOT_MATCHLIST:
        j = i;
        break;
    }
    id = Hdr->rootids[i];
    if (id==0) {
      (void) fprintf (stdout,"DdLoad Warning: NULL root found in file\n");
      f = NULL;
    }
    else if (id<0) 
      f = Cudd_Not(pnodes[-id]);
    else
      f = pnodes[id];
    proots[i] = f;
    cuddRef(f);
  }

  for(i=2; i<=Hdr->nnodes; ++i) { 
    f = pnodes[i];
    Cudd_RecursiveDeref(dd, f);
  }

  /*
   * now free everithing was allocated within this function
   */

load_end:

  DddmpFreeHeader(Hdr);

  DDDMP_FREE(pnodes);
  DDDMP_FREE(pvars1byte);
  DDDMP_FREE(pvars2byte);

   /* variable names are not freed because they were shared with varnames */
  DDDMP_FREE(sortedvarnames);

  DDDMP_FREE(permsupport);
  DDDMP_FREE(convertids);
  DDDMP_FREE(invconvertids);
  DDDMP_FREE(invauxids);

  *pproots = proots;
  return nroots;

failure:

  if (close_fp)
    fclose (fp);

  nroots = 0; /* return 0 on error ! */

  DDDMP_FREE(proots);

  goto load_end; /* this is done to free memory */

}

/**Function********************************************************************

  Synopsis    [Reads a the header of a dump file representing the argument 
               BDDs.]

  Description [Reads the header of a dump file. Builds a Dddmp_Hdr_t struct
  containing all infos in the header, for next manipulations.
  ]

  SideEffects [none]

  SeeAlso     []

******************************************************************************/
Dddmp_Hdr_t *
DddmpBddReadHeader (
  char *file		  /* file name */,
  FILE *fp                /* file pointer */
)
{
Dddmp_Hdr_t *Hdr;
char buf[DDDMP_MAXSTRLEN];
int i;
int close_fp = 0;

  if (fp == NULL) {
    fp = fopen (file, "r");
    if (fp == NULL) {
      (void) fprintf (stdout,"DdReadHeader: Error opening %s\n",file);
      goto failure;
    }
    close_fp = 1;
  }

  /* START HEADER */

  Hdr = DDDMP_ALLOC(Dddmp_Hdr_t,1);
  if (Hdr == NULL) return NULL;

  Hdr->ver = NULL;
  Hdr->mode = 0;
  Hdr->dd_type = DDDMP_BDD;
  Hdr->varinfo = DDDMP_VARIDS;
  Hdr->dd = NULL;
  Hdr->nnodes = 0;
  Hdr->nvars = 0;
  Hdr->nsuppvars = 0;
  Hdr->varnames = NULL;
  Hdr->ids = NULL;
  Hdr->permids = NULL;
  Hdr->auxids = NULL;
  Hdr->nroots = 0;
  Hdr->rootids = NULL;
  Hdr->rootnames = NULL;

  while (fscanf(fp,"%s",buf)!=EOF) {

    /* comment */
    if (buf[0] == '#') {
      fgets(buf,DDDMP_MAXSTRLEN,fp);
      continue;
    }

    if (buf[0] != '.') {
      (void) fprintf (stdout,"DdReadHeader Error at\n%s\n", buf);
      (void) fprintf (stdout,"line must begin with '.' or '#'\n");
      goto failure;
    }

    if matchkeywd(buf, ".ver") {    
      /* this not checked so far: only read */
      if (fscanf (fp, "%s",buf)==EOF) {
        (void) fprintf (stdout,
          "DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      if ((Hdr->ver=DddmpStrDup(buf))==NULL) {
        (void) fprintf (stdout,"DdReadHeader: Error allocating memory\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".add") {    
      Hdr->dd_type = DDDMP_ADD;
      continue;
    }
    if matchkeywd(buf, ".bdd") {    
      Hdr->dd_type = DDDMP_BDD;
      continue;
    }

    if matchkeywd(buf, ".mode") {    
      if (fscanf (fp, "%s", buf)==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      Hdr->mode = buf[0];
      continue;
    }
    if matchkeywd(buf, ".varinfo") {    
      if (fscanf (fp, "%d", &(Hdr->varinfo))==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".dd") {    
      if (fscanf (fp, "%s", buf)==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      if ((Hdr->dd=DddmpStrDup(buf))==NULL) {
        (void) fprintf (stdout,"DdReadHeader: Error allocating memory\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".nnodes") {
      if (fscanf (fp, "%d", &(Hdr->nnodes))==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".nvars") {   
      if (fscanf (fp, "%d", &(Hdr->nvars))==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".nsuppvars") {
      if (fscanf (fp, "%d", &(Hdr->nsuppvars))==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".varnames") {
      Hdr->varnames = DddmpStrArrayRead(fp,Hdr->nsuppvars);
      if (Hdr->varnames == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".ids") {
      Hdr->ids = DddmpIntArrayRead(fp,Hdr->nsuppvars);
      if (Hdr->ids == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".permids") {
      Hdr->permids = DddmpIntArrayRead(fp,Hdr->nsuppvars);
      if (Hdr->permids == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".auxids") {
      Hdr->auxids = DddmpIntArrayRead(fp,Hdr->nsuppvars);
      if (Hdr->auxids == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".nroots") {
      if (fscanf (fp, "%d", &(Hdr->nroots))==EOF) {
        (void) fprintf (stdout,"DdReadHeader: Error reading file\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".rootids") {
      Hdr->rootids = DddmpIntArrayRead(fp,Hdr->nroots);
      if (Hdr->rootids == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".rootnames") {
      Hdr->rootnames = DddmpStrArrayRead(fp,Hdr->nroots);
      if (Hdr->rootnames == NULL) {
        goto failure;
      }
      else
        continue;
    }

    if matchkeywd(buf, ".nodes") {
      fgets(buf,DDDMP_MAXSTRLEN,fp);
      break;
    }

  }

  /* END HEADER */

  return Hdr;

failure:

  if (close_fp)
    fclose (fp);

  DddmpFreeHeader(Hdr);

  return NULL;

}


/**Function********************************************************************

  Synopsis    [Reads a the header of a dump file representing the argument 
               BDDs.]

  Description [Reads a dump file representing the argument BDDs. The header is
  common to both text and binary mode. The node list is either 
  in text or binary format. A dynamic vector of DD pointers 
  is allocated to support conversion from DD indexes to pointers.
  Several criteria are supported for variable match between file
  and dd manager. Several changes/permutations/compositions are allowed
  for variables while loading DDs. Variable of the dd manager are allowed 
  to match with variables on file on ids, permids, varnames, 
  varauxids; also direct composition between ids and 
  composeids is supported. More in detail:
  <ol>
  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  allows the loading of a DD keeping variable IDs unchanged
  (regardless of the variable ordering of the reading manager); this
  is useful, for example, when swapping DDs to file and restoring them
  later from file, after possible variable reordering activations.
  
  <li> varmatchmode=DDDMP_VAR_MATCHPERMIDS <p>
  is used to allow variable match according to the position in the ordering.
  
  <li> varmatchmode=DDDMP_VAR_MATCHNAMES <p>
  requires a non NULL varmatchnames parameter; this is a vector of
  strings in one-to-one correspondence with variable IDs of the
  reading manager. Variables in the DD file read are matched with
  manager variables according to their name (a non NULL varnames
  parameter was required while storing the DD file).

  <li> varmatchmode=DDDMP_VAR_MATCHIDS <p>
  has a meaning similar to DDDMP_VAR_MATCHNAMES, but integer auxiliary
  IDs are used instead of strings; the additional non NULL
  varmatchauxids parameter is needed.

  <li> varmatchmode=DDDMP_VAR_COMPOSEIDS <p>
  uses the additional varcomposeids parameter is used as array of
  variable ids to be composed with ids stored in file.
  </ol>

  In the present implementation, the array varnames (3), varauxids (4)
  and composeids (5) need to have one entry for each variable in the 
  DD manager (NULL pointers are allowed for unused variables
  in varnames). Hence variables need to be already present in the 
  manager. All arrays are sorted according to IDs.
  ]

  SideEffects [A vector of pointers to DD nodes is allocated and freed.]

  SeeAlso     [Dddmp_cuddBddArrayStore]

******************************************************************************/
void
DddmpFreeHeader (
  Dddmp_Hdr_t *Hdr   /* pointer to header */
)
{

  DDDMP_FREE(Hdr->ver);
  DDDMP_FREE(Hdr->dd);
  DddmpStrArrayFree(Hdr->varnames,Hdr->nsuppvars);
  DDDMP_FREE(Hdr->ids);
  DDDMP_FREE(Hdr->permids);
  DDDMP_FREE(Hdr->auxids);
  DDDMP_FREE(Hdr->rootids);
  DddmpStrArrayFree(Hdr->rootnames,Hdr->nroots);

  DDDMP_FREE(Hdr);

}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/


/**Function********************************************************************

  Synopsis    [Performs the recursive step of Dddmp_bddStore.]

  Description [Stores a node to file in either test or 
  binary mode.<l>
  In text mode a node is represented (on a text line basis) as
  <UL>
  <LI> node-index \[var-extrainfo\] var-index Then-index Else-index
  </UL>
  
  where all indexes are integer numbers and var-extrainfo (optional
  redundant field) is either an integer or a string (variable name).
  Node-index is redundant (due to the node ordering) but we keep it
  for readability.<p>

  In binary mode nodes are represented as a sequence of bytes,
  representing var-index, Then-index, and Else-index in an optimized
  way.  Only the first byte (code) is mandatory. Integer indexes are
  represented in absolute or relative mode, where relative means
  offset wrt. a Then/Else node info.  Suppose Var(NodeId),
  Then(NodeId) and Else(NodeId) represent infos about a given node.<p>

  The generic "NodeId" node is stored as 

  <UL>
  <LI> code-byte
  <LI> \[var-info\]
  <LI> \[Then-info\]
  <LI> \[Else-info\]
  </UL>

  where code-byte contains bit fields

  <UL>
  <LI>Unused  : 1 bit
  <LI>Variable: 2 bits, one of the following codes
    <UL>
    <LI>DDDMP_ABSOLUTE_ID   var-info = Var(NodeId) follows
    <LI>DDDMP_RELATIVE_ID   Var(NodeId) is represented in relative form as
    var-info = Min(Var(Then(NodeId)),Var(Else(NodeId))) -Var(NodeId)
    <LI>DDDMP_RELATIVE_1    No var-info follows, because
    Var(NodeId) = Min(Var(Then(NodeId)),Var(Else(NodeId)))-1
    <LI>DDDMP_TERMINAL      Node is a terminal, no var info required
    </UL>
  <LI>T       : 2 bits, with codes similar to V
    <UL>
    <LI>DDDMP_ABSOLUTE_ID   Then-info = Then(NodeId) follows
    <LI>DDDMP_RELATIVE_ID   Then(NodeId) is represented in relative form as
    Then-info = Nodeid-Then(NodeId)
    <LI>DDDMP_RELATIVE_1    No info on Then(NodeId) follows, because
    Then(NodeId) = NodeId-1
    <LI>DDDMP_TERMINAL      Then Node is a terminal, no info required (for
    BDDs)
    </UL>
  <LI>Ecompl  : 1 bit, if 1 means complemented edge
  <LI>E       : 2 bits, with codes and meanings as for the Then edge
  </UL>
var-info, Then-info, Else-info (if required) are represented as unsigned 
integer values on a sufficient set of bytes (MSByte first).
              ]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
static int
StoreNodeRecur(
DdManager *dd     /* dd manager */,
DdNode *f         /* dd node to be stored */,
int mode          /* store mode */,
int *supportids   /* internal ids for variables */,
char **varnames   /* names of variables: to be stored with nodes */,
int *outids       /* output ids for variables */,
FILE *fp          /* store file */
)
{
  DdNode      *T, *E;
  int idf, idT, idE, vf, vT, vE;
  int retval, diff, var;
  struct binary_dd_code code;
  int nvars = dd->size;

#ifdef DDDMP_DEBUG
  assert(!Cudd_IsComplement(f));
  assert(f!=NULL);
  assert(supportids!=NULL);
#endif

  /* If already visited, nothing to do. */
  if (DddmpVisited(f))
    return(1);

  /* Mark node as visited. */
  DddmpSetVisited(f);

  if (Cudd_IsConstant(f)) {
    /* Check for special case: don't recur */
    idf = DddmpReadNodeIndex(f);
  }
  else {

#ifdef DDDMP_DEBUG
    /* BDDs! Only one constant supported */
    assert(!cuddIsConstant(f));
#endif

    /* 
     * Recursive call for Then edge
     */
    T = cuddT(f);
#ifdef DDDMP_DEBUG
    /* ROBDDs! No complemented Then edge */
    assert(!Cudd_IsComplement(T)); 
#endif
    /* recur */
    retval = StoreNodeRecur(dd,T,mode,supportids,varnames,outids,fp);
    if (retval != 1) return(retval);

    /* 
     * Recursive call for Else edge
     */
    E = Cudd_Regular(cuddE(f));
    retval = StoreNodeRecur(dd,E,mode,supportids,varnames,outids,fp);
    if (retval != 1) return(retval);

    /* 
     * Obtain nodeids and variable ids of f, T, E 
     */

    idf = DddmpReadNodeIndex(f);
    vf = f->index;

    idT = DddmpReadNodeIndex(T);
    if Cudd_IsConstant(T)
      vT = nvars;
    else
      vT = T->index;

    idE = DddmpReadNodeIndex(E);
    if Cudd_IsConstant(E)
	vE = nvars;
    else
	vE = E->index;

  }

  switch (mode) {

    case DDDMP_MODE_TEXT:

      if (Cudd_IsConstant(f)) {
        if (f == Cudd_ReadOne(dd)) {
	    if ((varnames != NULL)||(outids != NULL))
            retval = fprintf (fp,"%d T 1 0 0\n",idf);
          else
            retval = fprintf (fp,"%d 1 0 0\n",idf);
	  }
        else if (f == Cudd_ReadZero(dd)) {
	    if ((varnames != NULL)||(outids != NULL))
            retval = fprintf (fp,"%d T 0 0 0\n",idf);
          else
            retval = fprintf (fp,"%d 0 0 0\n",idf);
	  }
        else {
          /* a constant node different from 1: an ADD constant */
	    if ((varnames != NULL)||(outids != NULL))
	      retval = fprintf (fp,"%d T %-9g 0 0\n",idf,(*Cudd_V(f)).get_min());
	    else
	      retval = fprintf (fp,"%d %-9g 0 0\n",idf, (*Cudd_V(f)).get_min());
        }
      }
	else {
        if (Cudd_IsComplement(cuddE(f)))
	    idE = -idE;
	  if (varnames != NULL) {   
          retval = fprintf (fp,"%d %s %d %d %d\n",
                            idf,varnames[vf],supportids[vf],idT,idE);
        }
	  else {
  	    if (outids != NULL) {   
            retval = fprintf (fp,"%d %d %d %d %d\n",
                            idf,outids[vf],supportids[vf],idT,idE);
          }
          else 
            retval = fprintf (fp,"%d %d %d %d\n",
                            idf,supportids[vf],idT,idE);
	  }
	}

	break;

    case DDDMP_MODE_BINARY:

	/* only integer ids used, varnames ignored */

      if (Cudd_IsConstant(f)) {
        if (f == Cudd_ReadOne(dd)) {
          /*
           * Terminal one is coded as DDDMP_TERMINAL, all other fields are 0
           */
          code.Unused = 0;
          code.V = DDDMP_TERMINAL;
          code.T = 0;
          code.E = 0;
          code.Ecompl = 0;
          retval = DddmpWriteCode (fp,code);
	    if (retval == EOF) return(0);
    	  }
   	  else {
          (void) fprintf (stdout,"DddmpStoreNodeRecur: Binary mode not supported cor ADDs\n");
           return (0);
     	  }
  	}
	else {
        /*
         * Non terminal: output variable id
         */
	  var = supportids[vf];
	  diff = (supportids[vT]<supportids[vE]) ? 
			  (supportids[vT]-var) : (supportids[vE]-var);
        code.V = DDDMP_ABSOLUTE_ID;
        if (diff <= var) {
	    if (diff == 1)
              code.V = DDDMP_RELATIVE_1;
	    else {
              code.V = DDDMP_RELATIVE_ID;
              var = diff;
	    } 
	  }

	  if (T == DD_ONE(dd))
            code.T = DDDMP_TERMINAL;
	  else {
	    /* compute displacement */
	    diff = idf-idT;
            code.T = DDDMP_ABSOLUTE_ID;
	    if (diff <= idT) {
	      if (diff == 1)
                code.T = DDDMP_RELATIVE_1;
	      else {
                code.T = DDDMP_RELATIVE_ID;
		idT = diff;
	      } 
	    }
	  }

	  if (E == DD_ONE(dd))
            code.E = DDDMP_TERMINAL;
	  else {
	    /* compute displacement */
	    diff = idf-idE;
            code.E = DDDMP_ABSOLUTE_ID;
	    if (diff <= idE) {
	      if (diff == 1)
                code.E = DDDMP_RELATIVE_1;
	      else {
                code.E = DDDMP_RELATIVE_ID;
		idE = diff;
	      } 
	    }
	  }
        if (Cudd_IsComplement(cuddE(f)))
          code.Ecompl = 1;
        else
          code.Ecompl = 0;

        retval = DddmpWriteCode (fp,code);
	  if (retval == EOF)
 	    return(0);

        if ((code.V == DDDMP_ABSOLUTE_ID) || 
            (code.V == DDDMP_RELATIVE_ID)) { 
          retval = DddmpWriteInt (fp,var);
	    if (retval == EOF)
	      return(0);
        }

        if ((code.T == DDDMP_ABSOLUTE_ID) || 
            (code.T == DDDMP_RELATIVE_ID)) { 
  	    retval = DddmpWriteInt(fp,idT);
	    if (retval == EOF)
	      return(0);
        }

        if ((code.E == DDDMP_ABSOLUTE_ID) || 
            (code.E == DDDMP_RELATIVE_ID)) { 
	    retval = DddmpWriteInt(fp,idE);
  	    if (retval == EOF)
	      return(0);
        }

	}

	break;

      default:
	return(0);
    }

    if (retval == EOF) {
	return(0);
    } else {
	return(1);
    }

} /* end of StoreNodeRecur */


/**Function********************************************************************

  Synopsis    [String compare for qsort]

  Description [String compare for qsort]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
static int
QsortStrcmp(
const void *ps1,
const void *ps2
)
{
  return strcmp (*((char**)ps1),*((char **)ps2));
}

/**Function********************************************************************

  Synopsis    [Performs binary search of a name within a sorted array]

  Description [Binary search of a name within a sorted array of strings.
               used when matching names of variables.
              ]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
static int
FindVarname(
char *name,
char **array,
int n
)
{
  int d, m, u, t;

  d = 0; u = n-1;

  while (u>=d) {
    m = (u+d)/2;
    t=strcmp(name,array[m]);
    if (t==0)
      return m;
    if (t<0)
      u=m-1;
    else
      d=m+1;
  }
  return -1;
}


/**Function********************************************************************

  Synopsis    [Duplicates a string]

  Description [Allocates memory and copies source string] 

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
static char *
DddmpStrDup (
  char *str       /* string to be duplicated */
)
{
  char *str2;

  str2 = DDDMP_ALLOC(char,strlen(str)+1);
  if (str2 != NULL) {
    strcpy (str2,str);
  }

  return str2;
}

/**Function********************************************************************
  Synopsis    [Duplicates an array of strings]
  Description [Allocates memory and copies source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static char **
DddmpStrArrayDup (
  char **array       /* array of strings to be duplicated */,
  int n              /* size of the array */
)
{
  char **array2;
  int i;

  array2 = DDDMP_ALLOC(char *, n);
  if (array2 == NULL) {
    (void) fprintf (stdout,"DddmpStrArrayDup: Error allocating memory\n");
    return NULL;
  }
  /*
   * initialize all slots to NULL for fair FREEing in case of failure
   */
  for (i=0; i<n; i++) 
    array2[i] = NULL;

  for (i=0; i<n; i++) { 
    if (array[i] != NULL) {
      if ((array2[i]=DddmpStrDup(array[i]))==NULL) {
        DddmpStrArrayFree(array2,n);
        return (NULL);
      }
    }
  }

  return array2;
}


/**Function********************************************************************
  Synopsis    [Inputs an array of strings]
  Description [Allocates memory and inputs source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static char **
DddmpStrArrayRead (
  FILE *fp           /* input file */,
  int n              /* size of the array */
)
{
  char buf[DDDMP_MAXSTRLEN];
  char **array;
  int i;

  assert(fp!=NULL);

  array = DDDMP_ALLOC(char *, n);
  if (array == NULL) {
    (void) fprintf (stdout,"DddmpStrArrayRead: Error allocating memory\n");
    return NULL;
  }
  /*
   * initialize all slots to NULL for fair FREEing in case of failure
   */
  for (i=0; i<n; i++) 
    array[i] = NULL;

  for (i=0; i < n; i++) { 
    if (fscanf (fp, "%s", buf)==EOF) {
      (void) fprintf (stdout,"DddmpStrArrayRead: Error reading file - EOF found\n");
      DddmpStrArrayFree(array,n);
      return (NULL);
    }
    if ((array[i]=DddmpStrDup(buf))==NULL) {
      DddmpStrArrayFree(array,n);
      return (NULL);
    }
  }

  return array;
}

/**Function********************************************************************
  Synopsis    [Outputs an array of strings]
  Description [Allocates memory and inputs source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static int
DddmpStrArrayWrite (
  FILE *fp           /* output file */,
  char **array       /* array of strings */,
  int n              /* size of the array */
)
{
  int i;

  assert(fp!=NULL);

  for (i=0; i<n; i++) { 
    if (fprintf(fp," %s", array[i]) == EOF) {
      (void) fprintf (stdout,"DddmpStrArrayWrite: Error writing to file\n");
      return EOF;
    }
  }

  return n;
}


/**Function********************************************************************
  Synopsis    [Frees an array of strings]
  Description [Frees memory for strings and the array of pointers] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static void
DddmpStrArrayFree (
  char **array       /* array of strings */,
  int n              /* size of the array */
)
{
  int i;

  if (array == NULL)
    return;

  for (i=0; i<n; i++) 
    DDDMP_FREE(array[i]);

  DDDMP_FREE(array);
}


/**Function********************************************************************
  Synopsis    [Duplicates an array of ints]
  Description [Allocates memory and copies source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static int *
DddmpIntArrayDup (
  int *array         /* array of ints to be duplicated */,
  int n              /* size of the array */
)
{
  int *array2;
  int i;

  array2 = DDDMP_ALLOC(int, n);
  if (array2 == NULL) {
    (void) fprintf (stdout,"DddmpIntArrayDup: Error allocating memory\n");
    return NULL;
  }

  for (i=0; i<n; i++) { 
    array2[i] = array[i];
  }

  return array2;
}


/**Function********************************************************************
  Synopsis    [Inputs an array of ints]
  Description [Allocates memory and inputs source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static int *
DddmpIntArrayRead (
  FILE *fp           /* input file */,
  int n              /* size of the array */
)
{
  int *array;
  int i;

  assert(fp!=NULL);

  array = DDDMP_ALLOC(int, n);
  if (array == NULL) {
    (void) fprintf (stdout,"DddmpIntArrayRead: Error allocating memory\n");
    return NULL;
  }

  for (i=0; i < n; i++) { 
    if (fscanf (fp, "%d", &array[i])==EOF) {
      (void) fprintf (stdout,"DddmpIntArrayRead: Error reading file - EOF found\n");
      DDDMP_FREE(array);
      return (NULL);
    }
  }

  return array;
}

/**Function********************************************************************
  Synopsis    [Outputs an array of strings]
  Description [Allocates memory and inputs source array] 
  SideEffects [None]
  SeeAlso     []
******************************************************************************/
static int
DddmpIntArrayWrite (
  FILE *fp           /* output file */,
  int  *array        /* array of ints */,
  int n              /* size of the array */
)
{
  int i;

  assert(fp!=NULL);

  for (i=0; i<n; i++) { 
    if (fprintf(fp," %d", array[i]) == EOF) {
      (void) fprintf (stdout,"DddmpIntArrayWrite: Error writing to file\n");
      return EOF;
    }
  }

  return n;
}

