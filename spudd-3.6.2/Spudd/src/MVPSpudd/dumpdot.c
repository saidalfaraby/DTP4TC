/***********************************************************************
Dumpdot: to dump out the proper dot files for Spudd output (with nice colors)
  $Author: jhoey $
  $Date: 2004/07/05 23:07:33 $
  $Revision: 1.8 $
  $Log: dumpdot.c,v $
  Revision 1.8  2004/07/05 23:07:33  jhoey
  fixed ALLACTIONS so now you can get only one action at each state
  this option (undefining ALLACTIONS) *must* be used for numactions > 50 or so

  Revision 1.7  2004/03/18 16:42:13  jhoey
  matlab interface to server

  Revision 1.6  2004/03/12 19:56:30  jhoey
  *** empty log message ***

  Revision 1.5  2003/11/14 01:17:34  jhoey
  *** empty log message ***

  Revision 1.4  2003/06/20 03:40:49  jhoey
  fixed a little error in roundOffMMA

  Revision 1.3  2003/05/08 23:47:01  jhoey
  changed so it gives all optimal actions

  Revision 1.2  2003/05/07 00:17:27  jhoey
  Now the policies include all possible actions

  Revision 1.1  2003/04/21 20:32:47  jhoey
  new MDP class I put this in a separate project, which I hope will
  become the main branch, leaving MDPSpudd behind

  Revision 2.0  2003/01/16 22:43:59  jhoey
  got everything working now with reading/writing, mutli-valued vars
  it also works with the original data style

  Revision 1.1  2002/12/05 06:43:07  jhoey
  added mutli-valued variables to PSpudd

  Revision 1.3  2002/09/25 20:35:20  jhoey
  small changes for gcc 3.1
  updated dumpdot to make red lines for else arcs

  Revision 1.7  2002/09/24 18:56:50  jhoey
  added in ALLACTIONS option. Must change the define in spudd.hh and in dumpdot.c

  Revision 1.6  2002/02/15 20:41:20  jhoey
  changed ratio to 1

  Revision 1.4  2001/04/07 00:00:07  staubin
  ADDING THE NEW SPUDD DIRECTORY!

  Revision 1.2  2000/10/05 21:26:09  jhoey
  added cvs authorship stuff at header

 
  modified version of Cudd_DumpDot to change the output a bit - 
   changes are flagged with **JH99

***********************************************************************/
#include "MDP.h"

//int aconvert(char * names, char **lnames, double lval, char *separator="\\n");

/* helper function to convert actions in code back to strings */
int aconvert(char * names, char **lnames, double lval, char *separator)
{
  /* here lval will be an integer of 1s and 0s where a 1 in position i means
     that action i is chosen at this leaf */
  /*lval is now a double from 1 to actionnumber*/
  int i,j=0,first=1;
  char temp[64];
  double tval = lval;
  //int ilval = ((int) lval);
  //fprintf(stderr,"lval %f ilval %d\n",lval,ilval);
  if (lval == 0) {
    strcpy(names,"none");
  } else {
#ifdef ALLACTIONS
    //    while (ilval >= 1) {
    //i = ilval%2;
    while (tval >= 1) {
      i = ((int) round(tval-10.0*floor(tval/10.0)))%2;
      //fprintf(stderr,"i %d tval %f j %d\n",i,tval,j);
      if (i) {
	if (!first)
	  strcat(names,separator);
	sprintf(temp,"%s",lnames[j]);
	strcat(names,temp);
	if (first) 
	  first = 0;
      }
      //ilval = ilval >> 1;
      tval = round(tval/2.0-0.0001);
      j++;
    }
#else
    sprintf(temp,"%s",lnames[(int) (lval-1)]);
    strcat(names,temp);
#endif
  }
  return 1;
}
/**Function********************************************************************

  Synopsis    [Writes a dot file representing the argument DDs.]

  Description [Writes a file representing the argument DDs in a format
  suitable for the graph drawing program dot.
  It returns 1 in case of success; 0 otherwise (e.g., out-of-memory,
  file system full).
  Cudd_DumpDot does not close the file: This is the caller
  responsibility. Cudd_DumpDot uses a minimal unique subset of the
  hexadecimal address of a node as name for it.
  If the argument inames is non-null, it is assumed to hold the pointers
  to the names of the inputs. Similarly for onames.
  Cudd_DumpDot uses the following convention to draw arcs:
    <ul>
    <li> solid black line: THEN arcs;
    <li> dotted line: complement arcs;
    <li> red line: regular ELSE arcs.
    </ul>
  The dot options are chosen so that the drawing fits on a letter-size
  sheet.
  ]

  SideEffects [None]

  SeeAlso     [Cudd_DumpBlif Cudd_PrintDebug Cudd_DumpDDcal
  Cudd_DumpDaVinci Cudd_DumpFactoredForm]

******************************************************************************/
int
My_DumpDot(
  DdManager *dd /* manager */,
  int  n /* number of output nodes to be dumped */,
  DdNode **f /* array of output nodes to be dumped */,
  char ** inames /* array of input names (or NULL) */,
  char ** onames /* array of output names (or NULL) */,
  char ** lnames /* array of leaf names or NULL */, 
  FILE * fp /* pointer to the dump file */)
{
    DdNode	*support = NULL;
    DdNode	*scan, *var0, *var6, *hucandg;
    int		*sorted = NULL;
    int		nvars = dd->size;
    st_table	*visited = NULL;
    st_generator *gen = NULL;
    int		retval;
    int		i, j;
    int		slots;
    DdNodePtr	*nodelist;
    long	refAddr, diff, mask;
    //char        namestr[MAXLIN];
    char        namestr[1024];

    /* Build a bit array with the support of f. */
    sorted = ALLOC(int,nvars);
    if (sorted == NULL) {
	dd->errorCode = CUDD_MEMORY_OUT;
	goto failure;
    }
    for (i = 0; i < nvars; i++) sorted[i] = 0;

    /* Take the union of the supports of each output function. */
    support = Cudd_VectorSupport(dd,f,n);
    if (support == NULL) goto failure;
    Cudd_Ref(support);
    scan = support;
    while (!cuddIsConstant(scan)) {
	sorted[scan->index] = 1;
	scan = cuddT(scan);
    }
    Cudd_RecursiveDeref(dd,support);
    support = NULL; /* so that we do not try to free it in case of failure */

    /* Initialize symbol table for visited nodes. */
    visited = st_init_table(st_ptrcmp, st_ptrhash);
    if (visited == NULL) goto failure;

    /* Collect all the nodes of this DD in the symbol table. */

    for (i = 0; i < n; i++) {
	retval = cuddCollectNodes(Cudd_Regular(f[i]),visited);
	if (retval == 0) goto failure;
    }

    /* Find how many most significant hex digits are identical
    ** in the addresses of all the nodes. Build a mask based
    ** on this knowledge, so that digits that carry no information
    ** will not be printed. This is done in two steps.
    **  1. We scan the symbol table to find the bits that differ
    **     in at least 2 addresses.
    **  2. We choose one of the possible masks. There are 8 possible
    **     masks for 32-bit integer, and 16 possible masks for 64-bit
    **     integers.
    */

    /* Find the bits that are different. */
    refAddr = (long) Cudd_Regular(f[0]);
    diff = 0;
    gen = st_init_gen(visited);
    if (gen == NULL) goto failure;
    while (st_gen(gen, (char **) &scan, NULL)) {
	diff |= refAddr ^ (long) scan;
    }
    st_free_gen(gen); gen = NULL;

    /* Choose the mask. */
    for (i = 0; (unsigned) i < 8 * sizeof(long); i += 4) {
	mask = (1 << i) - 1;
	if (diff <= mask) break;
    }

    /* Write the header and the global attributes. */
    retval = fprintf(fp,"digraph \"DD\" {\n");
    if (retval == EOF) return(0);
    retval = fprintf(fp,
		     "size = \"7.5,10\"\nratio=1.0;\ncenter = true;\nedge [dir = none];\n");
    if (retval == EOF) return(0);

    /*** JH99

    retval = fprintf(fp,"{ node [shape = plaintext];\n");
    if (retval == EOF) goto failure;
    retval = fprintf(fp,"  edge [style = invis];\n");
    if (retval == EOF) goto failure;

    retval = fprintf(fp,"  \"CONST NODES\" [style = invis];\n");
    if (retval == EOF) goto failure;
    
    for (i = 0; i < nvars; i++) {
        if (sorted[dd->invperm[i]]) {
	    if (inames == NULL) {
		retval = fprintf(fp,"\" %d \" -> ", dd->invperm[i]);
	    } else {
		retval = fprintf(fp,"\" %s \" -> ", inames[dd->invperm[i]]);
	    }
            if (retval == EOF) goto failure;
        }
    }
    
    retval = fprintf(fp,"\"CONST NODES\"; \n}\n");
    if (retval == EOF) goto failure;
 
    */
   
    /* Write the output node subgraph. */
    retval = fprintf(fp,"{ rank = same; node [shape = box, style=filled, color=forestgreen]; edge [style = invis];\n");
    if (retval == EOF) goto failure;
    for (i = 0; i < n; i++) {
	if (onames == NULL) {
	    retval = fprintf(fp,"\"F%d\"", i);
	} else {
	    retval = fprintf(fp,"\"  %s  \"", onames[i]);
	}
	if (retval == EOF) goto failure;
	if (i == n - 1) {
	    retval = fprintf(fp,"; }\n");
	} else {
	    retval = fprintf(fp," -> ");
	}
	if (retval == EOF) goto failure;
    }

    /* Write rank info: All nodes with the same index have the same rank. */
    for (i = 0; i < nvars; i++) {
        if (sorted[dd->invperm[i]]) {
	    retval = fprintf(fp,"{ rank = same; node [shape=ellipse, style=filled, color=cornflowerblue];");
	    if (retval == EOF) goto failure;

	    /***JH99
		if (inames == NULL) {
		retval = fprintf(fp,"\" %d \";\n", dd->invperm[i]);
		} else {
		retval = fprintf(fp,"\" %s \";\n", inames[dd->invperm[i]]);
		}
		
		if (retval == EOF) goto failure;
	    */
	    nodelist = dd->subtables[i].nodelist;
	    slots = dd->subtables[i].slots;
	    for (j = 0; j < slots; j++) {
		scan = nodelist[j];
		while (scan != NULL) {
		    if (st_is_member(visited,(char *) scan)) {
		      /*****JH99
			retval = fprintf(fp,"\"%lx\";\n", (mask & (long) scan) / sizeof(DdNode));*/
			
		      retval = fprintf(fp,"\"%lx\" [label=\" %s \"];\n",(mask & (long) scan) / sizeof(DdNode),inames[dd->invperm[i]]);
			if (retval == EOF) goto failure;
		    }
		    scan = scan->next;
		}
	    }
	    retval = fprintf(fp,"}\n");
	    if (retval == EOF) goto failure;
	}
    }

    /* All constants have the same rank. */
    /** JH99
    retval = fprintf(fp,
	"{ rank = same; \"CONST NODES\";\n{ node [shape = box, color=goldenrod]; ");
	*/
    retval = fprintf(fp,
	"{ rank = same; \n{ node [shape = box, style=filled, color=goldenrod]; ");
    if (retval == EOF) goto failure;
    nodelist = dd->constants.nodelist;
    slots = dd->constants.slots;
    for (j = 0; j < slots; j++) {
	scan = nodelist[j];
	while (scan != NULL) {
	    if (st_is_member(visited,(char *) scan)) {
		retval = fprintf(fp,"\"%lx\";\n", (mask & (long) scan) / sizeof(DdNode));
		if (retval == EOF) goto failure;
	    }
	    scan = scan->next;
	}
    }
    retval = fprintf(fp,"}\n}\n");
    if (retval == EOF) goto failure;

    /* Write edge info. */
    /* Edges from the output nodes. */

    for (i = 0; i < n; i++) {
	if (onames == NULL) {
	    retval = fprintf(fp,"\"F%d\"", i);
	} else {
	    retval = fprintf(fp,"\"  %s  \"", onames[i]);
	}
	if (retval == EOF) goto failure;
	if (Cudd_IsComplement(f[i])) {
	    retval = fprintf(fp," -> \"%lx\" [style = dotted, color=navyblue];\n",
		(mask & (long) f[i]) / sizeof(DdNode));
	} else {
	    retval = fprintf(fp," -> \"%lx\" [style = solid,color=indianred];\n",
		(mask & (long) f[i]) / sizeof(DdNode));
	}
	if (retval == EOF) goto failure;
    }


    /* Edges from internal nodes. */
    for (i = 0; i < nvars; i++) {
        if (sorted[dd->invperm[i]]) {
	    nodelist = dd->subtables[i].nodelist;
	    slots = dd->subtables[i].slots;
	    for (j = 0; j < slots; j++) {
		scan = nodelist[j];
		while (scan != NULL) {
		    if (st_is_member(visited,(char *) scan)) {
			retval = fprintf(fp,
			    "\"%lx\" -> \"%lx\";\n",
			    (mask & (long) scan) / sizeof(DdNode),
			    (mask & (long) cuddT(scan)) / sizeof(DdNode));
			if (retval == EOF) goto failure;
			if (Cudd_IsComplement(cuddE(scan))) {
			    retval = fprintf(fp,
				"\"%lx\" -> \"%lx\" [style = dotted];\n",
				(mask & (long) scan) / sizeof(DdNode),
				(mask & (long) cuddE(scan)) / sizeof(DdNode));
			} else {
			    retval = fprintf(fp,
					     "\"%lx\" -> \"%lx\" [color=red];\n",  // can add style=dashed
				(mask & (long) scan) / sizeof(DdNode),
				(mask & (long) cuddE(scan)) / sizeof(DdNode));
			}
			if (retval == EOF) goto failure;
		    }
		    scan = scan->next;
		}
	    }
	}
    }

    /* Write constant labels. */
    nodelist = dd->constants.nodelist;
    slots = dd->constants.slots;
    for (j = 0; j < slots; j++) {
	scan = nodelist[j];
	while (scan != NULL) {
	    if (st_is_member(visited,(char *) scan)) {
	      if (lnames == NULL) {
		char *str = (*cuddV(scan)).toString();
		retval = fprintf(fp,"\"%lx\" [label = \"%s \"];\n",
				 (mask & (long) scan) / sizeof(DdNode), str);
		free(str);
	      } else {
		strcpy(namestr,"");
		if (aconvert(namestr, lnames, (*cuddV(scan)).get_min())) {
		  retval = fprintf(fp,"\"%lx\" [label = \" %s \"];\n",
				   (mask & (long) scan) / sizeof(DdNode), namestr);
		}
	      }
	      if (retval == EOF) goto failure;
	    }
	    scan = scan->next;
	}
    }

    /* Write trailer and return. */
    retval = fprintf(fp,"}\n");
    if (retval == EOF) goto failure;

    st_free_table(visited);
    FREE(sorted);
    return(1);

failure:
    if (sorted != NULL) FREE(sorted);
    if (support != NULL) Cudd_RecursiveDeref(dd,support);
    if (visited != NULL) st_free_table(visited);
    return(0);

} /* end of Cudd_DumpDot */
