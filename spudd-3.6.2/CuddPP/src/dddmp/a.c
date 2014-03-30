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
DdNode *f, *T, *E;
struct binary_dd_code code;
char buf[DDDMP_MAXSTRLEN];
int varinfo;
int id, size, maxv;
int nnodes, i, j, k, nsuppvars, nroots, maxaux, 
    var, vT, vE, idT, idE;
int  *permsupport = NULL;
int  *ids = NULL;
int  *permids = NULL;
int  *auxids = NULL;
int  *convertids = NULL;
int  *invconvertids = NULL;
int  *rootids = NULL;
int  *invauxids = NULL;
char rmode[3];
char *ddname = NULL;
char **varnames = NULL;
char **sortedvarnames = NULL;
char **rootnames = NULL;
int  nvars, nddvars;
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

  /* START HEADER */

  while (fscanf(fp,"%s",buf)!=EOF) {

    /* comment */
    if (buf[0] == '#') {
      fgets(buf,DDDMP_MAXSTRLEN,fp);
      continue;
    }

    if (buf[0] != '.') {
      (void) fprintf (stdout,"DdLoad Error at\n%s\n", buf);
      (void) fprintf (stdout,"line must begin with '.' or '#'\n");
      goto failure;
    }

    if matchkeywd(buf, ".ver") {    
      /* this not checked so far: only read */
      if (fscanf (fp, "%*s")==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".mode") {    
      if (fscanf (fp, "%s", rmode)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      if (rmode[0] != mode) {
        if (mode == DDDMP_MODE_DEFAULT) {
          mode = rmode[0];
        }
        else {
          (void) fprintf (stdout,"DdLoad Error: mode mismatch\n");
          goto failure;
        }
      }
      continue;
    }
    if matchkeywd(buf, ".varinfo") {    
      if (fscanf (fp, "%d", &varinfo)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".dd") {    
      if (fscanf (fp, "%s", buf)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      ddname = DDDMP_ALLOC(char,sizeof(buf)+1);
      strcpy(ddname,buf);
      continue;
    }
    if matchkeywd(buf, ".nnodes") {
      if (fscanf (fp, "%d", &nnodes)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }
    if matchkeywd(buf, ".nvars") {   
      if (fscanf (fp, "%d", &nvars)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      if (dd->size != nvars) {
        printf ("Warning: total number of dd manager vars doesn't match with writing manager\n");
        printf ("DDM: %d / BDD: %d\n", dd->size, nvars);
      }
      continue;
    }
    if matchkeywd(buf, ".nsuppvars") {
      if (fscanf (fp, "%d", &nsuppvars)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".varnames") {
      varnames = DDDMP_ALLOC(char *,nsuppvars);
      if (varnames == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nsuppvars; i++) { 
        if (fscanf (fp, "%s", buf)==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file - EOF found\n");
          goto failure;
        }
        varnames[i] = DDDMP_ALLOC(char,strlen(buf)+1);
        if (varnames[i] == NULL) {
          (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
          goto failure;
        }
        strcpy(varnames[i],buf);
      }
      continue;
    }

    if matchkeywd(buf, ".ids") {
      ids = DDDMP_ALLOC(int,nsuppvars);
      if (ids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nsuppvars; i++) { 
        if (fscanf (fp, "%d", &ids[i])==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file\n");
          goto failure;
        }
      }
      continue;
    }

    if matchkeywd(buf, ".permids") {
      permids = DDDMP_ALLOC(int,nsuppvars);
      if (permids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nsuppvars; i++) { 
        if (fscanf (fp, "%d", &permids[i])==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file\n");
          goto failure;
        }
      }
      continue;
    }

    if matchkeywd(buf, ".auxids") {
      auxids = DDDMP_ALLOC(int,nsuppvars);
      if (auxids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nsuppvars; i++) { 
        if (fscanf (fp, "%d", &auxids[i])==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file\n");
          goto failure;
        }
      }
      continue;
    }

    if matchkeywd(buf, ".nroots") {
      if (fscanf (fp, "%d", &nroots)==EOF) {
        (void) fprintf (stdout,"DdLoad: Error reading file\n");
        goto failure;
      }
      continue;
    }

    if matchkeywd(buf, ".rootids") {
      rootids = DDDMP_ALLOC(int,nroots);
      if (rootids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nroots; i++) {
        if (fscanf (fp, "%d", &rootids[i])==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file\n");
          goto failure;
        }
      }
      continue;
    }

    if matchkeywd(buf, ".rootnames") {
      rootnames = DDDMP_ALLOC(char *,nroots);
      if (rootnames == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i < nroots; i++) { 
        if (fscanf (fp, "%s", buf)==EOF) {
          (void) fprintf (stdout,"DdLoad: Error reading file\n");
          goto failure;
        }
        rootnames[i] = DDDMP_ALLOC(char,strlen(buf)+1);
        if (rootnames[i] == NULL) {
          (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
          goto failure;
        }
        strcpy(rootnames[i],buf);
      }
      continue;
    }

    if matchkeywd(buf, ".nodes") {
      if (fgets(buf,999,fp)== NULL) {
        (void) fprintf (stdout,"DdLoad: Error reading file\n");
        goto failure;
      }
      break;
    }

  }

  /* END HEADER */

  /*
   * for each variavle in the support, the relative position in the ordering
   * (within the support only) is computed
   */

  permsupport = DDDMP_ALLOC(int,nsuppvars);
  if (permsupport == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }
  for (i=0,k=0; i < nvars; i++) { 
    for (j=0; j < nsuppvars; j++) { 
      if (permids[j] == i) {
        permsupport[j] = k++;
      }
    }
  }
  assert (k==nsuppvars);

  if (varnames != NULL) {
    /*
     *  Varnames are sorted for binary search
     */
    sortedvarnames = DDDMP_ALLOC(char *,nsuppvars);
    if (sortedvarnames == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
    for (i=0; i < nsuppvars; i++) {
      if (varnames[i] == NULL) {
        (void) fprintf (stdout,"DdLoad Error: support variable name missing in file\n");
        goto failure;
      } 
      sortedvarnames[i] = varnames[i];
    }    
    
    qsort((void *)sortedvarnames,nsuppvars,sizeof(char *),QsortStrcmp);
    
  }

  /*
   * convertids is the array used to vonvert variable ids from positional (shrinked)
   * ids used within the DD file. Positions in the file are from 0 to nsuppvars-1.
   */ 

  convertids = DDDMP_ALLOC(int,nsuppvars);
  if (convertids == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }

  again_matchmode:
  switch (varmatchmode) {
    case DDDMP_VAR_MATCHIDS:
      for (i=0; i<nsuppvars; i++)
        convertids[permsupport[i]] = ids[i];
      break;
    case DDDMP_VAR_MATCHPERMIDS:
      for (i=0; i<nsuppvars; i++)
        convertids[permsupport[i]] = Cudd_ReadInvPerm(dd,permids[i]);
      break;
    case DDDMP_VAR_MATCHAUXIDS:
      if (auxids == NULL) {
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
      for (i=0; i<nsuppvars; i++)
        invauxids[varmatchauxids[ids[i]]] = ids[i];
      /* generate convertids array */
      for (i=0; i<nsuppvars; i++) {
        if ((auxids[i]>maxaux) || (invauxids[auxids[i]]<0)) {
          (void) fprintf (stdout,
                   "DdLoad Error: auxid %d not found in DD manager. ID matching forced (%d)\n", 
                   auxids[i], i);
          (void) fprintf (stdout,"Beware of possible overlappings with other variables\n"); 
          convertids[permsupport[i]]=i;
        }
        else
          convertids[permsupport[i]] = invauxids[auxids[i]];
      }
      break;
    case DDDMP_VAR_MATCHNAMES:
      if (varnames == NULL) {
        (void) fprintf (stdout,"DdLoad Error: variable names matching requested\n");
        (void) fprintf (stdout,"but .varnames not found in BDD file\n");
        (void) fprintf (stdout,"Matching IDs forced.\n");
        varmatchmode = DDDMP_VAR_MATCHIDS;
        goto again_matchmode;
      }
      /* generate invaux array */
      invauxids = DDDMP_ALLOC(int,nsuppvars);
      if (invauxids == NULL) {
        (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
        goto failure;
      }
      for (i=0; i<nsuppvars; i++)
        invauxids[i] = -1;
      for (i=0; i<nddvars; i++) {
        if (varmatchnames[i]==NULL) {
          (void) fprintf (stdout,"DdLoad Warning: NULL match variable name (id: %d). Ignored.\n",
                                 i);
        }
        else
          if ((j=FindVarname(varmatchnames[i],sortedvarnames,nsuppvars))>=0) {
            assert(j<nsuppvars);
            invauxids[j] = i;
          }
      }
      /* generate convertids array */
      for (i=0; i<nsuppvars; i++) {
        assert (varnames[i] != NULL);
        j=FindVarname(varnames[i],sortedvarnames,nsuppvars);
        assert((j>=0)&&(j<nsuppvars));
        if (invauxids[j]<0) {
          (void) fprintf (stdout,
              "DdLoad Error: varname %s not found in DD manager. ID matching forced (%d)\n", 
               varnames[i],i);
          convertids[permsupport[i]]=i;
        }
        else
          convertids[permsupport[i]] = invauxids[j];
      }
      break;
    case DDDMP_VAR_COMPOSEIDS:
      for (i=0; i<nsuppvars; i++)
        convertids[permsupport[i]] = varcomposeids[ids[i]];
      break;
  }

  maxv= -1;
  for (i=0; i<nsuppvars; i++)
    if (convertids[i] > maxv)
      maxv = convertids[i];
 
  invconvertids = DDDMP_ALLOC(int,maxv+1);
  if (invconvertids == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }
  for (i=0; i<=maxv; i++)
    invconvertids[i]= -1;
  for (i=0; i<nsuppvars; i++)
    invconvertids[convertids[i]] = i;

  pnodes = DDDMP_ALLOC(DdNode *,(nnodes+1));
  if (pnodes == NULL) {
    (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
    goto failure;
  }

  if (nsuppvars < 256) {
    pvars1byte = DDDMP_ALLOC(unsigned char,(nnodes+1));
    if (pvars1byte == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
  }
  else if (nsuppvars < 0xffff) {
    pvars2byte = DDDMP_ALLOC(unsigned short,(nnodes+1));
    if (pvars2byte == NULL) {
      (void) fprintf (stdout,"DdLoad: Error allocating memory\n");
      goto failure;
    }
  }
  else {
    printf ("DdLoad Error: more than %d variables. Not supported.\n", 0xffff);
    goto failure;
  }

  pnodes[1] = DD_ONE(dd);

  for (i=1; i<=nnodes; i++) {

    if (feof(fp)) {
      (void) fprintf (stdout,"DdLoad Error: unexpected EOF while reading DD nodes\n");
      goto failure;
    }

    switch (mode) {

      case DDDMP_MODE_TEXT:

        switch (varinfo) {
          case DDDMP_VARIDS:
          case DDDMP_VARPERMIDS:
          case DDDMP_VARAUXIDS:
          case DDDMP_VARNAMES:
            if (fscanf(fp,"%d %*s %d %d %d\n", 
                           &id, &var, &idT, &idE) < 4) {
              (void) fprintf (stdout,"DdLoad: Error reading nodes in text mode\n");
              goto failure;
            }
            break;
          case DDDMP_VARDEFAULT:
            if (fscanf(fp,"%d %d %d %d\n", 
                           &id, &var, &idT, &idE) < 4) {
              (void) fprintf (stdout,"DdLoad: Error reading nodes in text mode\n");
              goto failure;
            }
            break;
        }
#ifdef DDDMP_DEBUG
        assert (id == i);
#endif
        if (var==1 && idT==0 && idE==0)
        {
          /* the 1 leaf */
          continue;
        }
        else
        {
#ifdef DDDMP_DEBUG
          assert (idT > 0);
#endif
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
          /* only 1 terminal presently supported: do nothing */    
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
        vT = nsuppvars;
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
      assert (vT<=nsuppvars);
#endif

#ifdef DDDMP_DEBUG
      assert(idE<i);
#endif
      E = pnodes[idE];
      if (cuddIsConstant(E))
        vE = nsuppvars;
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
      assert (vE<=nsuppvars);
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
      assert (var<nsuppvars);
#endif

      break;

    }

    if (pvars1byte != NULL)
      pvars1byte[i] = (unsigned char) var;
    else if (pvars2byte != NULL)
      pvars2byte[i] = (unsigned short) var;

    var = convertids[var]; 
    pnodes[i] = Cudd_bddIte(dd, Cudd_bddIthVar(dd,var), T, E);
    cuddRef(pnodes[i]);

  }

  fgets(buf, 999,fp);
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
          if (strcmp(rootmatchnames[i],rootnames[j])==0)
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
    id = rootids[i];
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

  for(i=2; i<=nnodes; ++i) { 
    f = pnodes[i];
    Cudd_RecursiveDeref(dd, f);
  }

  /*
   * now free everithing was allocated within this function
   */

load_end:

  DDDMP_FREE(pnodes);
  DDDMP_FREE(pvars1byte);
  DDDMP_FREE(pvars2byte);

  DDDMP_FREE(ddname);
  if (varnames!=NULL)
    for (i=0;i<nsuppvars;i++)
      DDDMP_FREE(varnames[i]);
  DDDMP_FREE(varnames);
   /* variable names are not freed because they were shared with varnames */
  DDDMP_FREE(sortedvarnames);

  DDDMP_FREE(ids);
  DDDMP_FREE(permids);
  DDDMP_FREE(auxids);
  DDDMP_FREE(rootids);

  if (rootnames!=NULL)
    for (i=0;i<nroots;i++)
      DDDMP_FREE(rootnames[i]);
  DDDMP_FREE(rootnames);

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


/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/


