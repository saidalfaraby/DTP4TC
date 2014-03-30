#include "MDP.h"

int cntError;
Pair avgError;
double tolerance;
Pair tempOne,tempTwo;
Pair *pZero, *pOne;
DdNode *One, *Zero;
//Global Manager - does the manager *have* to be global? 
DdManager *gbm(NULL);

// for parser hack
MDP *__theMDP;

MDP::MDP() {
  init();
}
// copy constructor
MDP::MDP(MDP *mdp) {
  init();
  // allocate all new structures
  int i,j;
  numvars = mdp->numvars;
  numorigvars = mdp->numorigvars;
  numactions = mdp->numactions;
  structCopy(mdp->vars,mdp->prime_vars,mdp->orig_vars);
  mvinput = mdp->mvinput;
  actionnames = new char*[numactions];
  for (i=0;i <numactions; i++) 
    actionnames[i] = strdup(mdp->actionnames[i]);
  lnames = new char*[numvars];
  for (i=0; i<numvars; i++)
    lnames[i] = strdup(mdp->lnames[i]);
  /*  onames = new char*[numorigvars];
  for (i=0; i<numorigvars; i++)
    onames[i] = strdup(mdp->onames[i]);
  */
  // make a copy of NewPrime
  NewPrime = (DdNode ***)malloc(numactions*(sizeof(DdNode **)));
  for(i=0;i<numactions;i++) {
    NewPrime[i] = (DdNode **)malloc(numorigvars*(sizeof(DdNode*)));
    for (j=0; j<numorigvars; j++) {
      NewPrime[i][j] = mdp->NewPrime[i][j];
      Cudd_Ref(NewPrime[i][j]);
    }
  }
  RewardD = mdp->RewardD;
  Cudd_Ref(RewardD);
  horizon = mdp->horizon;
  discount_factor = mdp->discount_factor;
  bigadd = mdp->bigadd;
  toleranceP = new Pair(*(mdp->toleranceP));
  discount = mdp->discount;
  Cudd_Ref(discount);
  approxMeth = mdp->approxMeth;

  Array = mdp->Array;
  Array1 = mdp->Array1;
  Array2 = mdp->Array2;
  ArraySum = mdp->ArraySum;
}
MDP::MDP(char *inf, double badd) {
  init();
  bigadd = badd;
  if (!readMDP(inf)) {
    fprintf(stderr,"MDP not allocated properly.  This is probably due to a bad filename.\n");
    exit(0);
  }
}

void MDP::init() {
  // check if gbm exists already
  if (gbm == NULL) {
    delete_gbm = true;
    gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,MAXMEM);
  } else {
    delete_gbm = false;
  }

  pZero = new Pair(0.0);
  pOne = new Pair(1.0);

  for (int i=0; i<MAXVARS; i++) {
    orig_vars[i].nvals = 0;
    orig_vars[i].nbvars = 0;
  }

  One = Cudd_addConst(gbm,pOne);
  Cudd_Ref(One);

  Zero = Cudd_addConst(gbm,pZero);
  Cudd_Ref(Zero);
  Cudd_AutodynDisable(gbm);

  Pair half;
  half.set(0.5);
  Half = Cudd_addConst(gbm,&half);
  Cudd_Ref(Half);

  Pair verysmall;
  verysmall.set(-500000000.0);
  VerySmall = Cudd_addConst(gbm,&verysmall);
  Cudd_Ref(VerySmall);

  tolerance = 0.1;
  bigadd = BIGADD;
  maxusememory = 0;

  OptimalPolicy=NULL;
  OptimalValue = NULL;
  ApproximateValue = NULL;
  ApproximatePolicy = NULL;
  ApproximatePolicyValue = NULL;
  ApproximateMidPointValue = NULL;
  toleranceP = NULL;
  for (int i = 0; i < MAXACT; i++) {
    Qf[i] = NULL;
    actionCost[i] = NULL;
    actionCostNoDummy[i] = NULL;
  }
}
// not checked yet
//Destructor
MDP::~MDP() { 
  int i,j;

  //Part of the cleanup method.........
  for(i=0;i<numvars;i++)
    {
      Cudd_RecursiveDeref(gbm,vars[i].add_var);
      Cudd_RecursiveDeref(gbm,prime_vars[numvars-i-1].add_var);
      Cudd_RecursiveDeref(gbm,Array1[i]);
      Cudd_RecursiveDeref(gbm,Array2[i]);
      
    }
  for (i = 0; i < MAXACT; i++) {
    if (actionCost[i]) {
      Cudd_RecursiveDeref(gbm, actionCost[i]);
    }
    if (actionCostNoDummy[i]) {
      Cudd_RecursiveDeref(gbm, actionCostNoDummy[i]);
    }
  }
  for(i=0;i<numactions;i++) {
#ifdef OPTIMIZED
    Cudd_RecursiveDeref(gbm,Allprime[i][0]);  
#endif
    if (mvinput) {
      for(j=0;j<numorigvars;j++) {
	Cudd_RecursiveDeref(gbm,NewPrime[i][j]);  
      }
    } else {
      for(j=0;j<numorigvars;j++) {
	Cudd_RecursiveDeref(gbm,NewPrime[i][j]);  
      }
    }
    
  }
  Cudd_RecursiveDeref(gbm,RewardD);
  if (VerySmall) {
    Cudd_RecursiveDeref(gbm,VerySmall);
  }
  if (Half) {
    Cudd_RecursiveDeref(gbm,Half);
  }
  Cudd_RecursiveDeref(gbm,discount);
  Cudd_RecursiveDeref(gbm,Zero);
  Cudd_RecursiveDeref(gbm,One);
  if (OptimalValue) {
    Cudd_RecursiveDeref(gbm, OptimalValue);
  }
  if (OptimalPolicy) {
    Cudd_RecursiveDeref(gbm, OptimalPolicy);
  }
  if (ApproximateValue) {
    Cudd_RecursiveDeref(gbm, ApproximateValue);
  }
  if (ApproximatePolicy) {
    Cudd_RecursiveDeref(gbm, ApproximatePolicy);
  }
  if (ApproximatePolicyValue) {
    Cudd_RecursiveDeref(gbm, ApproximatePolicyValue);
  }
  if (ApproximateMidPointValue) {
    Cudd_RecursiveDeref(gbm, ApproximateMidPointValue);
  }
  for (i = 0; i < numactions; i++) {
    if (Qf[i]) {
      Cudd_RecursiveDeref(gbm, Qf[i]);
    }
  }

  if(Array)
    free(Array);
  
  if(Array1)
    free(Array1);
    
  if(Array2)
    free(Array2);
  
  if(ArraySum)
    free(ArraySum);
   
  if(Allprime){
    for(i=0;i<numactions;i++)
      free(Allprime[i]);
    free(Allprime);
  }
  if(NewPrime){
    int lim = numactions;
    if (mvinput) {
      lim = MAXACT;
    }
    for(i=0;i<lim;i++)
      free(NewPrime[i]);
    free(NewPrime);
  }
  if (varnames) {
    for (int i = 0; i < 2*numvars; i++) {
      free(varnames[i]);
    }
    free(varnames);
  }
  if (actionnames) {
    for (int i = 0; i < numactions; i++) {
      free(actionnames[i]);
      free(actionlist[i].name);
    }
    free(actionnames);
  }
  if (lnames) {
    for (int i = 0; i < numvars+3; i++) {
      free(lnames[i]);
    }
    free(lnames);
  }
  if (toleranceP) {
    delete toleranceP;
  }
  int curr_index = 0;
  i = 0;
  while (curr_index < numvars) {
    free(orig_vars[i].name);
    for (int j = 0; j < orig_vars[i].nvals; j++) {
      free(orig_vars[i].valname[j]);
    }
    for (int k = 0; k < orig_vars[i].nbvars; k++) {
      free(vars[curr_index].name);
      free(prime_vars[curr_index].name);
      curr_index++;
    }
    i++;
  }

  delete pZero;
  delete pOne;


  if (delete_gbm) {
    Cudd_Quit(gbm);
  }

}
//reads MDP from Dual Optimal ADD file
int MDP::readMDP(char *infile, bool dualOptimal) {
  if (dualOptimal) {
    //read dual optimal
    return 0;
  } else {
    return readMDPSpec(infile);
  }
}
int MDP::getNumActions() {
  return numactions;
}
int MDP::getNumvars() {
  return numvars;
}
int MDP::getNumOrigvars() {
  return numorigvars;
}
onum *MDP::getOrigVars() {
  return orig_vars;
}
rnum *MDP::getVars() {
  return vars;
}
char *MDP::getActionName(int i) {
  return actionnames[i];
}
DdNode *MDP::getNewPrime(int act, int var) {
  return NewPrime[act][var];
}
void MDP::setHorizon(int h)
{
  horizon = h;
}

double MDP::getDiscount() {
  if (!discount || !Cudd_IsConstant(discount)) {
    return 1;
  }
  return Cudd_V(discount)->get_max();
}

void MDP::setNewPrime(int act, int var, DdNode *newNewPrime) {
  Cudd_RecursiveDeref(gbm,NewPrime[act][var]);
  NewPrime[act][var] = newNewPrime;
  Cudd_Ref(NewPrime[act][var]);
}
void decodevvals(int *varvals, onum *orig_vars, int numorigvars, int state)
{
  int k;
  k=numorigvars-1;
  while (k>=0) {
    varvals[k] = state%(orig_vars[k].nvals);
    state = state/(orig_vars[k].nvals);
    k--;
  }
}
void MDP::decodevarvals(int *varvals, int state) {
  decodevvals(varvals,orig_vars,numorigvars,state);
}
void MDP::flattenT(double ***T, int & ns, int & na) {
  int i,j,a;
  Pair rval;
  int numstates(1);
  int *varvals, *varvalsp;
  varvals = new int[numvars];
  varvalsp = new int[numvars];
  for (i=0; i<numorigvars; i++) 
    numstates *= orig_vars[i].nvals;
  *T = new double*[numactions];
  buildAllPrime();
  for (a=0; a<numactions; a++) {
    (*T)[a] = new double[numstates*numstates];

    for (i=0; i<numstates; i++) {
      // decode values of variables from i
      decodevarvals(varvals,i);
      for (j=0; j<numstates; j++) {
	decodevarvals(varvalsp,j);
	
	
	getVal(gbm, Allprime[a][0], rval, varvals,varvalsp);
	
	(*T)[a][i*numstates+j] =  rval.get_min(); 

      }
    }
  }
  ns = numstates;
  na = numactions;
}

int MDP::parseInput() {
  int i;
  /* Parse the input file */  
  numvars =0;
  numactions=0;
  numorigvars = 0;
  //yydebug = 1;
  // only works for mvinput right now
  __theMDP = this;
  if (mvinput) {
    return !yyparse();
  } else {
    if (!pparse(stderr,&discount_factor,&horizon, 
	       &tolerance,&numvars,&numactions,vars,prime_vars,actionlist)) { 
      fprintf(stderr,"parse error\n");
      Cudd_Quit(gbm);
      return 0;
    }
    numorigvars = numvars;

    // copy vars into origvars
    for (i=0; i<numvars; i++) {
      orig_vars[i].name = strdup(vars[i].name);
      orig_vars[i].nvals = 2;
      orig_vars[i].valname[0] = "true";
      orig_vars[i].valname[1] = "false";
      orig_vars[i].nbvars = 1;
      orig_vars[i].var1index = i;
    }
  }
  return 1;
}
// constructs a new ADD variable for the novar original multi-valued variable
void MDP::newADDVar() {
  char tmp[MAXLEN];
  
  // figure out how many new variables to add
  // I realise this could be done by taking logs and ceil, but this is safer
  int i,bvars(1),bnvals(2);
  while (bnvals < orig_vars[numorigvars].nvals) {
    bvars++;
    bnvals = bnvals*2;
  }
  orig_vars[numorigvars].nbvars = bvars;
  orig_vars[numorigvars].var1index = numvars;
  for (i=0; i<bvars; i++) {
    vars[numvars].orig_num = numorigvars;
    vars[numvars].number = numvars;
    prime_vars[numvars].number = numvars;
    //fprintf(stderr,"numorigvars %d %s\n",numorigvars,orig_vars[numorigvars].name);
    sprintf(tmp,"%s%d",orig_vars[numorigvars].name,i);

    vars[numvars].name = strdup(tmp);
    sprintf(tmp,"%sr",vars[numvars].name);

    prime_vars[numvars].name = strdup(tmp);
    prime_vars[numvars].add_var = Cudd_addNewVar(gbm);
    Cudd_Ref(prime_vars[numvars].add_var);
    vars[numvars].add_var = Cudd_addNewVar(gbm);
    Cudd_Ref(vars[numvars].add_var);
    numvars++;
  }
}
//reads MDP from mdp file 
int MDP::readMDPSpec(char *infile) {
  int i;
  char firstchar;

  yyin = fopen(infile,"r");

  if (!yyin) {
    return 0;
  }

  setbuf(stdout,0); 
  

  // check for boolean input file
  mvinput=true;
  fscanf(yyin,"%c",&firstchar);
  if (firstchar == 'v') 
    mvinput = false; 

  // reset the input file
  rewind(yyin);

  if (!parseInput()) {
    fclose(yyin);
    return 0;
  }

  fclose(yyin);

  /* convert the pruning percentage to error by using the span and discount */
  extR = get_extent(RewardD);

  
  /* CREATING THE CONSTANTS FOR THE DISCOUNTING AND THE PRIMING OF THE ACTION */
  /* 
     //done in the parser now
  Pair discountPair;
  discountPair.set(discount_factor);
  discount = Cudd_addConst(gbm,&discountPair);
  Cudd_Ref(discount);
  */
  /* random seed */
  long curtime = time(NULL);
  unsigned seed = (unsigned) curtime % INT_MAX;
  srand(seed);

  //Cudd_EnableReorderingReporting(gbm);

  // ALLOCATE ALL MEMORY
  allocateMemory();

#ifdef OPTIMIZED
  // don't need these unless we call primeReward
  buildAllPrime();
#endif

  // set up some variable names
  for (i=0; i<numvars; i++) {
    strcpy(varnames[2*i+1],vars[i].name);
    strcpy(varnames[2*i],prime_vars[i].name);
    strcpy(lnames[i],vars[i].name);
  }

  sprintf(lnames[numvars],"VALUE");
  sprintf(lnames[numvars+1],"ACTION");
  sprintf(lnames[numvars+2],"REWARD");
  
  for (i=0; i<numactions; i++) {
    strcpy(actionnames[i],actionlist[i].name);
  }
  return 1;
}

void MDP::doshit(int *teststate, int sv, bool doq)  
{
  int a,j;
  DdNode *temp, *res, *res2, *temp1;
  double *v = new double[numactions];
  DdNode **pval = new DdNode*[numactions];
  DdNode **pprob = new DdNode*[numactions];
  if (doq) {
    for (a=0; a< numactions; a++) {
      temp = buildCubeOrig(teststate,false);
      res = Cudd_addApply(gbm,Cudd_addTimes,Qf[a],temp);
      Cudd_Ref(res);
      Cudd_RecursiveDeref(gbm,temp);
      // sum out unprimed variables
      temp = sumOutAllOrigVars(res);
      Cudd_RecursiveDeref(gbm,res);
      res = temp;
      v[a] = (*Cudd_V(res)).get_min();
      Cudd_RecursiveDeref(gbm,res);
    }
  }
  // now check the distributions over next states for sv only
  for (a=0; a< numactions; a++) {
    temp = buildCubeOrig(teststate,false);
    res = Cudd_addApply(gbm,Cudd_addTimes,NewPrime[a][sv],temp);
    Cudd_Ref(res);
    Cudd_RecursiveDeref(gbm,temp);
    // sum out unprimed variables
    temp = sumOutAllOrigVars(res);
    Cudd_RecursiveDeref(gbm,res);
    // and sum out all other primed variables - these must all be further 
    // down in the order than sv
    for (j=sv+1; j<numorigvars; j++) {
      if (ovarsupport[a][sv][j] >= 2) {
	res = sumOutPrime(temp,orig_vars+j,prime_vars);
	Cudd_RecursiveDeref(gbm,temp);
	temp = res;
      }
    }
    pprob[a] = temp;
  }
  // now check the and the values of those distributions over next states for sv only
  // get expected value of the future distributed over sv
  for (a=0; a< numactions; a++) {
    temp1 = Cudd_addSwapVariables(gbm,OptimalValue,Array1,Array2,numvars); 
    Cudd_Ref(temp1); 

    // we're only interested in this state
    temp = buildCubeOrig(teststate,false);
    for (j=0; j<numorigvars; j++) {
	// only interested in test state
        res = Cudd_addApply(gbm,Cudd_addTimes,NewPrime[a][j],temp);
	Cudd_Ref(res);
	res2 = Cudd_addApply(gbm,Cudd_addTimes,res,temp1);
	Cudd_Ref(res2);
	Cudd_RecursiveDeref(gbm,temp1);
	Cudd_RecursiveDeref(gbm,res);
	temp1 = res2;
	if (j != sv) {
	  res2 = sumOutPrime(temp1,orig_vars+j,prime_vars);
	  Cudd_RecursiveDeref(gbm,temp1);
	  temp1 = res2;
	}
    }
    Cudd_RecursiveDeref(gbm,temp);
    // sum out unprimed variables
    temp = sumOutAllOrigVars(temp1);
    Cudd_RecursiveDeref(gbm,temp1);
    pval[a] = temp;
  } 
  fprintf(stderr,"state: \n");
  for (j=0; j<numorigvars; j++) 
    fprintf(stderr,"%s : %s\n",orig_vars[j].name,orig_vars[j].valname[teststate[j]]);
  double *vprob = new double[orig_vars[sv].nvals];
  double *vval = new double[orig_vars[sv].nvals];
  int acttaken = consultOptimalPolicy(teststate);
  fprintf(stderr,"spudd's choice of action is ");
  printAction(acttaken,stderr);
  fprintf(stderr,"\n");
  for (a=0; a< numactions; a++) {
    if (doq)
      fprintf(stderr,"action %s (%d) value %f\n",actionnames[a],a,v[a]);
    else 
      fprintf(stderr,"action %s (%d) \n",actionnames[a],a);
    
    printDDSingleVariable(pprob[a],vprob,sv);
    printDDSingleVariable(pval[a],vval,sv);
    fprintf(stderr,"%s state : probability   expected value\n",orig_vars[sv].name);
    for (j=0; j<orig_vars[sv].nvals; j++) 
      fprintf(stderr,"%s : %f %f\n",orig_vars[sv].valname[j],vprob[j],vval[j]);
    Cudd_RecursiveDeref(gbm,pprob[a]);
    Cudd_RecursiveDeref(gbm,pval[a]);
  }
  delete [] pprob;
  delete [] pval;
  delete [] v;
  delete [] vprob;
  delete [] vval;

}
// does shit for more than one variable
void MDP::doshit2(int *teststate, int nsv, int *sv)  
{
  int a,j,i,k;
  DdNode *temp, *res, *res2;
  DdNode ***pval = new DdNode**[nsv];
  DdNode ***pprob = new DdNode**[nsv];
  DdNode ***npdsv = new DdNode**[nsv];
  // first, get the newprime distributions over each var in sv summed over teststate
  for (j=0; j<nsv; j++) {
    npdsv[j] = new DdNode*[numactions];
    for (a=0; a<numactions; a++) {
      temp = buildCubeOrig(teststate,false);
      res = Cudd_addApply(gbm,Cudd_addTimes,NewPrime[a][sv[j]],temp);
      Cudd_Ref(res);
      Cudd_RecursiveDeref(gbm,temp);
      // sum out unprimed variables
      temp = sumOutAllOrigVars(res);
      Cudd_RecursiveDeref(gbm,res);
      npdsv[j][a] = temp;
    }
  }
  // now check the distributions over next states for sv only
  // start at the lst one 
  int *newdep = new int[nsv];
  for (j=nsv-1; j>=0; j--) {
    pprob[j] = new DdNode*[numactions];
    for (a=0; a< numactions; a++) {
      // and sum out all other primed variables - these must all be further 
      // down in the order than sv
      temp = npdsv[j][a];
      Cudd_Ref(temp);
      for (k=0; k<nsv; k++) 
	newdep[k] = 0;
      for (i=j+1; i<nsv; i++) {
	if (newdep[i] || ovarsupport[a][sv[j]][sv[i]] >= 2) {
	  // see if this variable depends on any higher ones
	  for (k=i+1; k<nsv; k++) {
	    if (ovarsupport[a][sv[i]][sv[k]] >= 2) 
	      newdep[k] = 1;
	  }
	  // multiply by npdsv[i][a]
	  res = Cudd_addApply(gbm,Cudd_addTimes,temp,npdsv[i][a]);
	  Cudd_Ref(res);
	  Cudd_RecursiveDeref(gbm,temp);
	  temp = res;
	  // now sum out i
	  res = sumOutPrime(temp,orig_vars+sv[i],prime_vars);
	  Cudd_RecursiveDeref(gbm,temp);
	  temp = res;
	}
      }
      pprob[j][a] = temp;
    }
  }
  delete [] newdep;
  // now check the and the values of those distributions over next states for sv only
  // get expected value of the future distributed over sv
  DdNode *ovs = Cudd_addSwapVariables(gbm,OptimalValue,Array1,Array2,numvars); 
  Cudd_Ref(ovs); 
  // get value for current state - over all primed variables of interest
  for (i=0; i<nsv; i++) 
    pval[i] = new DdNode*[numactions];
  for (a=0; a< numactions; a++) {
    temp = ovs;
    Cudd_Ref(temp); 
    // get optimal value at primed state predicted by this action's effect
    for (j=0; j<nsv; j++) {
      res = Cudd_addApply(gbm,Cudd_addTimes,npdsv[j][a],temp);
      Cudd_Ref(res);
      temp = res;
    }
    // get the value for each primed variable of interest
    for (i=0; i<nsv; i++) {
      // sum out all other primed variables
      res2 = temp;
      Cudd_Ref(res2);
      for (j=0; j<nsv; j++) {
	if (j != i) {
	  res = sumOutPrime(res2,orig_vars+sv[j],prime_vars);
	  Cudd_RecursiveDeref(gbm,res2);
	  res2 = res;
	}
      }
      pval[i][a] = res2;
    }
  }
  fprintf(stderr,"state: \n");
  for (j=0; j<numorigvars; j++) 
    fprintf(stderr,"%s : %s\n",orig_vars[j].name,orig_vars[j].valname[teststate[j]]);
  double *vprob = new double[MAXVALS];
  double *vval = new double[MAXVALS];
  int acttaken = consultOptimalPolicy(teststate);
  fprintf(stderr,"spudd's choice of action is ");
  printAction(acttaken,stderr);
  fprintf(stderr,"\n");
  for (i=0; i<nsv; i++) {
    fprintf(stderr,"variable %d: %s\n",sv[i],orig_vars[sv[i]].name);
    for (a=0; a< numactions; a++) {
      fprintf(stderr,"action %s (%d) \n",actionnames[a],a);
      printDDSingleVariable(pprob[i][a],vprob,sv[i]);
      printDDSingleVariable(pval[i][a],vval,sv[i]);
      fprintf(stderr,"state  probability   expected value\n");
      for (j=0; j<orig_vars[sv[i]].nvals; j++) 
	fprintf(stderr,"%s : %f %f\n",orig_vars[sv[i]].valname[j],vprob[j],vval[j]);
      Cudd_RecursiveDeref(gbm,pprob[i][a]);
      Cudd_RecursiveDeref(gbm,pval[i][a]);
    }
  }
  for (j=0; j<nsv; j++) {
    delete [] pprob[j];
    delete [] npdsv[j];
    delete [] pval[j];
  }
  delete [] pprob;
  delete [] npdsv;
  delete [] pval;
  delete [] vprob;
  delete [] vval;
  
}
// returns the value of the dd for all the value of ovar
// assumes that only ovar is in dd
void MDP::printDDSingleVariable(DdNode *dd, double *vv, int ovar)
{
  DdNode *temp;
  for (int j=0; j<orig_vars[ovar].nvals; j++) {
    temp = restrictVal(gbm,dd,prime_vars,numvars,orig_vars,numorigvars,ovar,j);
    Cudd_Ref(temp);
    if (Cudd_IsConstant(temp)) {
      vv[j] = (*Cudd_V(temp)).get_min();
    } else {
      vv[j] = 0;
    }
  }

}
int MDP::generatePolicy(int rMeth, int rAMeth, double tol, int badd, float prun, int max_size, 
			int aMeth, int ltype, int tMode, int sFlag, bool dFlag) {
  tolerance = tol;
  bigadd = badd;
  prune = prun;
  Max_Size = max_size;
  approxMeth = aMeth;
  reorderMeth = rMeth;
  reorderApplyMeth = rAMeth;
  limit_type = ltype;
  toleranceMode = tMode;
  shuffleFlag = sFlag;

  if (toleranceP) {
    delete toleranceP;
  }

  toleranceP = new Pair(tolerance);

  if (shuffleFlag) {
    //fprintf(stderr,"Shuffling order...");
    shuffleRandom(); 
    //fprintf(stderr,"done\n");
  }
  
  // reroder according to reward
  reorder(RewardD);

 
  //get Optimal or approximate value function
  getValueFunction();

  return 1;
}

void MDP::getValueFunction() {
  
  DdNode *tempM, *maxADD;
  Pair maxOPTP, maxPAIR;
  double maxApp,minApp;

  // call Spudd() here - this actually gets the value function
  // according to the current set of parameters
  //Spudd();
  // now if OptimalValue exists, start from there!
  Spudd(OptimalValue);

  switch (approxMeth) 
    {
    case APPROX_NONE:      
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      // also sets ApproximateMidPointValue to be the mid-point aDD
      // and ApproximateValue to be the approximate value function (ranged)
      maxApp = (*Cudd_V(Cudd_addFindMax(gbm,ApproximateMidPointValue))).get_max();
      minApp = (*Cudd_V(Cudd_addFindMin(gbm,ApproximateMidPointValue))).get_min();
      
      maxPAIR.set(1/(maxApp-minApp));
      maxADD = Cudd_addConst(gbm,&maxPAIR);  
      Cudd_Ref(maxADD);
      
      tempM = Cudd_addApply(gbm,Cudd_addTimes,ApproximateMidPointValue,maxADD);
      Cudd_Ref(tempM);
      
      Cudd_RecursiveDeref(gbm,ApproximateMidPointValue);
      ApproximateMidPointValue = tempM;
      Cudd_Ref(ApproximateMidPointValue);
      Cudd_RecursiveDeref(gbm,tempM);
      Cudd_RecursiveDeref(gbm,maxADD);      
      break;
    default:
      break;
    }
}
DdNode *MDP::modifiedPolicyIteration(bool fromOpt)
{
  approxMeth = APPROX_NONE;
  toleranceP = new Pair(tolerance);

  DdNode *tmp;
  // assumes OptimalPolicy and OptimalValue already exist
  if (fromOpt) {
    tmp = OptimalPolicy;
    Cudd_Ref(tmp);
  } else {
    tmp = One;
    Cudd_Ref(tmp);
    OptimalValue = RewardD;
    Cudd_Ref(OptimalValue);
  }
  DdNode *tmp2 = modPolicyIteration(&OptimalValue,tmp,5);
  Cudd_RecursiveDeref(gbm,tmp);

  // testing
  //tmp = Cudd_addApply(gbm,BinaryAction,tmp2,One);
  //Cudd_Ref(tmp);

  OptimalPolicy = Cudd_addApply(gbm,BinaryAction,tmp2,One);
  Cudd_Ref(OptimalPolicy);
  
  Cudd_RecursiveDeref(gbm,tmp2);
  return NULL; 
}
// picks an action from each leaf of the policy
// assumes OptimalPolicy exists
void MDP::pickPolicyAction() 
{
  // newPol has a single action (a number from 0...numactions) at each leaf of the policy
  DdNode *newPol = Cudd_addApply(gbm,PickAction,OptimalPolicy,One);
  Cudd_Ref(newPol);

  Cudd_RecursiveDeref(gbm,OptimalPolicy);
  
  // convert back to binary format
  OptimalPolicy = Cudd_addApply(gbm,BinaryAction,newPol,One);
  Cudd_Ref(OptimalPolicy);
  Cudd_RecursiveDeref(gbm,newPol);
}

// modified policy iteration
// numiters is the number of VI iterations to do at each MPI iteration
DdNode *MDP::modPolicyIteration(DdNode **valD, DdNode *actD, int numiters)
{
  // first, restrict each action in actD to a single value
  DdNode *newPol, *Vpast, *tmp;

  // newPol has a single action (a number from 0...numactions) at each leaf of the policy
  newPol = Cudd_addApply(gbm,PickAction,actD,One);
  Cudd_Ref(newPol);
  
  

  // loop until value function converges
  int counter = -1;
  int iters(0);
  bool converged = false;
  while (!converged) {
    Vpast = *valD;
    Cudd_Ref(Vpast);
    
    if (counter == -1) {
      fprintf(stderr,"improving policy %d\n",iters++);
      newPol = improvePolicy(valD);
    } else {
      fprintf(stderr,"evaluating %d\n",counter);
      tmp = evaluatePolicy(newPol, *valD);
      Cudd_RecursiveDeref(gbm,*valD);
      *valD = tmp;
    }

    counter++;
    if (counter == numiters)
      counter = -1;

    converged = Cudd_EqualSupNorm(gbm,*valD,Vpast,toleranceP,0);
    Cudd_RecursiveDeref(gbm,Vpast);
  }

  return newPol;
}

// evaluates (backs up) policy actD for 1 iteration
// starting with value function valD
// returns the result pre-reffed
DdNode *MDP::evaluatePolicy(DdNode *actD, DdNode *valD)
{
  int a;
  DdNode *tmp,*tmp2,*tmp3;
  DdNode *newValD = valD;
  Cudd_Ref(newValD);

  DdNode **Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));
  for(a=0;a<numactions;a++) {
    Vcurrent[a] = valD;
    Cudd_Ref(Vcurrent[a]);
  }  

  tmp = Cudd_addSwapVariables(gbm,newValD,Array1,Array2,numvars);
  Cudd_Ref(tmp);
  
  for (a=0; a<numactions; a++) {
    tmp2 = multiplySumSet(tmp,NewPrime[a]);
    tmp3 = Cudd_addApply(gbm,Cudd_addTimes,tmp2,discount);
    Cudd_Ref(tmp3);
    Cudd_RecursiveDeref(gbm,tmp2);
    
    Cudd_RecursiveDeref(gbm,Vcurrent[a]);
    Vcurrent[a] = Cudd_addApply(gbm,Cudd_addPlus,tmp3,actionCost[a]);
    Cudd_Ref(Vcurrent[a]);
    Cudd_RecursiveDeref(gbm,tmp3);
  }
  tmp2 = actionMerge(Vcurrent, actD, numactions);
  Cudd_Ref(tmp2);
  
  newValD = Cudd_addApply(gbm,Cudd_addPlus,tmp2,RewardD);
  Cudd_Ref(newValD);
  Cudd_RecursiveDeref(gbm,tmp2);
  Cudd_RecursiveDeref(gbm,tmp);

  for(a=0;a<numactions;a++) 
    Cudd_RecursiveDeref(gbm,Vcurrent[a]);
  free(Vcurrent);

  return newValD;
}
// backs up value function valD once
// returns maximal actions (policy) and replaces valD
// with new value function
// returned policy is pre-reffed
DdNode *MDP::improvePolicy(DdNode **valD) {
  int i;
  DdNode *tmp, *tmp2;
  
  DdNode *nnValD = Cudd_addSwapVariables(gbm,*valD,Array1,Array2,numvars); 
  Cudd_Ref(nnValD);
  DdNode *actionDD = Zero;
  Cudd_Ref(actionDD);

  DdNode *newValD = VerySmall;
  Cudd_Ref(newValD);

  for(i=0;i<numactions;i++) { 
      
    // do the sum_t(Pr(s,a,t)*V(t)) 
    tmp = multiplySumSet(nnValD,NewPrime[i]);
    tmp2 = Cudd_addApply(gbm,Cudd_addTimes,tmp,discount);
    Cudd_Ref(tmp2);
      
    Cudd_RecursiveDeref(gbm,tmp);
    tmp = Cudd_addApply(gbm,Cudd_addPlus,tmp2,RewardD);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(gbm,tmp2);
      
    // add the action costs. 
    tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,actionCost[i]);
    Cudd_Ref(tmp2);

    Cudd_RecursiveDeref(gbm,tmp);
    
    /* maximize over actions */
    tmp=Cudd_addApply(gbm,Cudd_addMaximum,newValD,tmp2);
    Cudd_Ref(tmp);

    // add this action to actionDD at each node where tmp2 == tmp
    // replace it when tmp2 > tmp
    recoverPolicy(tmp2,tmp,newValD,&actionDD,i);
    Cudd_RecursiveDeref(gbm,tmp2);
    Cudd_RecursiveDeref(gbm,newValD);
    newValD = tmp;
  } // end of loop over actions
  // newValD is the final value function
  // actionDD is final policy
  Cudd_RecursiveDeref(gbm,*valD);
  Cudd_RecursiveDeref(gbm,nnValD);
  *valD = newValD;
  
  // pick one action from actionDD
  tmp = Cudd_addApply(gbm,PickAction,actionDD,One);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,actionDD);
  actionDD = tmp;
  return actionDD;
}

/*
Loop like policy Iteration....
Backs up the value function valD using the
policy specified by actD
*/
DdNode *MDP::evalPolicy(DdNode *valD, DdNode *actD)
{
  int i,test,counter;
  DdNode **Vcurrent;
  DdNode *Vpast,*temp,*temp1, *temp2;
  DdNode *Merged_Add;
  
  Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));
  
  test = 0;
    
  for(i=0;i<numactions;i++) {
    Vcurrent[i] = valD;
    Cudd_Ref(Vcurrent[i]);
  }  
  
  Vpast = valD;
  Cudd_Ref(Vpast);
  
  Merged_Add = valD;
  Cudd_Ref(Merged_Add);
  
  counter = 0;
  while ((horizon == -1 && test==0) || (horizon >= 0 && counter < horizon))
  { 
    // KEEP TRACK OF THE NUMBER OF ITERATIONS 
      counter = counter + 1;
      //fprintf(stderr,"iteration %d\n",counter);
      // KEEPING A COPY OF THE LAST ADD TO COMPARE FOR CONVERGENCE 
      Cudd_RecursiveDeref(gbm,Vpast);
      Vpast = Merged_Add;
      Cudd_Ref(Vpast);

      temp = Cudd_addSwapVariables(gbm,Merged_Add,Array1,Array2,numvars);
      Cudd_Ref(temp);
      for(i=0;i<numactions;i++){
#ifdef OPTIMIZED
	temp1 = primeReward(temp,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
	temp1 = multiplySumSet(temp,NewPrime[i]);
#endif
          
	temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp1);
	  
	temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
	Cudd_Ref(temp1);
	Cudd_RecursiveDeref(gbm,temp2);
	
	Cudd_RecursiveDeref(gbm,Vcurrent[i]);
	Vcurrent[i] = temp1;
	
      }
      Cudd_RecursiveDeref(gbm,temp);
      Merged_Add = actionMerge(Vcurrent, actD, numactions);
      Cudd_Ref(Merged_Add);
      // COMPARE THE TWO MOST RECENT ADDs IN ORDER TO VERIFY THE CONVERGENCE 
      //test = Merged_Add.EqualSupNorm(Vpast,toleranceP,0);
      test= Cudd_EqualSupNorm(gbm,Merged_Add,Vpast,toleranceP,0);
    }
  
  for (i=0; i<numactions; i++)
    Cudd_RecursiveDeref(gbm, Vcurrent[i]);
  free(Vcurrent);
  
  return Merged_Add;
}

void MDP::writePolicyToFile(char *outputfile, bool dotFlag) {
  char basepath[MAXLINE];
  int numwrite;
  switch (approxMeth) 
    {
    case APPROX_NONE: 
      //save the value add and the policy to a file. 

      // debugging
#ifdef DEBUG
      Cudd_PrintDebug(gbm,OptimalPolicy,4,100);
      Cudd_PrintDebug(gbm,OptimalValue,4,100);
#endif

       if (dotFlag) {
       	sprintf(basepath,"%s-OPT",outputfile);
       	writeVPDot(OptimalValue,OptimalPolicy,basepath);
       }
      sprintf(basepath,"%s-OPTDual.ADD",outputfile);
      numwrite = writeDualOptimal(basepath,OptimalPolicy,OptimalValue,vars,prime_vars,orig_vars);
      
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      //save the value add and the policy to a file.
      if (dotFlag) {
	sprintf(basepath,"%s-APP",outputfile);
	writeVPDot(ApproximateValue,ApproximatePolicy,basepath);
      }
      sprintf(basepath,"%s-APPDual.ADD",outputfile);
      numwrite = writeDualOptimal(basepath,ApproximatePolicy,ApproximateValue,vars,prime_vars,orig_vars);
      break;
    }      
}


// start policy Querying server
void MDP::startServer() {
  if (mvinput) {
    switch (approxMeth) 
      {
      case APPROX_NONE:      
	policyServe2(gbm, OptimalPolicy, OptimalValue, numorigvars);
	break;
      case APPROX_ROUNDOFF:
      case APPROX_ALLPAIRS:
	policyServe2(gbm, ApproximatePolicy, ApproximateValue, numorigvars);

	break;
      }    
  }
}

int MDP::consultOptimalPolicy(DdNode *pol, DdNode *val, int *varvals)
{
  Pair dval, aval;
  getAV(gbm,pol,val, aval,dval, varvals, vars, numvars, orig_vars, numorigvars);
  int bestact;
  if ((bestact = pickAction((int) (aval.get_min()))) < 0) {
    fprintf(stderr,"error in action\n");
    exit(0);
  }
  return bestact;
  
}
int MDP::consultOptimalPolicy(int *varvals) {
  Pair dval, aval;
  getAV(gbm,OptimalPolicy, OptimalValue, aval,dval, varvals, vars, numvars, orig_vars, numorigvars);
  int bestact;
  if ((bestact = pickAction((int) (aval.get_min()))) < 0) {
    fprintf(stderr,"error in action\n");
    exit(0);
  }
  return bestact;
}  

void MDP::policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars) {
  Pair dval, aval;
  Client *spudd_c;
  Server *spudd_s;
  // connect to client and server
  if (!connectToSupervisor(&spudd_c,&spudd_s)) {
    fprintf(stderr,"Error - could not connect to supervisor\n");
    exit(0);
  }

  int *varvals = new int[numovars];
  int bestact;
  char namestr[1024];
  while (true) {
    if (spudd_c->isDataAvailableToRead()) {
      //fprintf(stderr,"got some\n");
      int size = numovars*sizeof(int);
      spudd_c->net_get_data (varvals, size);
      
      //fprintf(stderr,"received ");
      //for (int i=0; i<numovars; i++)
      //fprintf(stderr,"%d ",varvals[i]);

      getAV(gb,act,val,aval,dval,varvals,vars,numvars,orig_vars,numorigvars);


      // send back the action
      if ((bestact = pickAction((int) (aval.get_min()-1))) < 0) {
	fprintf(stderr,"error in action\n");
	exit(0);
      }
      strcpy(namestr,"");
      aconvert(namestr,actionnames,aval.get_min(),"/");
      //fprintf(stderr,"   sending %d which is %s\n",bestact,namestr);
      
      spudd_s->net_write_data(&bestact,sizeof(int));
    } else {
      usleep(1000);
    }
  }
}

void MDP::policyServe2(DdManager *gb, DdNode *act, DdNode *val, int numovars) {
  Pair dval, aval;
  Client *spudd_c;
  Server *spudd_s;

  int *varvals = new int[numovars];
  int bestact;
  char namestr[1024];
  while (true) {

    // connect to client and server
    if (!connectToSupervisor(&spudd_c,&spudd_s)) {
      fprintf(stderr,"Error - could not connect to supervisor\n");
      exit(0);
    }

    if (spudd_c->isDataAvailableToRead()) {
      fprintf(stderr,"got some\n");
      int size = numovars*sizeof(int);
      spudd_c->net_get_data (varvals, size);
      
      fprintf(stderr,"received ");
      for (int i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);

      getAV(gb,act,val,aval,dval,varvals,vars,numvars,orig_vars,numorigvars);


      // send back the action
      if ((bestact = pickAction((int) (aval.get_min()-1))) < 0) {
	fprintf(stderr,"error in action\n");
	exit(0);
      }
      strcpy(namestr,"");
      aconvert(namestr,actionnames,aval.get_min(),"/");
      fprintf(stderr,"   sending %d which is %s\n",bestact,namestr);
      
      spudd_s->net_write_data(&bestact,sizeof(int));
    } else {
      usleep(1000);
    }
    delete spudd_c;
    delete spudd_s;
  }
}
void MDP::reorder(DdNode *accordingTo) {
  /* do an initial reordering - minspan just uses the reward function */
  switch (reorderMeth)
    {
    case REORDER_SIFT:
    case REORDER_SP_SIFT:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_SIFT, 0))
	fprintf(stderr,"*****ERROR: Reorder sift failed \n");
      break;
    case REORDER_RANDOM:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_RANDOM, 0))
	fprintf(stderr,"*****ERROR: Reorder random failed \n");
      break;
    case REORDER_MINSPAN:
      if (!reorderMinSpan(accordingTo))
	fprintf(stderr,"******ERROR: Reorder minSpan failed \n");
      break;
    case REORDER_EXACT:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_EXACT, 0))
	fprintf(stderr,"*****ERROR: Reorder exact failed\n");
      break;
    case REORDER_NONE:
      break;
    }


}  

void MDP::allocateMemory() {
  int i, j, k;
  DdNode ***ActionsPrimed;

  /* WE ALLOCATE MEMORY OF THE DATA STRUCTURES USED IN THE PROGRAM*/
  //VECTORS............
  Array  = (DdNode **)malloc(numvars*(sizeof(DdNode *)));
  Array1 = (DdNode **)malloc(numvars*(sizeof(DdNode *)));   //for vars
  Array2 = (DdNode **)malloc(numvars*(sizeof(DdNode *)));   //for primes vars
  ArraySum = (DdNode **)malloc(numvars*(sizeof(DdNode *)));
  
  //WE initialize Array1 and Array2
  for(i=0;i<numvars;i++)
    {
      Array1[i] = vars[i].add_var;
      Cudd_Ref(Array1[i]);
      Array2[i] = prime_vars[i].add_var;
      Cudd_Ref(Array2[i]);
    }
  
  //MATRICES............
  Allprime = (DdNode ***)malloc(numactions*(sizeof(DdNode**)));
  for(i=0;i<numactions;i++)
    Allprime[i] = (DdNode **)malloc(numvars*(sizeof(DdNode*)));

  
  varnames = (char **) malloc(2*numvars*(sizeof(char *)));
  for(i=0;i<2*numvars; i++)
    varnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
     
  actionnames = (char **) malloc(numactions*(sizeof(char *)));
  for(i=0;i<numactions; i++)
    actionnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
  
  lnames = (char **) malloc((numvars+3)*sizeof(char *));
  for (i=0; i<(numvars+3); i++) 
    lnames[i] = (char *) malloc(MAXLINE*(sizeof(char)));
  

 
  // extra stuff to do if old style input
  if (!mvinput) {
    cout << "\Not mvinput!";
    ActionsPrimed = (DdNode ***)malloc(numactions*(sizeof(DdNode **)));
    for(i=0;i<numactions;i++)
      ActionsPrimed[i] = (DdNode **)malloc(numvars*(sizeof(DdNode*)));
    
    NewPrime = (DdNode ***)malloc(numactions*(sizeof(DdNode **)));
    for(i=0;i<numactions;i++)
      NewPrime[i] = (DdNode **)malloc(numvars*(sizeof(DdNode*)));
    // WE build the ADD for 1-Action tree
    for(i=0;i<numactions;i++)
      {
	for(j=0;j<numvars;j++)
	  {
	    ActionsPrimed[i][j] = Cudd_addApply(gbm,Cudd_addMinus,One,Actions[i][j]);
	    Cudd_Ref(ActionsPrimed[i][j]);
	  } 
      }
    
    // We build the new prime ADD, i.e., the ADD rooted at a primed variable with Actions to the right 
    // and ActionsPrimed to the left
    // We do that for each variable within each action
    
    for(i=0;i<numactions;i++)
      {
	for(j=0;j<numvars;j++)
	  {
	    NewPrime[i][j] = Cudd_addIte(gbm,prime_vars[j].add_var,Actions[i][j],ActionsPrimed[i][j]);
	    Cudd_Ref(NewPrime[i][j]);
	  }
      }

    for(i=0;i<numactions;i++) {
      for (j=0; j<numvars; j++) 
	Cudd_RecursiveDeref(gbm,ActionsPrimed[i][j]);
      free(ActionsPrimed[i]);
    }
    free(ActionsPrimed);
    
  }


}
void MDP::buildAllPrime() {
  // Now, build up the Allprime trees - We want to build a minimum number of allprime
  // trees for each action, such that each contains less than bigadd nodes
  // We multiplys all the NewPrimed ADD until they reach a upper limit size...then we build another ADD.
  // We repeat that for each action.
  int i,j,k;
  int flag;
  DdNode *temp1, *temp2;
  int dagCount;
  
  for (i=0; i<numactions;i++){ 
    
    j=0;
    k=0;
    //JH09/05/00list[i][k] = Cudd_NodeReadIndex(prime_vars[0].add_var)+1;
    allPlist[i][k] = Cudd_NodeReadIndex(prime_vars[0].add_var);
    k++;
    while (j<numorigvars) {
      flag = 1;
      temp1 = One; 
      Cudd_Ref(temp1); 

      while ((flag==1) && j < numorigvars) {
	
	temp2 = Cudd_addApply(gbm,Cudd_addTimes, NewPrime[i][j], temp1);
	Cudd_Ref(temp2);

	dagCount = Cudd_DagSize(temp2);
	
	if(Cudd_ReadMemoryInUse(gbm)> MAXMEMHARD) {
	  fprintf(stderr,"Memory required exceeds availabe %ld bytes. \n Try a smaller bigadd limit constant\n",MAXMEMHARD);
	  Cudd_Quit(gbm);
	  exit(0);
	}

	if (dagCount < bigadd) {                   /*we're still ok*/
	  Cudd_RecursiveDeref(gbm,temp1);
	  temp1 = temp2;
	  //Cudd_Ref(temp1);
	  //Cudd_RecursiveDeref(gbm,temp2);
	  j++; 
	  
	} else {
	  fprintf(stderr,"The %d allprime diagram got too big - %d nodes %f\n",k,dagCount,bigadd);
	  flag = 0;
	  /** NEW FEB 4 **/
	  Cudd_RecursiveDeref(gbm,temp2);
	}
      }
      
      fprintf(stderr,"Assigned complete action diagram %d for action %d with %d nodes \n",k-1,i,dagCount);
      
      Allprime[i][k-1] = temp1;
      // we want the index of the first binary variable which corresponds to the
      // jth  original variable, primed
      if (mvinput) {
	allPlist[i][k] = Cudd_NodeReadIndex(prime_vars[orig_vars[j-1].var1index].add_var)+2*orig_vars[j-1].nbvars;
      } else {
	allPlist[i][k] = Cudd_NodeReadIndex(prime_vars[j-1].add_var)+2;
      }
      k++;
    }
  }
  
}
void MDP::approximate() {

  int dagCount, dagLeafCount;
  double max_error,tmp_error,tol;
  // approximation stuff
  /* Nodes in the Value ADD*/
  dagCount = Cudd_DagSize(Merged_Add); // + Cudd_CountLeaves(Merged_Add);
  dagLeafCount = Cudd_CountLeaves(Merged_Add);
  tmp_error = get_error(Merged_Add);
    

  
  //fprintf(stderr,"nodes in value ADD after summation %d\n",dagCount);

  tol = -1.0;
  if (limit_type == LIMIT_ERROR && (counter > -1)) //We always prune 
    tol = ERROR_TOL;  //tmp_error > Max_Error) //&&	!close_enough(tmp_error,Max_Error,ERROR_TOL)) 
  else if (limit_type == LIMIT_SIZE && 
	   dagCount > Max_Size) // &&  !close_enough((double) dagCount, Max_Size, SIZE_TOL))
    tol = SIZE_TOL;
  //fprintf(stdout,"Dead nodes at end of loop - before approx: %d\n",Cudd_ReadDead(gbm));
  

  // approximate
  if (tol > 0) {

    /* Nodes in the Value ADD*/
    dagCount = Cudd_DagSize(Merged_Add); //+ Cudd_CountLeaves(Merged_Add);
    //fprintf(stderr,"nodes in value ADD after reorder %d\n",dagCount);
    
    // We need to change tolerance if it is in the sliding mode
    switch (toleranceMode) 
      {
      case TOLERANCE_SLIDING:
	Max_Error = prune*extR*((pow(discount_factor,((double)counter +1))-1)/(discount_factor-1));
	//printf("we change the tolerance\n");
	break;
      default:
	Max_Error = prune*extR;
	break;
      }
    
    /* Nodes in the Value ADD*/
    dagCount = Cudd_DagSize(Merged_Add); // + Cudd_CountLeaves(Merged_Add);
    max_error = get_error(Merged_Add);
    //fprintf(stderr,"nodes in value ADD before approx %d\n",dagCount);
    
    if(limit_type != LIMIT_SIZE || dagCount > Max_Size) {
      // Our two approximation techniques
      switch (approxMeth) 
	{
	case APPROX_ROUNDOFF:
	  size_approx();
	  break;
	case APPROX_ALLPAIRS:
	  allPairsApprox(&Merged_Add);
	  break;
	default:
	  break;
	}
    }
    dagCount = Cudd_DagSize(Merged_Add); //+ Cudd_CountLeaves(Merged_Add) ;
    //fprintf(stderr,"nodes in value ADD after approx %d\n",dagCount);
    
  }
  // end of approximate stuff
}

int MDP::Spudd(DdNode *Reward)
{
  int test;
  int i;
  double temps0,temps1;
  long int usememory=0;
  DdNode *actionDD, *MidPoint_Add;
  DdNode *temp1 = NULL, *temp2 = NULL;
  DdNode **Vcurrent,*Vpast;
  //DdNode *Reward;

  /* Structure to get the time of execution*/
  struct tms _t;
  long clk_tck = sysconf(_SC_CLK_TCK);
  
  //FILE *value,*action,*both,*stats,*log,*log2,*reward,*actions;

  if (Reward == NULL) {
    Reward = RewardD;
  }
  Cudd_Ref(Reward);

  Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));

  for(i=0;i<numactions;i++) {
    Vcurrent[i] = RewardD;
    Cudd_Ref(Vcurrent[i]);
  }  

  // Convergence test
  test = 0;
  //Number of iteration to convergence
  counter = 0;
  
  Vpast = RewardD;
  Cudd_Ref(Vpast);
  
  Merged_Add = Reward;
  Cudd_Ref(Merged_Add);
    
  
  actionDD = Zero;	/* Stupid ref to maintain loop invariant */
  Cudd_Ref(actionDD);

  MidPoint_Add = VerySmall;
  Cudd_Ref(MidPoint_Add);
  
  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  
  int lastiteration=0;
  if (horizon == 1)
    lastiteration = 1;
  /* We get the starting up time */
  times(&_t);
  temps0 = 1.0*_t.tms_utime/clk_tck;

  // swap action costs to primed variables
#ifdef COSTSPRIMED
  for(i=0;i<numactions;i++) { 
    temp1 = Cudd_addSwapVariables(gbm,actionCost[i],Array1,Array2,numvars); 
    Cudd_Ref(temp1); 
    Cudd_RecursiveDeref(gbm,actionCost[i]);
    actionCost[i]= temp1;
  }
#endif
  /* VALUE ITERATION LOOP COMMENCES */
  while ((horizon == -1 && (test==0 || lastiteration == 1)) || (horizon >= 0 && counter < horizon)){ 

    /* KEEP TRACK OF THE NUMBER OF ITERATIONS */
    counter = counter + 1;
    fprintf(stderr,"\niteration %d\n",counter);
    
    /* update past Value function for convergence check */
    Cudd_RecursiveDeref(gbm,Vpast); 
    Vpast = Merged_Add; 
    Cudd_Ref(Vpast); 

    // print out stats
    //printStats(One,Merged_Add,stderr);
    /* Reward is Merged_Add with variables primed */
    Cudd_RecursiveDeref(gbm,Reward); 
    Reward = Cudd_addSwapVariables(gbm,Merged_Add,Array1,Array2,numvars); 
    Cudd_Ref(Reward); 

    Cudd_RecursiveDeref(gbm,actionDD); 
    actionDD = Zero;	
    Cudd_Ref(actionDD);

    /* Merged_Add maintains the current Value function*/
    Cudd_RecursiveDeref(gbm,Merged_Add);
    Merged_Add = VerySmall;
    Cudd_Ref(Merged_Add);
 
    /*Loop over the states */
    for(i=0;i<numactions;i++) { 
      
      /* do the sum_t(Pr(s,a,t)*V(t)) */
#ifdef OPTIMIZED
      temp1 = primeReward(Reward,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
      temp1 = multiplySumSet(Reward,NewPrime[i]);
#endif
      //fprintf(stderr,"iteration %d action %d\n",counter,i);
      //Cudd_PrintDebug(gbm,temp1,4,100);
      /* multiply by the discount factor */
      temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
      Cudd_Ref(temp2);

      /* add the reward to get the Vcurrent for this action*/
      Cudd_RecursiveDeref(gbm,temp1);
      temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
      Cudd_Ref(temp1);
      
      Cudd_RecursiveDeref(gbm,temp2);
      
      // add the action costs. 
      temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp1,actionCost[i]);
      Cudd_Ref(temp2);

      Cudd_RecursiveDeref(gbm,temp1);
      

      /* maximize over actions */
      if (approxMeth == APPROX_NONE) {
	temp1=Cudd_addApply(gbm,Cudd_addMaximum,Merged_Add,temp2);
      } else {
	temp1=Cudd_addApply(gbm,My_addMaximum,Merged_Add,temp2);
      }
      Cudd_Ref(temp1);
      
#ifdef DEBUGP
      FILE *ffileh = fopen("shit.dot","w");
      DumpDot_p(gbm, temp1, vars, prime_vars, numvars, orig_vars, numorigvars, NULL, ffileh);
      fclose(ffileh);
#endif
      
      // get the policy during the final iteration
      if (lastiteration) {
	if (approxMeth == APPROX_NONE) {
	  // add this action to actionDD at each node where temp2 == temp1
	  // replace it when temp2 > temp1
	  //fprintf(stdout,"calling recover policy at %d\n",i);
	  recoverPolicy(temp2,temp1,Merged_Add,&actionDD,i);
	  // save the Qfunction
	  Qf[i] = temp2;
	  Cudd_Ref(temp2);
	  Cudd_RecursiveDeref(gbm,temp2);
	} else {
	  // get MidPoint ADD
	  // temp3 is the midpoint Q add for this action
	  DdNode *temp3 = Cudd_addApply(gbm, addMean, temp2, One);
	  Cudd_Ref(temp3);
	  Cudd_RecursiveDeref(gbm,temp2);
	  
	  // get new MidPoint_Add as those maxima
	  temp2=Cudd_addApply(gbm,Cudd_addMaximum,MidPoint_Add,temp3);
	  Cudd_Ref(temp2);

	  // add this action to actionDD policy if  temp3 > temp2 (current MidPoint_Add) 
	  recoverPolicy(temp3,temp2,MidPoint_Add,&actionDD,i);
	  Cudd_RecursiveDeref(gbm,temp3);
	  Cudd_RecursiveDeref(gbm,MidPoint_Add);
	  MidPoint_Add = temp2;
	}
      } else {
	Cudd_RecursiveDeref(gbm,temp2);
      }

      Cudd_RecursiveDeref(gbm,Merged_Add);
      Merged_Add = temp1;
      

      usememory = Cudd_ReadMemoryInUse(gbm);
      if (maxusememory < usememory)
	maxusememory = usememory;
      if (usememory > MAXMEMHARD) {
	fprintf(stderr,"Memory required exceeds availabe %d bytes. \nTry a smaller bigadd limit constant\n",MAXMEMHARD);
	Cudd_Quit(gbm);
	return 0;
      }

    } // end of loop over actions

    //test = Cudd_DagSize(Merged_Add);
    //fprintf(stderr,"size of function is now %d\n",test);
    //Cudd_PrintInfo(gbm,stderr);
      
    
    //fprintf(stderr,"memory in use is %ld\n",usememory);

      // attempt a reorder?
    if (reorderApplyMeth ==REORDERAPPLY_ALL || (reorderApplyMeth == REORDERAPPLY_FIRSTFIVE && counter < 5)) 
      reorder(Merged_Add);
    
    // approximation takes place here
    // could modify this so approximate takes Merged_Add as an argument - then 
    // Merged_Add would not have to be a member variable
    if (!lastiteration) 
      approximate();

    /* COMPARE THE TWO MOST RECENT ADDs IN ORDER TO VERIFY THE CONVERGENCE  */
    if (approxMeth == APPROX_NONE) {
      test = Cudd_EqualSupNorm(gbm,Merged_Add,Vpast,toleranceP,0);
    } else {
      DdNode *testD = Cudd_addApply(gbm,Convergence_Test,Vpast,Merged_Add);
      Cudd_Ref(testD);
      if((*Cudd_V(testD)).get_min() == 1.0)
	test = 1;
      Cudd_RecursiveDeref(gbm,testD);
    }
    if (test>0 || lastiteration ||  (horizon >= 0 && counter == horizon-1))
      lastiteration++;
  }   /* END OF VALUE ITERATION */
  
  /* We get the time of completion*/
  times(&_t);
  temps1 = 1.0*_t.tms_utime/clk_tck;
  
  temps = temps1- temps0;
  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  
  // possibly deref the old optimal policy
  if (approxMeth == APPROX_NONE) {
    if (OptimalPolicy) 
      Cudd_RecursiveDeref(gbm,OptimalPolicy);
    if (OptimalValue) 
      Cudd_RecursiveDeref(gbm,OptimalValue);

    // set global OptimalPolicy and OptimalValue to new values
    OptimalPolicy = actionDD;
    Cudd_Ref(OptimalPolicy);
    OptimalValue = Merged_Add;
    Cudd_Ref(OptimalValue);
  } else {
    if (ApproximateMidPointValue)
      Cudd_RecursiveDeref(gbm,ApproximateMidPointValue);
    ApproximateMidPointValue = MidPoint_Add;
    Cudd_Ref(ApproximateMidPointValue);

    if (ApproximatePolicy)
      Cudd_RecursiveDeref(gbm,ApproximatePolicy);
    ApproximatePolicy = actionDD;
    Cudd_Ref(ApproximatePolicy);

    if (ApproximateValue)
      Cudd_RecursiveDeref(gbm,ApproximateValue);
    ApproximateValue =  Merged_Add;
    Cudd_Ref(ApproximateValue);
  }  

  //i added the lines until the return --jlb 7/20/10
  for(i=0;i<numactions;i++) {
    Cudd_RecursiveDeref(gbm, Vcurrent[i]);
  }  
  free(Vcurrent);
  Cudd_RecursiveDeref(gbm, Vpast);
  Cudd_RecursiveDeref(gbm, Reward);
  Cudd_RecursiveDeref(gbm, actionDD);
  Cudd_RecursiveDeref(gbm, MidPoint_Add);
  Cudd_RecursiveDeref(gbm, Merged_Add);
  //success
  return 1;
}
void MDP::printOptimalStats(char *statsfilename) {
  FILE *stats;
  stats = fopen(statsfilename,"w");
  printStats(OptimalPolicy,OptimalValue,stats);
  fclose(stats);
}
void MDP::printApproximateStats(char *statsfilename) {
  FILE *stats;
  stats = fopen(statsfilename,"w");
  printStats(ApproximatePolicy,ApproximateValue,stats);
  fclose(stats);
}
void MDP::printStats(DdNode *actionDD, DdNode *valueDD, FILE *stats) {
  
  int  dagCount,treeCount,treeLeafCount,dagLeafCount;
  
  fprintf(stats,"\n\n-------------------- Spudd stats  --------------------\n\n");
  fprintf(stats,"Discount factor is : %f \n Tolerance is: %f \n horizon is: %f \n The BIGADD limit is set to: %lf \n Target maximum memory use: %d bytes \n Hard Limit on memory use %d bytes\n",discount_factor, tolerance, horizon, bigadd,MAXMEM, MAXMEMHARD);
  
  fprintf(stats," Iterations to convergence %d\n",counter);
  fprintf(stats," Final execution time: %8.4f  seconds\n",temps);
  fprintf(stats," Maximum memory usage: %ld bytes\n",maxusememory);

  // nodes in the Action ADD 
  dagCount = Cudd_DagSize(actionDD);
  dagLeafCount = Cudd_CountLeaves(actionDD);

  //Nodes in the equivalent action tree 
  treeCount = count_internal_nodes_tree(actionDD);
  treeLeafCount = count_leaves_tree(actionDD);

  fprintf(stats,"Number of nodes in the action ADD: %d  internal nodes   %d leaves %d total nodes\n",
	  dagCount-dagLeafCount,dagLeafCount,dagCount);
  fprintf(stats,"Number of nodes in the equivalent tree: %d internal nodes  %d leaves  %d total nodes\n",
	  treeCount,treeLeafCount,treeCount+treeLeafCount); 
  
  // Nodes in the Value ADD
  dagCount = Cudd_DagSize(valueDD);
  dagLeafCount = Cudd_CountLeaves(valueDD);

  //Nodes in the equivalent value tree 
  treeCount = count_internal_nodes_tree(valueDD);
  treeLeafCount = count_leaves_tree(valueDD);

  fprintf(stats,"Number of nodes in the value ADD: %d  internal nodes   %d leaves %d total nodes\n",
	   dagCount-dagLeafCount,dagLeafCount,dagCount);   
  fprintf(stats,"Number of nodes in the value tree: %d internal nodes  %d leaves  %d total nodes\n",
	  treeCount,treeLeafCount,treeCount+treeLeafCount);
  
  fprintf(stats,"\n\n-------------------- other info from dd Manager --------------------\n\n");
  Cudd_PrintInfo(gbm,stats);
}

// currently dysfunctional
void MDP::printAdd(DdNode *add, FILE *fp)
{
  char tabstop[1024];
  sprintf(tabstop,"");
  //printDdNode(gbm,add,vars,prime_vars,numvars,orig_vars,numorigvars,fp,tabstop);
}

// multiply reward by each primes[i] (using multiplySum)
// and sum out all primed variables
// return is pre-reffed
DdNode *MDP::multiplySumSet(DdNode *reward, DdNode **primes) {
  int j;
  DdNode *temp1, *temp;
  temp1 = reward;
  Cudd_Ref(temp1);
  for (j=0; j<numorigvars; j++) {
    temp = multiplySum(temp1,primes[j],prime_vars[j].add_var,orig_vars+j);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;

#ifdef DEBUGP
    FILE *fileh = fopen("shit.dot","w");
    DumpDot_p(gbm, temp1, vars, prime_vars, numvars, orig_vars, numorigvars, NULL, fileh);
    fclose(fileh);
#endif

  }
  // it is still possible that there are some primed variables
  // since they can be in the CPTs - sum them out!
  // no this does not work because you have to sum out a parent
  // before a child node
  /*
  if (temp = sumPrimes(temp1,prime_vars,numvars)) {
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;
  }
  */
  return temp1;
}

// multiply reward by prime and sum out varsum
// return is pre-reffed
DdNode *MDP::multiplySum(DdNode *reward, DdNode *prime, DdNode *varsum, onum *ovar) {
  
  DdNode *temp, *temp1;

  temp1 = Cudd_addApply(gbm,Cudd_addTimes,prime,reward);
  Cudd_Ref(temp1);

  if (mvinput) {
    temp = sumOutPrime(temp1,ovar,prime_vars);
  } else {
    DdNode *cube = Cudd_addIte(gbm,varsum,One,Zero);
    Cudd_Ref(cube);
    
    temp = Cudd_addExistAbstract(gbm,temp1,cube);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,cube);
  }
  // JH 26/05/04 added this...
  Cudd_RecursiveDeref(gbm,temp1);
  return temp;
}
// old way
// mutliply reward by each dual action diagram and sum out subtrees
DdNode *MDP::newPrimeReward(DdNode *reward, DdNode **NewPrime) {
  
  int j;
  DdNode *temp, *temp1;

  temp1 = reward;
  Cudd_Ref(temp1);
  for (j=0; j<numorigvars; j++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,NewPrime[j],temp1);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;

    if (mvinput) {
      temp = sumOutPrime(temp1,orig_vars+j,prime_vars);
    } else {
      DdNode *cube = Cudd_addIte(gbm,prime_vars[j].add_var,One,Zero);
      Cudd_Ref(cube);
      
      temp = Cudd_addExistAbstract(gbm,temp1,cube);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,cube);
    }
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;
  }
  return temp1;
}

/*
This does the multiplication of Pr(s,a,t) by V(t) summed over t
in the value iteration loop. Pr(s,a,t) is in AllPrime as an array of
the chunks of the complete action diagram. 
It does so by descending V(t) over the 'chunked'
levels to the lowest one, where it multiplies by the corresponding
chunk in Allprime, followed by a sum over the primed variables in 
the chunk. See section 4.2 in the Tech Report 
*/
DdNode  *MDP::primeReward(DdNode *reward, DdNode **AllPrime, rnum *prime_vars, 
                     int n, int *aplist, int numvars, int act)
{
  int i,test;
  DdNode *temp,*temp1,*tempT, *tempE;
  test = 0;

  /* Get the index of the root node of the reward tree */
  if ((Cudd_IsConstant(reward) == 1)) {
    temp = reward;
    Cudd_Ref(temp);
  } else {
    i = Cudd_NodeReadIndex(reward);
    
    //JH09/05/00if (i < aplist[n]) {
    // if we're at the end - do the multiplication
    if (i >= aplist[n]) {
      //JH09/05/00      if (aplist[n+1] == numvars) {
      if (aplist[n+1] == 2*numvars) {
        temp = reward;
        Cudd_Ref(temp);
      } else{
        temp = primeReward(reward, AllPrime, prime_vars, n+1, aplist, numvars,act);
      }
     
      temp1 = Cudd_addApply(gbm,Cudd_addTimes,AllPrime[n],temp);
      Cudd_Ref(temp1);
      
      Cudd_RecursiveDeref(gbm,temp);

      temp = sumSubtrees(temp1,aplist[n],aplist[n+1],prime_vars,numvars);
      Cudd_RecursiveDeref(gbm,temp1);
      
    } else {
      // descend to the next level
      tempT = primeReward(Cudd_Then(reward), AllPrime, prime_vars, n, aplist, numvars,act); 
      tempE = primeReward(Cudd_Else(reward), AllPrime, prime_vars, n, aplist, numvars,act); 
      //     temp = Cudd_addIte(gbm, prime_vars[2*numvars-1-i].add_var, tempT, tempE);
      // recombine the results
      temp = Cudd_addIte(gbm,prime_vars[i/2].add_var,tempT,tempE);
      Cudd_Ref(temp);
      /* deref these because I think that Cudd_addIte does its own reffing! */
      Cudd_RecursiveDeref(gbm, tempT);
      Cudd_RecursiveDeref(gbm, tempE);
    }
  }
  return temp;
}

// Vmerged is max(Vcurrent,Vprior)
// *actionAdd is the current policy
// recoverPolicy assigns 2^off to each state for which Vcurrent > Vprior
//                       actionDD+2^off to each state for which Vcurrent == Vprior
//                       actionDD       to each state for which Vcurrent < Vprior
// 1) get ADD X_1 which is 1 wherever Vcurrent > Vprior, 0 elsewhere
// 2) X_2 = (1-X_1)*actionDD   has old actions everywhere where Vcurrent <= Vprior
// 3) get ADD X_3 which is 1 whereever Vcurrent >= Vprior, 0 elsewhere
// 4) X_4 = X_3*2^off
// 5) new actionDD is X_4 + X_2
// This function is used to add to the policy *actionAdd the off action which contributed 
// 
void MDP::recoverPolicy(DdNode *Vcurrent, DdNode *Vmerged, DdNode *Vprior, DdNode **actionAdd, double off)
{

  DdNode *temp,*temp2,*temp3, *X_4, *X_2, *offset;
  Pair offPair;
  double poff;

  // get X_4
  // subtract the Vcurrent from Vmerged - the nodes which match are 0 in the result
  temp = Cudd_addApply(gbm,Cudd_addMinus,Vcurrent,Vmerged);
  Cudd_Ref(temp);

  
  //turn the resulting add into a BDD by thresholding at 0 - now the 1 nodes of the BDD are the 
  // states for which the kth action was the one that contributed to Vmerged - that is
  // the states for which Vcurrent >= Vprior
  temp2 = Cudd_addBddThreshold(gbm,temp,pZero);
  Cudd_Ref(temp2);

   // turn that BDD back into a 0-1 ADD
  Cudd_RecursiveDeref(gbm,temp);
  temp = Cudd_BddToAdd(gbm,temp2);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,temp2);

  // temp is X_3: a 0-1 ADD with all 1's are states where the
  //kth action generated the value in Vmerged
  //  offset the value by 1e+off (this is to keep track of ALL the possible actions)
  
#ifdef ALLACTIONS
  poff = pow(2.0,off);
  offPair.set(poff);
  offset = Cudd_addConst(gbm,&offPair);  
  Cudd_Ref(offset);
#else
  offPair.set(off+1);
  offset = Cudd_addConst(gbm,&offPair);
  Cudd_Ref(offset);
#endif

  //X_4 is the same structure as temp but has the
  //number 2e+k instead of the 1 nodes
  X_4 = Cudd_addApply(gbm,Cudd_addTimes,temp,offset);
  Cudd_Ref(X_4);
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,offset);

  //temp = Cudd_addApply(gbm,Cudd_addPlus,temp2,*actionAdd);
  //  temp = Cudd_addApply(gbm,MyPlus,temp2,*actionAdd);
  //Cudd_Ref(temp);
  //Cudd_RecursiveDeref(gbm,temp2);
    
  //fprintf(stdout,"X_4 is\n");
  //Cudd_PrintDebug(gbm,X_4,4,100);

  // now, get X_2
  // first get X_1 - an ADD which is 1 wherever Vcurrent > Vprior
  // 
  temp2 = Cudd_addApply(gbm,Cudd_addMinus,Vcurrent,Vprior);
  Cudd_Ref(temp2);

  // turn into a BDD with 1s wherever temp2 is different from 0 (0s wherever Vcurrent and Vprior agree)
  // must use addBddPattern here because Vcurrent may be less than *or* greater than Vprior (so addBddThreshold will not work)
  //temp3 = Cudd_addBddPattern(gbm,temp2);

  // turn into a BDD with 1s wherever temp2 is strictly greater than 0 (where Vcurrent > Vprior
#ifdef ALLACTIONS
  temp3 = Cudd_addBddStrictThreshold(gbm,temp2,pZero);
#else
  temp3 = Cudd_addBddThreshold(gbm,temp2,pZero);
#endif
  Cudd_Ref(temp3);
  Cudd_RecursiveDeref(gbm,temp2);
  
  // turn back to a ADD
  temp2 = Cudd_BddToAdd(gbm,temp3);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp3);
  
  // we actually want the negation of this (1s wherever Vcurrent <= Vprior)
  temp3 = Cudd_addApply(gbm,Cudd_addMinus,One,temp2);
  Cudd_Ref(temp3);
  Cudd_RecursiveDeref(gbm,temp2);
  
  // get the actions in actionAdd at each state where Vcurrent and Vprior agree
  X_2 = Cudd_addApply(gbm,Cudd_addTimes,temp3,*actionAdd);
  Cudd_Ref(X_2);
  Cudd_RecursiveDeref(gbm,temp3);
  Cudd_RecursiveDeref(gbm,*actionAdd);
  
  //fprintf(stdout,"X_2 is\n");
  //  Cudd_PrintDebug(gbm,X_2,4,100);
  // X_2  has actionAdd at each leaf where Vprior was the same as Vcurrent
  // so both actions need to be remembered
  temp3 = Cudd_addApply(gbm,Cudd_addPlus,X_2,X_4);
  Cudd_Ref(temp3);
  Cudd_RecursiveDeref(gbm,X_2);
  Cudd_RecursiveDeref(gbm,X_4);
  *actionAdd = temp3;
  //*actionAdd = temp;
}
int  MDP::close_enough(double val, double limit, double tol) 
{
  if(fabs(limit-val) < tol)
    return 1;
  else 
    return 0;
}



DdNode * MDP::roundOffMMA(DdNode *res)
{
  DdGen *gen;
  DdNode *node, *tmp, *tmp2, *dRes, *thisLeaf;
  Pair tmpError;
  
  
  //tmpError = new Pair();

  dRes = Zero;
  Cudd_Ref(dRes);

  Cudd_ForeachNode(gbm,res,gen,node) {
    if (Cudd_IsConstant(node)) {
      avgError.set_max(-1.0*LARGEFLAG);
      avgError.set_min(1.0*LARGEFLAG);
      
      
      tmp = Cudd_addApply(gbm, addAddExact, res, node);
      Cudd_Ref(tmp);
     
      //tmp2 = Cudd_addApply(gbm, getTotalSpan, tmp, Merged_Add);
      ///NEED TO MODIFY THIS BACK TO MYADDAPPLY!!!!!!!
      tmp2 = My_addApply(getTotalSpan, tmp, Merged_Add);
      Cudd_Ref(tmp2);
      
      tmpError = avgError;
      thisLeaf = Cudd_addConst(gbm,&tmpError);
      Cudd_Ref(thisLeaf);
      
      Cudd_RecursiveDeref(gbm,tmp2);
      tmp2 = Cudd_addApply(gbm, Cudd_addTimes, thisLeaf, tmp);
      Cudd_Ref(tmp2);
      
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,thisLeaf);
       
      tmp = Cudd_addApply(gbm, Cudd_addPlus, dRes, tmp2);
      Cudd_Ref(tmp);
      Cudd_RecursiveDeref(gbm,tmp2);
      Cudd_RecursiveDeref(gbm,dRes);
      dRes = tmp;
    }
  }
  return (dRes);
}

void  MDP::size_approx() 
{
  DdNode *res, *tmp;
  DdNode *dRoe;
  double oldroe = 0.0, olderoe;
  
  int count, oldcount;
  int iter1(0),iter2(0),maxiters(100);
  Pair pRoe;
    
  //just to make sure that it is bigger
  double roe = get_extent(Merged_Add)*2.0;

  if (limit_type == LIMIT_ERROR) {
    pRoe.set(Max_Error);
    fprintf(stderr,"Limit Error - Max_Error: %f \n",Max_Error);
  } else if (limit_type == LIMIT_SIZE) {
    pRoe.set(roe);
    // double precision digit limit
    // this should be a little bigger than unique->epsilon in the dd manager
    Max_Error = 2.0*(Cudd_ReadEpsilon(gbm)->get_min());
    //fprintf(stderr,"Limit Size - Max_Size: %d\n",Max_Size);
  }

  dRoe = Cudd_addConst(gbm,&pRoe);
  Cudd_Ref(dRoe);

  cntError = 0;

  res = Cudd_addApply(gbm, myRoundOff,Merged_Add, dRoe);
  Cudd_Ref(res);
  count = Cudd_DagSize(res);
  
  oldcount = -1;
  if (limit_type == LIMIT_ERROR) {
    /* replace Merged_Add with rounded off version */
    tmp = roundOffMMA(res);
    Cudd_RecursiveDeref(gbm,Merged_Add);
    Merged_Add = tmp;
    Cudd_RecursiveDeref(gbm,res);
  } else if (limit_type == LIMIT_SIZE) {

    /* while ADD is too small       - decrease the roe - increase the ADD size */
    // the maxiters thing is only to avoid infite loops which used to happen - 
    // probably should not be there
    oldroe = 0.0;
    while (iter1 < maxiters && count < Max_Size && (fabs(oldroe-roe) > Max_Error)) {
      oldcount = count;
      //fprintf(stderr,"count %d Max_Size %d  error %e max_error %e oldroe %e roe %e ",
      //      count,Max_Size,fabs(oldroe-roe),Max_Error,oldroe,roe);
      roe = roe - (roe - oldroe)/2.0;
      //fprintf(stderr,"new roe %e\n",roe);
      pRoe.set(roe);
      
      Cudd_RecursiveDeref(gbm,dRoe);
      dRoe = Cudd_addConst(gbm,&pRoe);
      Cudd_Ref(dRoe);
      
      cntError = 0;
      avgError.set(0.0);
      
      Cudd_RecursiveDeref(gbm,res);
      res = Cudd_addApply(gbm, myRoundOff,Merged_Add, dRoe);
      Cudd_Ref(res);
      
      count = Cudd_DagSize(res);
      if (count > Max_Size) 
	oldroe = roe*2.0 - oldroe;
      
      /* while the ADD is too big - increase the roe - decrease the ADD size */
      iter2=0;
      while (iter2 < maxiters && count > Max_Size & (fabs(oldroe-roe) > Max_Error)) {
	/* replace Merged_Add with rounded off version 
	   since we know it will be smaller than this */
	//fprintf(stderr,"--------count %d Max_Size %d  error %e max_error %e oldroe %e roe %e ",
	//	count,Max_Size,(oldroe-roe)/2.0,Max_Error,oldroe,roe);
	
	roe = roe + (oldroe - roe)/2.0;
	pRoe.set(roe);
	//fprintf(stderr,"new roe %f\n",roe);
	
	Cudd_RecursiveDeref(gbm,dRoe);
	dRoe = Cudd_addConst(gbm,&pRoe);
	Cudd_Ref(dRoe);


	Cudd_RecursiveDeref(gbm,res);
	res = Cudd_addApply(gbm, myRoundOff,Merged_Add, dRoe);
	Cudd_Ref(res);
 
	
	count = Cudd_DagSize(res);
	if (count <= Max_Size) {
	  olderoe = roe;
	  oldroe = roe*2.0 - oldroe;
	  //fprintf(stderr,"**** resetting count %d Max_Size %d  oldroe %e roe %e\n ",
	  //  count,Max_Size,oldroe,roe);	
	}
	else if (fabs(roe-oldroe) < Max_Error) {
	  roe = olderoe;
	  while (count > Max_Size) {
	    //fprintf(stderr,"?????????? count %d Max_Size %d  error %e max_error %e oldroe %f roe %f ",
	    //	    count,Max_Size,(oldroe-roe)/2.0,Max_Error,oldroe,roe);
	    pRoe.set(roe);
	    
	    Cudd_RecursiveDeref(gbm,dRoe);
	    dRoe = Cudd_addConst(gbm,&pRoe);
	    Cudd_Ref(dRoe);
	    
	    Cudd_RecursiveDeref(gbm,res);
	    res = Cudd_addApply(gbm, myRoundOff,Merged_Add, dRoe);
	    Cudd_Ref(res);
	    count = Cudd_DagSize(res);
	    roe += tolerance;
	  }
	  oldroe = roe;
	}
	iter2++;
	iter1=0;
      }
    }
    //fprintf(stderr,"roundoffMMA\n");
    tmp = roundOffMMA(res);
    Cudd_RecursiveDeref(gbm,Merged_Add);
    Merged_Add = tmp;
    Cudd_RecursiveDeref(gbm,res);
  }
}
int  MDP::reorderMinSpan(DdNode *addBase)
{
  double *infor;
  double s,st,se;
  int i,nt,ne,res;
  int *phase;
  int *reordList;
  DdNode *vRest, *oCube;
  DdNode *cVars;

  infor = (double *) malloc(numvars*sizeof(double));
  reordList = (int *) malloc(2*numvars*sizeof(int));
  phase = (int *) malloc(sizeof(int));

  for (i=0; i<2*numvars;i++)
    reordList[i] = i;
  
  s = get_extent(addBase);

  for (i=0;i<numvars;i++) {

    cVars = vars[i].add_var;
    Cudd_Ref(cVars);

    /* set this variable to 1 */
    phase[0] = 1;
    oCube = Cudd_addComputeCube(gbm,&cVars,phase,1);
    Cudd_Ref(oCube);


    vRest = Cudd_addRestrict(gbm,addBase,oCube);
    Cudd_Ref(vRest);
    Cudd_RecursiveDeref(gbm,oCube);

    st = get_extent(vRest);
    nt = Cudd_CountLeaves(vRest);
    Cudd_RecursiveDeref(gbm,vRest);

    /* set it to 0 */
    phase[0] = 0;
    oCube = Cudd_addComputeCube(gbm,&cVars,phase,1);
    Cudd_Ref(oCube);
    Cudd_RecursiveDeref(gbm,cVars);

    vRest = Cudd_addRestrict(gbm, addBase,oCube);
    Cudd_Ref(vRest);
    Cudd_RecursiveDeref(gbm,oCube);
    se = get_extent(vRest);
    ne = Cudd_CountLeaves(vRest);
    Cudd_RecursiveDeref(gbm,vRest);
    
    //infor[i] = (2*s-st-se)*nt/(s*(nt+ne));
    //infor[i] = -1.0*(log(1-st/s)+log(1-se/s))*nt/(nt+ne);
    infor[i] = 0;
    if (st > 0) 
      infor[i] += (nt/nt+ne)*log(st/s);
    else
      infor[i] -= 10000.0;
    if (se > 0)
      infor[i] += (ne/nt+ne)*log(se/s);
    else
      infor[i] -=  10000.0;

  }
  /* sort infor list and build reorder list (indices of variables which had
     smallest errors */

  infoSorter(infor,reordList,numvars);
  res = Cudd_ShuffleHeap(gbm,reordList);
  free(infor);
  free(reordList);
  free(phase);
  return(res);
  
}



//This approximation method finds the 2 pairs with minimal span and approximates them
// by the global span of the two
void  MDP::allPairsApprox(DdNode **theADD){
  DdGen *gen;
  DdNode *node,*Approx,*temp;
  list_cst *const_list,*ins,*curr,*suivant;
  double min,max,min_span,min_temp;
  double current_min;
  list_cst *elem1, *elem2;
  int count;
  int typeFlag;

  //const_list  is a pointer to a constant node list 
  //At first the list is Null
  const_list = NULL;
  
  //We need to find all the constant nodes so that we include them in the constant node list
  Cudd_ForeachNode(gbm,*theADD,gen,node) {
    if (Cudd_IsConstant(node)) { // I want to put the pointer in a list by order of minimum
      
      //we allocate the memory for the new constant node in the list
      ins = (list_cst *)malloc(sizeof(list_cst));
      ins->next = NULL;
      ins->add = node;
 
      Cudd_Ref(ins->add);
     
      //we now insert in the list
      //we begin at the start of the list
      curr = const_list;
      suivant = const_list;
      
      //if the list is empty, we put this element at the beginning
      if (curr == NULL)
	{
	  const_list = ins;
	}
      //we put this element at his good position:(increasing order of the min)
      else{
	suivant = curr->next;
	//if the new constant has a smaller min than the current one, put it at the start of the list.
	if( (*Cudd_V(curr->add)).get_min() >= (*Cudd_V(ins->add)).get_min()){
	  ins->next = curr;
	  const_list = ins;
	}else{
	  //we put the new constant at the good position in the list
	  while(suivant !=NULL){
	    if((*Cudd_V(suivant->add)).get_min() < (*Cudd_V(ins->add)).get_min())
	      {
		curr = suivant;
		suivant = suivant->next;
	      }else{
		break;
	      }
	  }
	  ins->next = suivant;
	  curr->next = ins;
	}
      }
    }
  }
  //we are done ordering the constant nodes by increasing size of their min.
  //now we need to find which pairs to approximate.
  if(const_list == NULL)
    printf("There is a problem with my list.....\n");
  
    
  
  nodeOne = Zero;
  Cudd_Ref(nodeOne);
  nodeTwo = Zero;
  Cudd_Ref(nodeTwo);
  //Need to optimize it so that I don't run through the whole list n^2 time......
  //We run through all the constant node to find the minimum span.
  /*############################

    MIGHT WANT TO KEEP MORE THAN ONE MINIMUM.....TO LOOK INTO MAYBE

    #############################

  */

  min_span = 0.0;
  count = Cudd_DagSize(*theADD);
  if(limit_type == LIMIT_ERROR)
    typeFlag = 1;
  else if(limit_type == LIMIT_SIZE)
    typeFlag = 0;

  while ( ((count > Max_Size) && (typeFlag == 0)) || (((min_span/2.0) <= Max_Error) && typeFlag == 1)) {
    
    //Set the span_min to a big number so that I find a new best span
    min_span = DBL_MAX;
    curr = const_list;  
    elem1 = NULL;
    elem2 = NULL;
    while(curr != NULL){
      current_min = (*Cudd_V(curr->add)).get_min(); 
      suivant = curr->next;
      
      while(suivant != NULL){
	if( ((typeFlag == 1) && ((2*Max_Error) > (*Cudd_V(suivant->add)).get_min() - current_min)) || ((typeFlag == 0) && (min_span > (*Cudd_V(suivant->add)).get_min() - current_min)))  {
	  if(min_span > (min_temp = length_span(curr->add,suivant->add)))
	    {
	      min_span = min_temp;
	      elem1 = curr;
	      elem2 = suivant;
	    }
	  suivant = suivant->next;
	}else{
	  break;
	}
      }
      curr = curr->next;
    }
    if(((min_span/2.0) > Max_Error) && (typeFlag ==1))
      break;
    
    Cudd_RecursiveDeref(gbm,nodeOne);
    nodeOne = elem1->add;
    Cudd_Ref(nodeOne);
    Cudd_RecursiveDeref(gbm,nodeTwo);
    nodeTwo = elem2->add;
    Cudd_Ref(nodeTwo);
    
    //We are done, now lets create the new constant that we will replace the two nodes with
    if((*Cudd_V(nodeOne)).get_min() <= (*Cudd_V(nodeTwo)).get_min())
      min = (*Cudd_V(nodeOne)).get_min(); 
    else
      min = (*Cudd_V(nodeTwo)).get_min(); 
    
    if((*Cudd_V(nodeOne)).get_max() <= (*Cudd_V(nodeTwo)).get_max())
      max = (*Cudd_V(nodeTwo)).get_max(); 
    else
      max = (*Cudd_V(nodeOne)).get_max(); 
    
    Pair newConstant;
    newConstant.set_min(min);
    newConstant.set_max(max);
    temp = Cudd_addConst(gbm,&newConstant);
    Cudd_Ref(temp);
  
    tempOne = *Cudd_V(nodeOne);
    tempTwo = *Cudd_V(nodeTwo);
    

    //Now we get the new theADD 
  
    Approx = Cudd_addApply(gbm,replacePairs,*theADD,temp);
    Cudd_Ref(Approx);
    Cudd_RecursiveDeref(gbm,*theADD);

    *theADD = Approx;
    //    Cudd_Ref(theADD);
    //Cudd_RecursiveDeref(gbm,Approx);

    count = Cudd_DagSize(*theADD); // + Cudd_CountLeaves(theADD);
  
    //Now update the list of constants.
    //We replace the element with the smallest min by the new span constant
    Cudd_RecursiveDeref(gbm,elem1->add);
    elem1->add = temp;

    //Now we need to remove the other element from the list
    removeElement(const_list,elem2);

  }
}
void MDP::compareSubPol(int ovar) {
  int n = orig_vars[ovar].nvals;
  double *mergeStates = new double[n*(n-1)/2];
  findLargestSpan(ovar,mergeStates);
  int i,j,k(0);
  for (i=0; i<n; i++) 
    for (j=i+1; j<n; j++) 
      fprintf(stderr,"%d %d %f\n",i,j,mergeStates[k++]);
  delete [] mergeStates;
    
}
// looks at the OptimalPolicy and computes the 'agreement' between the policies for
// all the pairs of values of original variable ovar
// puts the results in mergeStates, which, if n = orig_vars[ovar].nvals
// should have lenth n*(n-1)/2
void MDP::findLargestSpan(int ovar, double *mergeStates) {
  int i,j, k(0);
  DdNode *pol1, *pol2;
  DdNode *tempres;
  int novvals = orig_vars[ovar].nvals;
  int *listovars = new int[numorigvars-1];
  double spanStates(1);
  for (i=0; i<numorigvars; i++) {
    if (i != ovar) {
      listovars[k++] = i;
      spanStates *= orig_vars[i].nvals;
    }
  }
  k=0;
  for (i=0; i<novvals; i++) { 
    pol1 = restrictVal(gbm,OptimalPolicy,vars,numvars,orig_vars,numorigvars,ovar,i);
    Cudd_Ref(pol1);
    for (j=i+1; j<novvals; j++) {
      pol2 = restrictVal(gbm,OptimalPolicy,vars,numvars,orig_vars,numorigvars,ovar,j);
      Cudd_Ref(pol2);
      tempres = Cudd_addApply(gbm,addBitwiseCompare,pol1,pol2);
      Cudd_Ref(tempres);
      //Cudd_PrintDebug(gbm,tempres,4,100);
      fprintf(stderr,"compared for %d vs %d\n",i,j);
      mergeStates[k] = countSpan(tempres, listovars, numorigvars-1, spanStates);
      k++;
      Cudd_RecursiveDeref(gbm,tempres);
      Cudd_RecursiveDeref(gbm,pol2);
    }
    Cudd_RecursiveDeref(gbm,pol1);    
  }
  delete [] listovars;
}
// compares policies for each of orig variable ovar's branches
// find all pairs of policies that agree
// pick the pair with the smallest value difference and put their indices in mergeStates
int MDP::compareSubPol(int ovar, int *mergeStates) {
  DdNode *pol1, *pol2;
  DdNode *val1, *val2;
  DdNode *tempres, *tempres2;
  int i,j;
  int novvals = orig_vars[ovar].nvals;
  int *compmat = new int[novvals*novvals];
  int nMerge = 0;
  for (i=0; i<novvals*novvals; i++)
    compmat[i] = 0;
  for (i=0; i<novvals; i++) { 
    pol1 = restrictVal(gbm,OptimalPolicy,vars,numvars,orig_vars,numorigvars,ovar,i);
    Cudd_Ref(pol1);
    for (j=i+1; j<novvals; j++) {
      pol2 = restrictVal(gbm,OptimalPolicy,vars,numvars,orig_vars,numorigvars,ovar,j);
      Cudd_Ref(pol2);
      tempres = Cudd_addApply(gbm,addBitwiseCompare,pol1,pol2);
      Cudd_Ref(tempres);
      fprintf(stderr,"compared for %d vs %d\n",i,j);
      Cudd_PrintDebug(gbm,tempres,4,100);
      if (tempres == One) {
	compmat[i*novvals+j] = 1;
	if (i != j) {
	  Cudd_PrintDebug(gbm,pol1,4,100);
	  Cudd_PrintDebug(gbm,pol2,4,100);
	}
      }
      Cudd_RecursiveDeref(gbm,tempres);
      Cudd_RecursiveDeref(gbm,pol2);
    }
    Cudd_RecursiveDeref(gbm,pol1);    
  }
  // now go through compmat and find the pair that agrees and 
  // is closest in value
  double closestVal;
  int firstone = 1;
  for (i=0; i<novvals; i++) 
    for (j=i+1; j<novvals; j++) {
      if (compmat[i*novvals+j]) {
	// restrict value to i, j
	val1 = restrictVal(gbm,OptimalValue,vars,numvars,orig_vars,numorigvars,ovar,i);
	Cudd_Ref(val1);
	val2 = restrictVal(gbm,OptimalValue,vars,numvars,orig_vars,numorigvars,ovar,j);
	Cudd_Ref(val2);
	// subtract restricted
	tempres = Cudd_addApply(gbm,Cudd_addMinus,val1,val2);
	Cudd_Ref(tempres);
	// take absolute value by squaring (no absolute value fuction?)
	tempres2 = Cudd_addApply(gbm,Cudd_addTimes,tempres,tempres);
	Cudd_Ref(tempres2);
	Cudd_RecursiveDeref(gbm,tempres);
	// get maximum leaf value 
	tempres = Cudd_addFindMin(gbm,tempres2);
	Cudd_Ref(tempres);
	Cudd_RecursiveDeref(gbm,tempres2);
	if (firstone || (*Cudd_V(tempres)).get_min() < closestVal) {
	  firstone = 0;
	  mergeStates[0] = i;
	  mergeStates[1] = j;
	  nMerge = 2;
	  closestVal = (*Cudd_V(tempres)).get_min();
	}
      }
    }
  return nMerge;
}

// compares the ADDs corresponding to orig variable ovar's branches in the value function
// returns 1 if any are equal sup.norm to within tol
int MDP::compareValueVals(int ovar, double tol, int **compmat, int *nov)
{
  int rres;
  *nov = orig_vars[ovar].nvals;
  int nnov = *nov;
  *compmat = new int[nnov*nnov];
  if (OptimalValue) 
    rres = compareVarVals(gbm, OptimalValue, vars, numvars, orig_vars, numorigvars, ovar, tol, *compmat) ;
  return rres;
}
int MDP::includedIn(DdManager *dd,int nsuppvars,char **Support){

  int i,j;

  for(i=0;i<numvars;i++){
    for(j=0;j<nsuppvars;j++){
      if(!strcmp(Support[j],vars[i].name))
	break;
      else if(j==nsuppvars-1)
	{
	Cudd_addIthVar(dd,vars[i].add_var->index);
	//fprintf(stderr,"We are adding a var:%d\n",vars[i].add_var->index);
	}
    }
  }
  return 1;
}


int MDP::fillRnumDdNode(DdManager *dd)
{
  int i;
  
  for(i=0;i<numvars;i++){
    //variables[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,2*i+1));
    //prime_variables[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,2*i));
    vars[i].add_var = Cudd_addIthVar(gbm,Cudd_ReadPerm(dd,vars[i].add_var->index)); //Cudd_ReadVars(dd,Cudd_ReadPerm(dd,vars[i].add_var->index));
    prime_vars[i].add_var = Cudd_addIthVar(dd,Cudd_ReadPerm(dd,prime_vars[i].add_var->index)); //Cudd_ReadVars(dd,Cudd_ReadPerm(dd,prime_vars[i].add_var->index));
  }
  return 1;
}

//We copy the read structures (structIn), into the one that we will use throught the rest of the program (structOut)
int MDP::structCopy(rnum *varsIn,rnum *primeVarsIn,onum *origVarsIn)
{

  int i,j,nvals;
  
  for(i=0;i<numvars;i++){
 
    vars[i].name = strdup(varsIn[i].name);
    vars[i].number = varsIn[i].number;
    vars[i].orig_num = varsIn[i].orig_num;
    /* why was this here? - */
    //vars[i].add_var = (DdNode *)malloc(sizeof(DdNode));
    //vars[i].add_var->index = varsIn[i].add_var->index; 
    vars[i].add_var = varsIn[i].add_var;
    Cudd_Ref(vars[i].add_var);

    prime_vars[i].name = strdup(primeVarsIn[i].name);
    prime_vars[i].number = primeVarsIn[i].number;
    prime_vars[i].orig_num = primeVarsIn[i].orig_num;
    //prime_vars[i].add_var = (DdNode *)malloc(sizeof(DdNode));
    //prime_vars[i].add_var->index = primeVarsIn[i].add_var->index;
    prime_vars[i].add_var = primeVarsIn[i].add_var;
    Cudd_Ref(prime_vars[i].add_var);
  }

  for(i=0;i<numorigvars;i++){
    orig_vars[i].name = strdup(origVarsIn[i].name);
    orig_vars[i].nvals = origVarsIn[i].nvals;
    
    nvals = origVarsIn[i].nvals;
    for(j=0;j<nvals;j++){
      orig_vars[i].valname[j] = strdup(origVarsIn[i].valname[j]);
    }
    
    orig_vars[i].nbvars = origVarsIn[i].nbvars;
    orig_vars[i].var1index = origVarsIn[i].var1index;
  }
  
  return 1;
}

//Saves a binary files of vars, prime_vars and orig_vars structures and appends the two Optimal ADDS (value and policy)
//NOTE: the binary version of dddmp_cuddAddArrayStore is not written yet....so we save in text format.
//NEED

int MDP::writeDualOptimal(char *binfilename, DdNode *action, DdNode *value, 
			  rnum *varsT,rnum *prime_varsT, onum *orig_varsT)
{
  FILE *fp;
  int i,j,res;
  DdNode **List;
  char *basepath = "./";
  size_t length;

  
  //now let's try to save them in one file as a vector of ADD
  List  = (DdNode **)malloc(2*(sizeof(DdNode *)));
  
  List[0] = action;
  Cudd_Ref(action);
  List[1] = value;
  Cudd_Ref(value);

  //Names of the ADDs in the file saved on disk
  char **addnames = (char **) malloc(2*(sizeof(char *)));
  for(i=0;i<2; i++)
    addnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
  strcpy(addnames[0],"OptimalPolicy");
  strcpy(addnames[1],"OptimalValue");
  

  int *varid = (int *) malloc(2*numvars*(sizeof(int)));
  for(i=0;i<numvars; i++){
    varid[2*i] = varsT[i].number;
    varid[2*i+1] = prime_varsT[i].number;
  }

  /* prepares binary file for writing */ 
  if ( (fp = fopen(binfilename, "wb")) != NULL ) 
    { 
      
    
      //I put all the structures in the bin file that I am opening
      
      //We write the number of variables in the domain......
      res = fwrite(&numvars, sizeof(int), 1, fp); 
      res += fwrite(&numorigvars,sizeof(int),1,fp);
      
      
      //We write the number of actions in the domain
      res += fwrite(&numactions , sizeof(int), 1, fp); 
      //write the action names
      for(i=0;i<numactions;i++){
      	length = strlen(actionnames[i]) + 1 ;
      	res += fwrite(&length,sizeof(size_t),1,fp);
      	res += fwrite(actionnames[i],sizeof(char),length,fp);
      }
      //I write every rnum_ structure...need to add the DdNode * later
      for(i=0;i<numvars;i++){
      	length = strlen(vars[i].name) + 1 ;
      	res += fwrite(&length, sizeof(size_t), 1, fp); 
      	res += fwrite(varsT[i].name, sizeof(char), length, fp); 
      	res += fwrite(&varsT[i].number, sizeof(int), 1, fp); 
      	res += fwrite(&varsT[i].orig_num, sizeof(int), 1, fp); 
      	//I am printing the index field of the DdNode * in the file now.....
      	res += fwrite(&varsT[i].add_var->index, sizeof(DdHalfWord), 1, fp);
	

      }

      //I write every rnum_ structure...need to add the DdNode * later
      for(i=0;i<numvars;i++){
      	length = strlen(prime_varsT[i].name) + 1;
      	res += fwrite(&length, sizeof(size_t), 1, fp); 
      	res += fwrite(prime_varsT[i].name, sizeof(char), length, fp); 
      	res += fwrite(&prime_varsT[i].number, sizeof(int), 1, fp); 
      	res += fwrite(&prime_varsT[i].orig_num, sizeof(int), 1, fp); 
      	//I am printing the index field of the DdNode * in the file now.....
      	res += fwrite(&prime_varsT[i].add_var->index, sizeof(DdHalfWord), 1, fp);
      }

      //I write every onum_ structure...need to add the DdNode * later
      for(i=0;i<numorigvars;i++){
      	length = strlen(orig_varsT[i].name) + 1;
      	res += fwrite(&length, sizeof(size_t), 1, fp); 
      	res += fwrite(orig_varsT[i].name, sizeof(char), length, fp); 
      	res += fwrite(&orig_varsT[i].nvals, sizeof(int), 1, fp); 
      	for(j=0;j<orig_varsT[i].nvals;j++){
      	  length = strlen(orig_varsT[i].valname[j]) + 1;
      	  res += fwrite(&length, sizeof(size_t), 1, fp); 
      	  res += fwrite(orig_varsT[i].valname[j], sizeof(char), length, fp); 
      	}
      	res += fwrite(&orig_varsT[i].nbvars, sizeof(int), 1, fp); 
      	res += fwrite(&orig_varsT[i].var1index, sizeof(int), 1, fp); 
      }
            

      //We write the ADD in TEXT format
      Dddmp_cuddAddArrayStore(gbm,"DualAdds",2,List,addnames,varnames,varid,DDDMP_MODE_TEXT,DDDMP_VARNAMES,basepath,fp);
   
      fclose(fp); 

    }
  
  
  if(addnames){
    for(i=0;i<2;i++){
      free(addnames[i]);
    }
    free(addnames);
  }
  
  Cudd_RecursiveDeref(gbm,action);
  Cudd_RecursiveDeref(gbm,value);

  free(varid);
  free(List);

  return res;
}


int MDP::readDualOptimal(char *filename) 
{
  DdNode **UPList= (DdNode **)malloc(2*(sizeof(DdNode *)));
  readDualOptimal(gbm,filename,&UPList);
  OptimalPolicy = UPList[0];
  Cudd_Ref(OptimalPolicy);
  OptimalValue = UPList[1];
  Cudd_Ref(OptimalValue);
}
// This function reads from a binary file with the rnum structures and the Dual ADDs
// It needs to do error checking for the case where it compares the variables with what is in the current gbm
int MDP::readDualOptimal(DdManager *dd,char *binfilename,DdNode ***DualOptimal)
{
  FILE *fp;
  int i,j,res,compres,tmp;
  
  int numa,numv,numov; //variables used until bounding the values to numvars, numactions, numorigvars
  char *temp;
  int *varid;
  char **actionREAD;
  char *basepath = "./";
  rnum varsREAD[MAXVARS];
  rnum prime_varsREAD[MAXVARS];
  onum orig_varsREAD[MAXVARS];
  size_t length;
  
  //variables to use in the HeaderLoad thing.....
  int nvars;
  int nsuppvars;
  int nRoots;
  char **orderedVarNames;
  char ** suppVarNames;
  int *varIDs;
  Dddmp_DecompType ddtype;


  char **addnames = (char **) malloc(2*(sizeof(char *)));
  for(i=0;i<2; i++)
    addnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
  
  strcpy(addnames[0],"OptimalPolicy");
  strcpy(addnames[1],"OptimalValue");
 

  /* prepares binary file for reading */ 
  if ( (fp = fopen(binfilename, "rb")) != NULL ) { 
    //read the number of variables in the domain
    res = fread(&numv,sizeof(int),1,fp);
    //read the number of original variables in the domain
    res += fread(&numov,sizeof(int),1,fp);
    //read the number of actions
    res += fread(&numa,sizeof(int),1,fp);
    
    //I am in the case where I don't even need to read the file
    if( ((numv != numvars) && (numvars != 0)) || ((numa != numactions) && (numactions != 0)) || ((numov != numorigvars) && (numorigvars != 0 ))) {
      //      fprintf(stderr,"******ERROR: PROBLEM WITH THE DOMAIN: not the same as the one uploaded already\n");
      return 0;
    }

    //read the actionnames and fill the structure with it
    actionREAD = (char**)malloc(numa *(sizeof(char *)));
    

    for(i=0;i<numa;i++){
      res += fread(&length, sizeof(size_t), 1, fp);
      actionREAD[i] = (char*)malloc(length);
      res += fread(actionREAD[i],sizeof(char),length,fp);
    }
      
    for(i=0;i<numv;i++){
      res += fread(&length, sizeof(size_t), 1, fp);
      varsREAD[i].name = (char*)malloc(length);
      res += fread(varsREAD[i].name, sizeof(char), length, fp); 
      res += fread(&varsREAD[i].number, sizeof(int),1, fp); 
      res += fread(&varsREAD[i].orig_num, sizeof(int), 1, fp); 
      //I am reading the index field of the DdNode * in the file now.....
      varsREAD[i].add_var = (DdNode *)malloc(sizeof(DdNode));
      res += fread(&varsREAD[i].add_var->index, sizeof(DdHalfWord), 1, fp);
    }

    for(i=0;i<numv;i++){
      res += fread(&length, sizeof(size_t), 1, fp); 
      prime_varsREAD[i].name = (char*)malloc(length);
      res += fread(prime_varsREAD[i].name, sizeof(char), length, fp); 
      res += fread(&prime_varsREAD[i].number, sizeof(int),1, fp); 
      res += fread(&prime_varsREAD[i].orig_num, sizeof(int), 1, fp); 
      //I am printing the index field of the DdNode * in the file now.....
      prime_varsREAD[i].add_var = (DdNode *)malloc(sizeof(DdNode));
      res += fread(&prime_varsREAD[i].add_var->index, sizeof(DdHalfWord), 1, fp);
    }
    
    for(i=0;i<numov;i++){
      res += fread(&length, sizeof(size_t), 1, fp); 
      orig_varsREAD[i].name = (char*)malloc(length);
      res += fread(orig_varsREAD[i].name, sizeof(char), length, fp);
      res += fread(&orig_varsREAD[i].nvals, sizeof(int),1, fp); 
      for(j=0;j<orig_varsREAD[i].nvals;j++){
	res += fread(&length, sizeof(size_t), 1, fp); 
	orig_varsREAD[i].valname[j] = (char*)malloc(length);
	res += fread(orig_varsREAD[i].valname[j], sizeof(char), length, fp); 
      } 
      res += fread(&orig_varsREAD[i].nbvars, sizeof(int), 1, fp); 
      res += fread(&orig_varsREAD[i].var1index, sizeof(int), 1, fp); 
    }

    //If the varsCOMP is NULL, then I need to write what I read from the file in it.
    //We do this before reading the ADDs since we might not need to read in the ADDs
        
    //This is the first domain that I enter
    if(numvars == 0){
	
      //      fprintf(stderr,"We are uploading the initial domain\n");
      //We just copy the read structures into the struct that we use
      numvars = numv;
      numorigvars = numov;
      numactions = numa;
      
      //Put the actionnames in the structure
      actionnames = (char**)malloc(numa *(sizeof(char *)));
      for(i=0;i<numa;i++)
	actionnames[i] = strdup(actionREAD[i]);
      
      //copy the structures that I just read into the permanent structures
      structCopy(varsREAD,prime_varsREAD,orig_varsREAD);
      
      //prepare to rewind the file
      //obtains the current value of the file position indicator for FILE *fp
      long filePosition = ftell(fp);
      
      //Here I read the header of the dump file to make sure that the manager had all the info needed.
      Dddmp_cuddHeaderLoad(&ddtype,&nvars,&nsuppvars,&nRoots,&orderedVarNames, &suppVarNames,&varIDs,binfilename, fp);
      
      //I include the variables that are not in the support.
      if(nsuppvars != numvars)
	includedIn(dd,nsuppvars,suppVarNames);
      
      //I go back to where I was before calling cuddHeaderLoad
      fseek(fp,filePosition,0);

    }else{   
      //There already is a domain that was uploaded, now I need to compare this new one
      //I am checking if the structures are the same, otherwise I need to return 0

      //      fprintf(stderr,"We are uploading a second structure\n");
      //compare actionnames
      for(i=0;i<numa;i++){
	if(strcmp(actionREAD[i],actionnames[i])){
	  fprintf(stderr,"THE TWO STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
	}
      }

      //compare varsREAD-vars
      if(!compareRnumStruct(vars,varsREAD,numvars))  {
	//fprintf(stderr,"THE TWO VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
      }

      //compare prime_varsREAD-prime_vars
      if(!compareRnumStruct(prime_vars,prime_varsREAD,numvars)) {
	//fprintf(stderr,"THE TWO PRIME_VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
      }


      //compare orig_varsREAD-orig_vars
      if(!compareOnumStruct(orig_vars,orig_varsREAD, numorigvars)) {
	//fprintf(stderr,"THE TWO ORIG_VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	return 0;
      }
      
    }
    varid = (int *) malloc(2*(numvars)*(sizeof(int)));
    for(i=0;i<numvars; i++){
      varid[2*i] = varsREAD[i].number;
      varid[2*i+1] = prime_varsREAD[i].number;
    }
    //Loading the array of ADDs that form the graph structure to look.
    Dddmp_cuddAddArrayLoad(dd,DDDMP_ROOT_MATCHNAMES,addnames,DDDMP_VAR_MATCHIDS,orderedVarNames,varid,varid,DDDMP_MODE_TEXT,basepath,fp,DualOptimal);
    fclose(fp); 
      
  }
  
  //We fill the Rnum structures with the (DdNode *) that we are looking for....seems to work
  // this doesn't work.
  fillRnumDdNode(dd);

  if(addnames){
    for(i=0;i<2;i++){
      free(addnames[i]);
    }
    free(addnames);
  }
  
  if(actionREAD){
    for(i=0;i<numa;i++)
      free(actionREAD[i]);
    free(actionREAD);
  }
  
  free(varid);
  
  return res;
}
void MDP::shuffleRandom() 
{
  /* generate random order */
  int i;
  double *infor = (double *) malloc(numvars*sizeof(double));
  int *reordList = (int *) malloc(2*numvars*sizeof(int));
  
  for (i=0; i<numvars; i++) {
    infor[i] = ((double) rand())/((double) RAND_MAX);
    //fprintf(stderr,"%f ",infor[i]);
  }
  //fprintf(stderr,"\n");
  
  infoSorter(infor,reordList,numvars);
  (void) Cudd_ShuffleHeap(gbm,reordList);    
  free(infor);
  free(reordList);
}
void MDP::writeOptimalVPDot(char *fileprefix) {
  writeVPDot(OptimalValue, OptimalPolicy, fileprefix);
}
void MDP::writeApproximateVPDot(char *fileprefix) {
  writeVPDot(ApproximateValue, ApproximatePolicy, fileprefix);
}
void MDP::writeDdDot(DdNode *thedd, char *filename) {
  FILE *fileh;
  fileh = fopen(filename,"w");
  if (mvinput) {
    DumpDot_p(gbm, thedd, vars, prime_vars, numvars, orig_vars, numorigvars, NULL, fileh);
  } else {
    Array[0] = thedd;  
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars],NULL,fileh);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(fileh);
}
void MDP::writeMDPDot(char *fileprefix) {
  int i,j;
  FILE *fileh;
  char fname[256];

  // reward function
  sprintf(fname,"%sreward.dot",fileprefix);
  fileh = fopen(fname,"w");
  if (mvinput) {
    DumpDot(gbm, RewardD, vars, numvars, orig_vars, numorigvars, NULL, fileh);
  } else {
    Array[0] = RewardD;  
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars],NULL,fileh);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(fileh);

  //Dual action diagrams
  for(i=0;i<numactions;i++) {
    for (j=0; j<numorigvars; j++) {
      sprintf(fname,"%sCPT%d-%d.dot",fileprefix,i,j);
      fileh = fopen(fname,"w");
      if (mvinput) {
	DumpDot_p(gbm, NewPrime[i][j], vars, prime_vars, numvars, orig_vars, numorigvars, NULL, fileh);
      } else {
	Array[0] = NewPrime[i][j];  
	Cudd_Ref(Array[0]);
	My_DumpDot(gbm,1,Array,varnames,&lnames[numvars],NULL,fileh);  
	Cudd_RecursiveDeref(gbm,Array[0]);
      }
      fclose(fileh);
    }
  }  
}
void MDP::writeVPDot(DdNode *valADD, DdNode *polADD,char *fileprefix) {
  FILE *value, *action;
  char fname[256];
  // write out data files
  sprintf(fname,"%svalue.dot",fileprefix);
  value = fopen(fname,"w");
  if (mvinput) {
    DumpDot(gbm, valADD, vars, numvars, orig_vars, numorigvars, NULL, value);
  } else {
    Array[0] = valADD;  
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars],NULL,value);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(value);
  
  sprintf(fname,"%spolicy.dot",fileprefix);
  
  action = fopen(fname,"w");
  
  if (mvinput) {
    DumpDot(gbm, polADD, vars, numvars, orig_vars, numorigvars, actionnames,action);
  } else {
    Array[0] = polADD;
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars+1],actionnames,action);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(action);
}


/* 
   Convergence test to be used in conjuncture with Cudd_addApply
   Determines if the pair ADD has converged.  
   We are using a very conservative criterion: if all intervals overlap
   or are within tolerance eps of each other, we have converged.
*/

DdNode * Convergence_Test(DdManager * gbm, DdNode ** f, DdNode ** g)
{
 
  DdNode *F, *G;
  DdNode *oneD,*zeroD;
  F = *f; G = *g;
    
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double min1 = (*Cudd_V(F)).get_min();
    double min2 = (*Cudd_V(G)).get_min();
    double max1 = (*Cudd_V(F)).get_max();
    double max2 = (*Cudd_V(G)).get_max();
    
    /*for the case where there is no approximation*/
    /*when there is approximation being done*/
    if( ((min1 <= min2) && (min2 <= max1)) || ((min2 <= min1) && (min1 <= max2)) || ((min2 >= max1) && ((min2 - max1) <= tolerance)) || ((max2 <= min1) && ((min1 - max2) <= tolerance))){
      oneD = Cudd_addConst(gbm,pOne);
      return (oneD);
    }else{
      zeroD = Cudd_addConst(gbm,pZero);
      return(zeroD);
    } 
  } 
  return(NULL);
} 
DdNode *addBitwiseCompare(DdManager *dd , DdNode **f, DdNode **g)
{
    DdNode *F, *G;

    F = *f; G = *g;

    if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
      double fmn = (*Cudd_V(F)).get_min();
      double gmn = (*Cudd_V(G)).get_min();
      double result(0);
      int fival = ((int) fmn);
      int gival = ((int) gmn);
      result = fival & gival;
      /*
	// have to do this with doubles for lots of actions
      while (result < 1 && fival >= 1 && gival >= 1) { 
	if (fival%2 == gival%2)
	  result = 1;
	fival = fival >> 1;
	gival = gival >> 1;
      }
      */
      if (result < 1) 
	return (Zero);
      else
	return (One);
      
    } 
    return (NULL);

}
DdNode * reducePrecision(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;
  
  if (F == DD_MINUS_INFINITY(dd)) return(G);
  if (G == DD_MINUS_INFINITY(dd)) return(F);
  
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double fmx = (*Cudd_V(F)).get_min();
    double gmx = (*Cudd_V(G)).get_min();
    double nv = (double) ((floor(fmx*gmx))/gmx);
    Pair newPair;
    newPair.set(nv);
    DdNode *res = Cudd_addConst(dd,&newPair);
    return (res);
  }
  return (NULL);
}

DdNode *My_addMaximum(DdManager * dd, DdNode ** f, DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;

    if (F == DD_MINUS_INFINITY(dd)) return(G);
    if (G == DD_MINUS_INFINITY(dd)) return(F);

    if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
        if ((*Cudd_V(F)) >= (*Cudd_V(G))) {
          return(F);
        } else if ((*Cudd_V(G)) > (*Cudd_V(F))) {
          return(G);
        } else {
          Pair newPair;
          double mn,mx;
          double fmx = (*Cudd_V(F)).get_max();
          double fmn = (*Cudd_V(F)).get_min();
          double gmx = (*Cudd_V(G)).get_max();
          double gmn = (*Cudd_V(G)).get_min();
          if (fmx >= gmx)
            mx = fmx;
          else 
            mx = gmx;
          if (fmn >= gmn)
            mn = fmn;
          else
            mn = gmn;

          newPair.set_both(mn,mx);
          DdNode *res = Cudd_addConst(dd,&newPair);
          return (res);
          
        }
    }
    return(NULL);

} /* end of My_addMaximum */


// other, non-MDP functions
//removes an element from the constant list
void  removeElement(list_cst *list,list_cst *element){
  list_cst *curr,*suivant;
  
  curr = list;
  suivant = list->next;

  while(suivant != NULL){

    if(suivant == element)
      {
	curr->next = suivant->next;
	Cudd_RecursiveDeref(gbm,suivant->add);
	free(suivant);
	break;
      }else{
	curr = suivant;
	suivant = suivant->next;
      }
  }
}

//takes in a ADD and a constant ADD node and replaces the 2 nodes that I want by the ADD constant
DdNode *replacePairs(DdManager *dd, DdNode **f, DdNode **newC)
{
  DdNode *F, *ROE;
  
  F = *f; ROE = *newC;
  //Look for a constant that equals 
  //PairOne and PairTwo are the 2 constants that need to be replaced by the newC constant node
  if (Cudd_IsConstant(F)){
    if( ((*Cudd_V(F)) == tempOne) || ((*Cudd_V(F)) == tempTwo) ) {
      return ROE;
    }else{
      return F;
    } 
  }else{
    return(NULL);
  }
}


/*expects two constants f,g, rounds off f to within trunc and returns the result */ 

DdNode * myRoundOff(DdManager *dd, DdNode **f, DdNode **roe)
{
  DdNode *tmp, *res, *F, *ROE;
  Pair value;

  F = *f; ROE = *roe;
  //divides so that it is between integers
  if (Cudd_IsConstant(F) && Cudd_IsConstant(ROE)) {
    tmp = Cudd_addDivide(dd, f, roe);
    Cudd_Ref(tmp);
    
    //since we are using integer round-off
    if ((floor((*Cudd_V(tmp)).get_min())) == (floor((*Cudd_V(tmp)).get_max()))) {
      avgError.span(*Cudd_V(F));
      cntError++;
      // why was this LARGEFLAG thing here????
      //value.set_min(2*LARGEFLAG + floor((*Cudd_V(tmp)).get_min()));
      //value.set_max(2*LARGEFLAG + ceil((*Cudd_V(tmp)).get_max()));
      value.set_min(floor((*Cudd_V(tmp)).get_min()));
      value.set_max(ceil((*Cudd_V(tmp)).get_max()));
      /* deal with the fence-sitters*/
      if(value.get_min() == value.get_max())
	value.set_max(value.get_max()+1);
      Cudd_RecursiveDeref(dd, tmp);
      res = Cudd_addConst(dd,&value);
      return (res);
    } else {
      Cudd_Deref(tmp);
      return (tmp);                // This added late thursday night in a rush - right?
    }
  }
  return (NULL);
}
DdNode * getTotalSpan(DdManager *dd, DdNode **f, DdNode **x)
{
  Pair pRes;
  DdNode *dMax, *dMin;
  
  
  if (Cudd_IsConstant(*f)){
    if ((*Cudd_V(*f)) == (*pOne)) {
      dMax = Cudd_addFindMax(gbm, *x);
      Cudd_Ref(dMax);
      dMin = Cudd_addFindMin(gbm, *x);
      Cudd_Ref(dMin);
      pRes.set_max((*(Cudd_V(dMax))).get_max());
      pRes.set_min((*(Cudd_V(dMin))).get_min());

      Cudd_RecursiveDeref(gbm,dMax);
      Cudd_RecursiveDeref(gbm,dMin);
      avgError.span(pRes);
      return(Cudd_addConst(dd,&pRes));
    }
    return(Zero);
  } 
  return (NULL);     
} 
DdNode *  My_cuddAddApplyRecur(DdNode * (*op)(DdManager *, DdNode **, DdNode **), 
			       DdNode * f, DdNode * g)
{
  DdManager *dd = gbm;
  DdNode *res,*fv, *fvn, *gv, *gvn, *T, *E;
  unsigned int ford, gord;
  unsigned int index;
  
  /* Check terminal cases. Op may swap f and g to increase the
   * cache hit ratio.
   */
  res = (*op)(dd,&f,&g);
  if (res != NULL) return(res);
  
  /* Recursive Step */
  ford = Cudd_ReadPerm(dd,f->index);
  gord = Cudd_ReadPerm(dd,g->index);
  if (ford <= gord) {
    index = f->index;
    fv = Cudd_T(f);
	fvn = Cudd_E(f);
    } else {
	index = g->index;
	fv = fvn = f;
    }
    if (gord <= ford) {
	gv = Cudd_T(g);
	gvn = Cudd_E(g);
    } else {
	gv = gvn = g;
    }

    T = My_cuddAddApplyRecur(op,fv,gv);
    if (T == NULL) 
      return(NULL);
    Cudd_Ref(T);

    E = My_cuddAddApplyRecur(op,fvn,gvn);
    if (E == NULL) {
      Cudd_RecursiveDeref(dd,T);
      return(NULL);
    }
    Cudd_Ref(E);

    //res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    // NEEDS TO BE FIXED - DON'T HAVE ACCESS TO VARS HERE ******************************
    //    res = Cudd_addIte(gbm,vars[(int) (index/2)].add_var,T,E);

    if (res == NULL) {
      Cudd_RecursiveDeref(dd, T);
      Cudd_RecursiveDeref(dd, E);
      return(NULL);
    }
    Cudd_Deref(T);
    Cudd_Deref(E);

    return(res);

} 

// My_addApply
DdNode * My_addApply( DdNode * (*op)(DdManager *, DdNode **, DdNode **),DdNode * f,  DdNode * g)
{
    DdNode *res;
    res = My_cuddAddApplyRecur(op,f,g);
    return(res);

}
 
//converts the ADD f to a 0-1 ADD by setting every leaf in f 
//which is equal to the constant c to 1, and everything else to 0 
//For use with addApply

DdNode * addAddExact(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *G, *res;
  F = *f;
  G = *g;
  double fval, gval;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    fval = (*Cudd_V(F)).get_min();
    gval = (*Cudd_V(G)).get_min();
    if (fval == gval) 
      res = Cudd_addConst(dd,pOne);
    else
      res = Cudd_addConst(dd,pZero);
    return (res);
  } 
  return (NULL);
}
Pair  get_span(DdNode *x) {
  Pair *res;
  res = (Pair *)malloc(sizeof(Pair));
  (*res).set_max((*Cudd_V(Cudd_addFindMax(gbm,x))).get_max());
  (*res).set_min((*Cudd_V(Cudd_addFindMin(gbm,x))).get_min());
  return *res;
}
double get_extent(DdNode *x) 
{
  double error;
  DdNode *mx,*mn;
  mx = Cudd_addFindMax(gbm,x);
  Cudd_Ref(mx);
  mn = Cudd_addFindMin(gbm,x);
  Cudd_Ref(mn);
  
  error = ((*(Cudd_V(mx))).get_max() - (*(Cudd_V(mn))).get_min());
  Cudd_RecursiveDeref(gbm,mx);
  Cudd_RecursiveDeref(gbm,mn);
  return error;
}

double get_error(DdNode *x){

  DdNode *errorAdd,*maxNode;
  double errorD;

  errorAdd = Cudd_addApply(gbm,getErrorAdd,x,One);
  Cudd_Ref(errorAdd);
  maxNode = Cudd_addFindMax(gbm,errorAdd);
  Cudd_Ref(maxNode);
  Cudd_RecursiveDeref(gbm,errorAdd);
  errorD = (*Cudd_V(maxNode)).get_max();
  //printf("Error:%f\n",errorD);
  Cudd_RecursiveDeref(gbm,maxNode);
  return errorD;
}

DdNode * getErrorAdd(DdManager *dd, DdNode **f, DdNode **x)
{
  Pair error;
  double errorD;
  DdNode *res;
  
  if (Cudd_IsConstant(*f)){
    errorD = (*Cudd_V(*f)).get_max() -(*Cudd_V(*f)).get_min();
    //printf("errorD:%f\n",errorD);
    error.set(errorD);
    res = Cudd_addConst(gbm,&error);
    return(res);
  } 
  return (NULL);     
} 


double  length_span(DdNode *x,DdNode *y){
  double length;
  double min,max;

  if((*Cudd_V(x)).get_min() <= (*Cudd_V(y)).get_min())
    min = (*Cudd_V(x)).get_min(); 
  else
    min = (*Cudd_V(y)).get_min(); 

  if((*Cudd_V(x)).get_max() <= (*Cudd_V(y)).get_max())
    max = (*Cudd_V(y)).get_max(); 
  else
    max = (*Cudd_V(x)).get_max(); 
      
  length = max-min;
  return length;
}




//returns 0 if there is a discrepancy between the structures..... 1 if everything is fine.
int compareRnumStruct(rnum varsP[],rnum varsT[], int nv)
{
  int i;

  for(i=0;i< nv;i++){

    if( strncmp(varsP[i].name,varsT[i].name,128) || (varsP[i].number != varsT[i].number) || (varsP[i].orig_num != varsT[i].orig_num) || (varsP[i].add_var->index != varsT[i].add_var->index)) 
      return 0;
  }
  return 1;
}

//Compare two structures of type onum....returns 1 if is it the same, 0 otherwise
//NOTE it assumes numvars is a global variable.....we might have to add it as a parameter if it is not the case
int compareOnumStruct(onum orig_varsP[],onum orig_varsT[], int nov)
{
  int i,j;

  for(i=0;i<nov;i++){
    if( strncmp(orig_varsP[i].name,orig_varsT[i].name,128) || (orig_varsP[i].nvals != orig_varsT[i].nvals) || (orig_varsP[i].nbvars != orig_varsT[i].nbvars) || (orig_varsP[i].var1index != orig_varsT[i].var1index)){ 
      return 0;
    }
    for(j=0;j<(orig_varsP[i].nvals);j++)
      if(strcmp(orig_varsP[i].valname[j],orig_varsT[i].valname[j])){
	return 0;
      }
  }
  return 1;
}


//return the maximum of two integers
int max(int i, int j)
{
  if( i< j) 
    return j;
  else
    return i;
}



// we want this function to mutliply 
DdNode * newPrimeTimes(DdManager *dd, DdNode **f, DdNode **g) {
  DdNode *F, *G;
  Pair pVal;
  
  
  F = *f; G = *g;
  return F;
}
DdNode * getAction(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *G, *F;
  double gval, fval;
  F = *f;
  G = *g;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    fval = (*Cudd_V(F)).get_min();
    gval = (*Cudd_V(G)).get_min();
    if (fval == gval)
      return(One);
    else
      return(Zero);
  }
  return(NULL);
}
DdNode * PickAction(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *res;
  double fval;
  Pair pVal;
  F = *f;
  if (Cudd_IsConstant(F)) {
    fval = (*Cudd_V(F)).get_min();
    pVal.set((double) pickAction(fval));
    res = Cudd_addConst(gbm,&pVal);
    return(res);
  }
  return(NULL);    
}
// converts a policy where each leaf is a numbered action 0....numactions
// back to the binary format where leaf in binary has a 1 at position i if action i is optimal
DdNode * BinaryAction(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *res;
  double fval;
  int ifval;
  Pair pVal;
  F = *f;
  if (Cudd_IsConstant(F)) {
    fval = (*Cudd_V(F)).get_min();
    ifval = int(fval);
    pVal.set((double) binaryAction(ifval));
    res = Cudd_addConst(gbm,&pVal);
    return(res);
  }
  return(NULL);    
}
DdNode * MyPlus(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  double gval,fval;
  Pair pVal;
  
  
  F = *f; G = *g;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    gval = (*Cudd_V(G)).get_min();
    fval = (*Cudd_V(F)).get_min();
    
    if (fval == 0.0) return(G);
    
    return(F);
  }
  return(NULL);

}

DdNode * addMean(DdManager * dd, DdNode ** f, DdNode ** g)
{
    DdNode *F, *res;
    double fmn,fmx;
    Pair pVal;


    F = *f; 
    if (Cudd_IsConstant(F)) {
      fmn = (*Cudd_V(F)).get_min();
      fmx = (*Cudd_V(F)).get_max();
      
      pVal.set((fmn+fmx)/2.0);
      res = Cudd_addConst(gbm,&pVal);
      return(res);
    }
    return(NULL);

}

// take all vcurrent and action diag and merges them by picking the right action
DdNode *actionMerge(DdNode **valD, DdNode *actD, int numact)
{
  int i;
  DdNode *pickAct, *temp, *temp1, *sum;
  Pair iPair;

  pickAct = Cudd_addConst(gbm,pOne);
  Cudd_Ref(pickAct);
  
  sum = Cudd_addConst(gbm,pZero);
  Cudd_Ref(sum);
  
  temp1 = Cudd_addConst(gbm,pOne);
  Cudd_Ref(temp1);

  for (i=0; i<numact; i++) {
    
    //iPair.set((double) i+1);
    iPair.set((double) i);
    Cudd_RecursiveDeref(gbm,pickAct);
    pickAct = Cudd_addConst(gbm,&iPair);
    Cudd_Ref(pickAct);

    temp = Cudd_addApply(gbm, addAddExact, actD, pickAct);
    Cudd_Ref(temp);
    
    Cudd_RecursiveDeref(gbm, temp1);
    temp1 = Cudd_addApply(gbm, Cudd_addTimes, temp,valD[i]);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm, temp);

    /** JH FEB 4 **/
    
    temp = Cudd_addApply(gbm, Cudd_addPlus, temp1, sum);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,sum);
    sum = temp;
  }

  Cudd_RecursiveDeref(gbm,temp1);
  Cudd_RecursiveDeref(gbm,pickAct);
    
  return sum;
}

// Wrapper to get around some problems in Cudd
DdNode *Cudd_Else(DdNode *foo)
{
        if (Cudd_IsComplement(foo)) {
                return Cudd_Not(Cudd_E(foo));
        } else {
                return Cudd_E(foo);
        }
}

// Wrapper to get around some problems in Cudd
DdNode *Cudd_Then(DdNode *foo)
{
        if (Cudd_IsComplement(foo)) {
                return Cudd_Not(Cudd_T(foo));
        } else {
                return Cudd_T(foo);
        }
}
/*
Sums over the primed variables in dd from index first to index last
*/
DdNode *sumSubtrees(DdNode *dd, int first, int last,  rnum *prime_vars, int numvars)
{
  DdNode *temp, *cube;
  int i,j;
  DdNode **ArraySum;
  ArraySum = (DdNode **)malloc(numvars*(sizeof(DdNode *)));
  j=0;
  //JH09/05/00  for (i=first-1; i>=last; i--) {
  for (i=first; i<last; i+=2) {
    ArraySum[j] = prime_vars[i/2].add_var;
    Cudd_Ref(ArraySum[j]);
    j++;
  }
  //JH09/05/00 cube = Cudd_addComputeCube(gbm,ArraySum,NULL,first-last);
  cube = Cudd_addComputeCube(gbm,ArraySum,NULL,(last-first)/2);
  Cudd_Ref(cube);

  j=0;
  for (i=first; i<last; i+=2) {
    Cudd_RecursiveDeref(gbm,ArraySum[j]);
    j++;
  }

  temp = Cudd_addExistAbstract(gbm,dd,cube);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,cube);
  free(ArraySum);
  return temp;
}

DdNode *MDP::sumOutAllOrigVars(DdNode *dd)
{
  int j;
  DdNode *temp, *res;
  res = dd;
  Cudd_Ref(res);
  for (j=0; j<numorigvars; j++) {
    temp = sumOutPrime(res,orig_vars+j,vars);
    Cudd_RecursiveDeref(gbm,res);
    res = temp;
  }
  return res;
}

// sums out primed variable ovar
DdNode *sumOutPrime(DdNode *dd, onum *ovar, rnum *prime_vars) 
{
  int first, last;
  first = ovar->var1index;
  last = ovar->var1index+ovar->nbvars;      
  return newSumSubtrees(dd,first,last,prime_vars);
}

/*
Sums over the primed variables in dd from index first to index last
for multi-valued diagrams
*/
DdNode *newSumSubtrees(DdNode *dd, int first, int last, rnum *prime_vars)
{

  DdNode *temp, *cube;
  int i,j;
  DdNode ** ArraySum = (DdNode **)malloc((last-first+1)*(sizeof(DdNode *)));

  j=0;
  for (i=first; i<last; i++) {
    ArraySum[j] = prime_vars[i].add_var;
    Cudd_Ref(ArraySum[j]);
    j++;
  }
  cube = Cudd_addComputeCube(gbm,ArraySum,NULL,(last-first));
  Cudd_Ref(cube);

  j=0;
  for (i=first; i<last; i++) {      // JH-26/05/04 i+=2) {
    Cudd_RecursiveDeref(gbm,ArraySum[j]);
    j++;
  }

  temp = Cudd_addExistAbstract(gbm,dd,cube);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,cube);
  free(ArraySum);
  return temp;
}

/*
Sums over the primed variables in dd with indices 
so its the sum_{i=1...numindices}(prime_vars[indices[i]])
for multi-valued diagrams
*/
DdNode *sumOutPrime(DdNode *dd, int *indices, int numindices, rnum *prime_vars, int numvars)
{
  DdNode *temp, *cube;
  int i;
  DdNode ** ArraySum = (DdNode **)malloc(numindices*(sizeof(DdNode *)));

  for (i=0; i<numindices; i++) {
    ArraySum[i] = prime_vars[indices[i]].add_var;
    Cudd_Ref(ArraySum[i]);
  }
  cube = Cudd_addComputeCube(gbm,ArraySum,NULL,numindices);
  Cudd_Ref(cube);

  for (i=0; i<numindices; i++) 
    Cudd_RecursiveDeref(gbm,ArraySum[i]);


  temp = Cudd_addExistAbstract(gbm,dd,cube);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,cube);
  free(ArraySum);
  return temp;
}
// sums out all the primed variables in toSum
DdNode * sumPrimes(DdNode *toSum, rnum *prime_vars, int numvars) {
  DdNode *temp;
  int i,j;
  int *indices, numindices, *allindices;
  // find indices of primed variables in support  of toSum
  allindices = Cudd_SupportIndex(gbm,toSum);
  
  indices = new int[2*numvars];
  // loop only over prime variables
  j=0;
  numindices = 0;
  for (i=0; i<2*numvars; i+=2) {
    if (allindices[i]) {
      indices[numindices] = j;
      numindices++;
    }
    j++;
  }
  if (numindices > 0) 
    temp = sumOutPrime(toSum,indices,numindices,prime_vars,numvars);
  else 
    temp = NULL;
  delete [] indices;
  /*old 
  Cudd_Ref(temp);
  return temp;
  sum = toSum;
  Cudd_Ref(sum);
  for (i=0; i<numvars; i++) {
    if (prime_var[i] is in the support of sum) {
      temp = newSumSubtrees(sum,i,i+1,prime_vars);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,sum);
      sum = temp;
    }
  }
  */
  return temp;
}

//Counts the number of internal nodes for the tree equivalent to the Ddnode *x
int count_internal_nodes_tree(DdNode *x)
{
  DdNode *tempE, *tempT;

  int counter=0,leftcounter=0,rightcounter=0;
  
  if (!Cudd_IsConstant(x)) {
    tempE = Cudd_Else(x);
    tempT = Cudd_Then(x);
    counter = 1;
    leftcounter = count_internal_nodes_tree(tempT);
    rightcounter = count_internal_nodes_tree(tempE);
  }

  counter += leftcounter+rightcounter;
  return counter;
}

//Counts the number of leaves for a Tree similar to the Ddnode *x
int count_leaves_tree(DdNode *x)
{
  DdNode *tempE, *tempT;

  int counter=0,leftcounter=0,rightcounter=0;
  
  if (Cudd_IsConstant(x)) {
    counter = 1;
  } else {
    tempE = Cudd_Else(x);
    tempT = Cudd_Then(x);

    leftcounter = count_leaves_tree(tempT);
    rightcounter = count_leaves_tree(tempE);
    counter = leftcounter+rightcounter;
  }
  return counter;
}
/* To dump the output as an HTML format*/
void dumpHtmlRew(FILE *log,char *outpath)
{
  char temp[2*MAXLEN];
 
  fprintf(log,"<html>\n\t<head>\n");
  fprintf(log,"<title>Reward</title>");
  fprintf(log,"</head>\n\n\t<body>\n\t<h1>Reward</h1>\n");
  strcpy(temp,"<img src=http://www.research.att.com/~north/cgi-bin/webdot.cgi/http://www.cs.ubc.ca/spider/jhoey/spudd/");
  strcat(temp,outpath);
  strcat(temp,"reward.dot.gif>\n");
  fprintf(log,temp);
  fprintf(log,"\n\n\t<hr>\n</body>\n</html>\n");
} 
void  infoSorter(double *infoList, int *reorderList, int numvars) 
{
  int i,j,minj;
  double maxI=-10000.0, minI;
  /* find maximum value for bookeeping*/
  for (i=0; i<numvars; i++) 
    if (infoList[i] > maxI)
      maxI = infoList[i];
 

  for (i=0; i<numvars; i++) {

    /* find minimum value */
    minI = maxI + 1.0;
    for (j=0; j<numvars; j++) 
      if (infoList[j] < minI) {
	minj = j;
	minI = infoList[j];
      }

    /* get this value out of the way*/
    infoList[minj] = maxI+1.0;

    /*update the unprimed variable order */

    reorderList[2*i+1] = 2*minj+1;
    reorderList[2*i] = 2*minj;
  }
}

// for connecting to supervisor
int connectToSupervisor     (Client **spudd_c, Server **spudd_s) {
  //start a server and wait for a connection
  *spudd_s = new Server (SPUDD_S,SPUDD_S+PORT_RANGE);
  if ( *spudd_s ){
    (*spudd_s)->startServer ( );
    if ( (*spudd_s)->waitForClient ( ) ) {
      //Connect to the server
      sleep(1);
      fprintf(stderr,"connecting to server - \n");
      *spudd_c = new Client(SPUDD_C,SPUDD_C+PORT_RANGE,SUPERVISOR_HOST_TMP );
      if ( *spudd_c && (*spudd_c)->startClient ( ) ) 
	return true;
    }
  }
  return false;
}
// compares ADDs corresponding to all values of orig. variable ovar
// with sup. norm tolerance tol
// fills the novar*novar input compmat with 1s for each values i,j that match
// returns 1 if any two match (if there is any off-diagonal 1s in compmat
int compareVarVals(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		   int ovar, double tol, int *compmat) 
{
  int i,j;
  int twoeq = 0;
  int novvals = ov[ovar].nvals;
  for (i=0; i<novvals*novvals; i++)
    compmat[i] = 0;
  for (i=0; i<novvals; i++) 
    for (j=i; j<novvals; j++) {
      compmat[i*novvals+j] = compareVarVals(dd, val, v, nvars, ov, novars, ovar, i, j, tol);
      if (j>i) 
	twoeq = twoeq + compmat[i*novvals+j];
    }
  return twoeq;
}

// compares the values of an ADD for two values of a variable 
// using supremum norm - returns result
int compareVarVals(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		  int ovar, int ovarval1, int ovarval2, double tol) {
  DdNode *val1, *val2;
  val1 = restrictVal(dd,val,v,nvars,ov,novars,ovar,ovarval1);
  Cudd_Ref(val1);
  val2 = restrictVal(dd,val,v,nvars,ov,novars,ovar,ovarval2);
  Cudd_Ref(val2);
  
  Pair * pTol = new Pair(tol);

  int res = Cudd_EqualSupNorm(dd,val1,val2,pTol,0);
  Cudd_RecursiveDeref(dd,val1);
  Cudd_RecursiveDeref(dd,val2);
  return res;
}
// returns the restricted ADD val to the value of the original variable ovar = ovarval
// doesn't change val
// doesn't ref the return result
DdNode *restrictVal(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, 
		    int novars,
		    int ovar, int ovarval) {

  int j,nbv,nbvals,tmp;
  // figure out the variable values for ovar=ovarval
  nbv = ov[ovar].nbvars;
  nbvals = int(pow(2.0,nbv));
  tmp = nbvals-ovarval-1;
  int *varass = new int[nbv];
  DdNode **varss = new DdNode*[nbv];
  for (j=nbv-1; j>=0; j--) {
    varss[j] = v[ov[ovar].var1index+j].add_var;
    Cudd_Ref(varss[j]);
    //i think this could be done faster
    //with bit moves but i'm trying not 
    //to break too much --jlb 7/6/10
    varass[j] = tmp%2;
    tmp = tmp/2;
  }
  // get cube for that ovar ovarval
  DdNode *oCube = Cudd_addComputeCube(dd,varss,varass,nbv);
  Cudd_Ref(oCube);

  // restrict the add to that value
  DdNode *restrictedVal = Cudd_addRestrict(dd,val,oCube);
  Cudd_RecursiveDeref(dd,oCube);

  for (j=0; j<nbv; j++) 
    Cudd_RecursiveDeref(dd,varss[j]);

  delete [] varass;
  delete [] varss;

  return restrictedVal;
}
void MDP::printAction(int act, FILE *fp)
{
  fprintf(fp,"%s",actionlist[act].name);
}
void MDP::printVarName(int v, FILE *fp)
{
  fprintf(fp,"%s",orig_vars[v].name);
}
void MDP::printVars(FILE *fp)
{
  int i,j;
  fprintf(fp,"(variables ");
  for (i=0; i<numorigvars; i++) {
    fprintf(fp,"(%s ",orig_vars[i].name);
    for (j=0; j<orig_vars[i].nvals; j++)
      fprintf(fp,"%s ",orig_vars[i].valname[j]);
    fprintf(fp,") ");
  }
  fprintf(fp,")");
}
void MDP::printMDPSpudd(FILE *fp) 
{
  int i,j;
  // testing only
  /*
  int *varvals = new int[numorigvars];
  for (i=0; i<numorigvars; i ++)
    varvals[i] = 0;
  tmp = buildCubeOrig(varvals, false);
  Cudd_Ref(tmp);
  for (i=0; i<23; i++) {
    for (j=0; j<numorigvars; j++) {
      tmp = Cudd_addApply(gbm,Cudd_addMinus,NewPrime[i][j],NewPrime[i+23][j]);
      Cudd_Ref(tmp);
      fprintf(stderr,"difference %d %d\n",i,j);
      Cudd_PrintDebug(gbm,tmp,4,100);
      Cudd_RecursiveDeref(gbm,tmp);
    }
  }
  */
  printVars(fp);
  fprintf(fp,"\nunnormalized\n");
  int **varddnames = new int*[numactions];
  int *costnames = new int[numactions];
  int dagCount = 100000;
  int numexpadds, numadds;
  numexpadds = dagCount;
  numadds = 0;
  rnum **addlist = new rnum*[numexpadds];
  
  for (i=0; i<numactions; i++) {
    varddnames[i] = new int[numorigvars];
    // first print all the found sub-adds
    // have to actually add a "j" tag on each one so they don't
    // get confused
    for (j=0; j<numorigvars; j++) {
      fprintf(stderr,"getting dds for action %d/%d var %d/%d --- current numadds %d\n",i,numactions,j,numorigvars,numadds);
      varddnames[i][j] = printDdNode(gbm,&addlist,numadds,numexpadds,NewPrime[i][j],vars,prime_vars,numvars,orig_vars,numorigvars,fp);
    }
    fprintf(stderr,"getting cost function %d\n",i);
    costnames[i] = printDdNode(gbm,&addlist,numadds,numexpadds,actionCostNoDummy[i],vars,prime_vars,numvars,orig_vars,numorigvars,fp);
  }
  fprintf(stderr,"getting reward function\n");
  int rewddname = printDdNode(gbm,&addlist,numadds,numexpadds,RewardDNoDummy,vars,prime_vars,numvars,orig_vars,numorigvars,fp);
  fprintf(stderr,"numadds found %d\n",numadds);
  for (int i=0; i<numadds; i++) {
    fprintf(fp,"dd dd%d\n",i);
    fprintf(fp,"%s",addlist[i]->name);
    fprintf(fp,"\nenddd\n");
  }
      
  
  for (i=0; i<numactions; i++) {
    fprintf(fp,"action ");
    printAction(i,fp);
    fprintf(fp,"\n");
    //printAdd(NewPrime[i][j],fp);
    for (j=0; j<numorigvars; j++) {
      printVarName(j,fp);
      fprintf(fp," (dd%d)\n",varddnames[i][j]-1);
    }
    // since its already negative, and the parser expects +ve costs
    // (which it negates) - we have to multiply by -1 here!
    fprintf(fp,"cost (dd%d)\n",costnames[i]-1);
    fprintf(fp,"endaction\n");
  }
  for (i=0; i<numadds; i++) 
    delete addlist[i];
  delete [] addlist;

  //  printReward(fp);
  fprintf(fp,"reward ");
  fprintf(fp,"(dd%d)\n",rewddname-1);
  printDiscount(fp);
  printTolHor(fp);
}

void MDP::printReward(FILE *fp)
{
  fprintf(fp,"reward ");
  printAdd(RewardD,fp);
}
void MDP::printNewPrime(FILE *fp)
{
  int i,j;
  for (i=0; i<numactions; i++) 
    for (j=0; j<numorigvars; j++) 
      printAdd(NewPrime[i][j],fp);
}
void MDP::printDiscount(FILE *fp)
{
  fprintf(fp,"discount %f\n",discount_factor);
}
void MDP::printTolHor(FILE *fp)
{
  if (horizon < 0) 
    fprintf(fp,"tolerance %f\n",tolerance);
  else
    fprintf(fp,"horizon %lf\n",horizon);
} 
// prints the names and values of the original vars in state
void MDP::printState(int *state, FILE *fp)
{
  int i,j;
  int maxstrlen = -1, slen;
  for (i=0; i<numorigvars; i++) {
    if ((slen = strlen(orig_vars[i].name)) > maxstrlen)
      maxstrlen = slen;
  }
  maxstrlen++;
  for (i=0; i<numorigvars; i++) {
    fprintf(fp,"%s",orig_vars[i].name);
    slen = strlen(orig_vars[i].name);
    for (j=0; j<maxstrlen-slen; j++)
      fprintf(fp," ");
    fprintf(fp,"%s",orig_vars[i].valname[state[i]]);
    fprintf(fp,"\n");
  }
  
  //  for (i=0; i<numorigvars; i++)
  //fprintf(fp,"%s\t",orig_vars[i].valname[state[i]]);
  fprintf(fp,"\n");
}
// updateQ updates the Q functions based on current transition function
// and reward function e.g. Dyna
void MDP::updateQ(DdNode **Q, int acttaken, int *varvals, int k)
{
  int j;
  DdNode *temp1, *temp2, *mQ;
  // recall old Q function
  Cudd_RecursiveDeref(gbm,VcurrentMin[acttaken]);
  VcurrentMin[acttaken] = Q[acttaken];
  Cudd_Ref(VcurrentMin[acttaken]);

  // loop through actions ap, constructing mQ = max_{ap} VcurrentMin(s',ap)
  temp1 = VerySmall;
  Cudd_Ref(temp1);

  for (j=0; j<numactions; j++) {
    temp2=Cudd_addApply(gbm,Cudd_addMaximum,temp1,VcurrentMin[j]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp2;
  }
  // change s -> s' in mQ 
  mQ = Cudd_addSwapVariables(gbm,temp1,Array1,Array2,numvars); 
  Cudd_Ref(mQ); 
  Cudd_RecursiveDeref(gbm,temp1);

  //   compute NewPrimed[actions] *mQ
  //   and sum out primed variables
  temp1 = multiplySumSet(mQ,NewPrime[acttaken]);

  // multiply by discount factor
  temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);
  
  // add reward - should be action dependent
  temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,temp2);

  // now temp1 is the new Q computed for *all* states for acttaken
  // now we need to replace the Q values in the real Q (the one passed in)
  // with these new values - only for the following variable assignments:
  // 1 - varvals
  // 2-k+1 - k other randomly chosen states
  // first - with varvals
  replaceVal(&Q[acttaken],temp1,varvals);

  int *otvv = new int[numorigvars];
  for (j=0; j < k; j++) {
    acttaken = (int) floor(((double) rand())/((double) RAND_MAX+1.0)*numactions);
    for (k=0; k < numorigvars; k++)
      otvv[k] = (int) floor(((double) rand())/((double) RAND_MAX+1.0)*orig_vars[k].nvals);
    // construct variable assignments for this value of orig_vars
    replaceVal(&Q[acttaken],temp1,otvv);
  }
  
  delete otvv;
}
// replaces dd1 by ddr for assignment of variables in varvals
// result it properly reffed
void MDP::replaceVal(DdNode **dd1, DdNode *ddr, int *varvals)
{
  DdNode *temp1, *temp2, *temp3;
  // construct a cube of varvals
  temp1 = buildCubeOrig(varvals, false);
  
  // mutiply this by ddr to get cube with ddr value at leaf
  temp2 = Cudd_addApply(gbm,Cudd_addTimes,ddr,temp1);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);

  // get inverse (1-cube) - so its zero everywhere where varvals = true
  temp1 = Cudd_addApply(gbm,Cudd_addMinus,One,temp2);
  Cudd_Ref(temp1);

  // multiply this by dd1 to wipe out the values of dd1
  temp3 = Cudd_addApply(gbm,Cudd_addTimes,*dd1,temp1);
  Cudd_Ref(temp3);
  Cudd_RecursiveDeref(gbm,temp1);

  // add the two together to get the result
  temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,temp3);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,temp2);
  Cudd_RecursiveDeref(gbm,temp3);
  Cudd_RecursiveDeref(gbm,*dd1);
  *dd1 = temp1;
}

// extracts the policy  \pi(s) = max_{a \in A} Q(s,A)
// returns the policy (pre-reffed)
// and value in ma
DdNode * MDP::extractQPolicy(DdNode **Q, DdNode **ma)
{
  int i;
  DdNode *act;
  DdNode *temp1;
  *ma = Zero;
  Cudd_Ref(*ma);
  act = Zero;	
  Cudd_Ref(act);

  for (i=0; i<numactions; i++) {
    // maximize 
    temp1 = Cudd_addApply(gbm,Cudd_addMaximum,Q[i],*ma);
    Cudd_Ref(temp1);
    recoverPolicy(Q[i],temp1,*ma,&act,i);
    Cudd_RecursiveDeref(gbm,*ma);
    *ma = temp1;
  }
  return act;
}
// gets action with maximal q value
// at state varvals
int MDP::getMaxQ(DdNode **Q, int *varvals)
{
  // check each Vcurrent and get one with maximum value
  int acttotake(-1);
  int i;

  Pair aval;
  double bestactsofar = -10000;
  for (i=0; i<numactions; i++) {
    getVal(gbm, Q[i], aval, varvals);
    if (aval.get_max() > bestactsofar) {
      bestactsofar = aval.get_max();
      acttotake = i;
    }
  }
  if (acttotake < 0) {
    acttotake = (int) floor(((double) rand())/((double) RAND_MAX+1.0)*numactions);
  }
  return acttotake;
}
// performs reinforcement learning loop for numiterations 
// using simulator as the environmental simulator and 
// reward function
// returns reward acheived - 0 if no reward was achieved
// (due to error, or a really shitty agent)
int MDP::reinforcementLearn(MDP *simulator, int numiterations)
{
  int i,j;
  int bigadd;
  int approxMeth, optimalFrom,limit_type,reorderMeth,reorderApplyMeth,shuffleFlag,dotFlag;

  // default constants
  tolerance = 0.1;
  bigadd = BIGADD;
  double Max_Error = MAX_ERROR;
  limit_type = LIMIT_SIZE;    /* can be LIMIT_SIZE or LIMIT_ERROR */
  reorderApplyMeth = REORDERAPPLY_ALL;        /* can be REORDERAPPLY_(FIRSTFIVE | ALL)*/
  reorderMeth = REORDER_NONE;        /* can be REORDER_(SIFT | RANDOM |MINSPAN | NONE)*/
  approxMeth = APPROX_NONE;
  optimalFrom = OPTIMAL_NONE;
  int acttaken;
  shuffleFlag = 0;
  dotFlag = 0;

  // second argument can be the -d <dual optimal> to read in (simulate)
  // or -l <log file> for log read in (no simulation)
  // otherwise, will generate and simulate

  // to do the learning, we must classify the support of each NewPrime diagram
  // and then build arrays of the unprimed variables in the support of each primed diagram
  // this must be done prior to reinforcement Learning (at beginning of trial)
  // need this for the simulator?
  //   simulator->classifySupport();
  classifySupport();
  
  // check if simulator agrees on this, at least
  if (numorigvars != simulator->getNumOrigvars() || numactions != simulator->getNumActions()) {
    fprintf(stderr,"wrong simulator");
    return 0;
  }
  // Vcurrent are the Q functions
  DdNode **Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));

  // build MDP's VCurrentMin - used for storing Q functions
  VcurrentMin = (DdNode **)malloc(numactions*(sizeof(DdNode*)));
  
  // RewardD is the reward function
  // NewPrimes are the transition functions
  DdNode ***counts = new DdNode **[numactions];
  DdNode **rewcounts = new DdNode *[numactions];
  DdNode *tmp, *ttmp;
  for (i=0; i<numactions; i++) {
    counts[i] = new DdNode*[numorigvars];
    rewcounts[i] = Zero;
    Cudd_Ref(rewcounts[i]);
    VcurrentMin[i] = Zero;
    Cudd_Ref(VcurrentMin[i]);
    Vcurrent[i] = Zero;
    Cudd_Ref(Vcurrent[i]);
    for (j=0; j<numorigvars; j++) {
      counts[i][j] = Zero;
      Cudd_Ref(counts[i][j]);
    }
  }
  
  // numorigvars is really the number of original vars
  int *varvals = new int[numorigvars];
  int *newvarvals = new int[numorigvars];
  int *exclude = new int[numorigvars];
  // initial state
  for (i=0; i<numorigvars; i++) {
    varvals[i] = 1;
    newvarvals[i] = 1;
    exclude[i] = 0;
  }
  // for coffee2.cost.dat only
  varvals[1] = 2;
  varvals[3] = 2;
  varvals[5] = 0;
  varvals[6] = 2;
  acttaken = 0;

  double instantR;
  Pair ir(instantR);
  DdNode *instantRConst = Cudd_addConst(gbm,&ir);
  Cudd_Ref(instantRConst);
  double beta(1.0); //(learning rate)
  int numiters(0);
  while (numiters < numiterations) { 
    fprintf(stderr,"%d %f ",numiters,beta);
    
    // simulate new state based on actual model
    instantR = simulator->simulater(varvals,acttaken,newvarvals,exclude);
    ir.set(instantR);
    Cudd_RecursiveDeref(gbm,instantRConst);
    instantRConst = Cudd_addConst(gbm,&ir);
    Cudd_Ref(instantRConst);
    
    /*
    fprintf(stderr,"varvals   : ");
    for (j=0; j<numorigvars; j++)
      fprintf(stderr,"%d ",varvals[j]);

    fprintf(stderr,"\nnewvarvals: ");
    for (j=0; j<numorigvars; j++)
      fprintf(stderr,"%d ",newvarvals[j]);
    fprintf(stderr,"\nreward: %f\n",instantR);
    */
    // add counts to counts matrix
    for (j=0; j < numorigvars; j++) {
      tmp =buildCPTCount(j,acttaken,newvarvals,varvals);
      ttmp = Cudd_addApply(gbm,Cudd_addPlus,counts[acttaken][j],tmp);
      Cudd_Ref(ttmp);
      Cudd_RecursiveDeref(gbm,counts[acttaken][j]);
      Cudd_RecursiveDeref(gbm,tmp);
      counts[acttaken][j] = ttmp;

      // transform counts matrix ->> transition matrix
      normalizeFunction(counts[acttaken][j], &NewPrime[acttaken][j],j);
    } 
    // also add to reward function for newstate only
    // THIS SHOULD BE for current state = but not for action independent reward
    tmp = buildCubeOrig(newvarvals);
    Cudd_Ref(tmp);
    ttmp = Cudd_addApply(gbm,Cudd_addTimes,instantRConst,tmp);
    Cudd_Ref(ttmp);
    Cudd_RecursiveDeref(gbm,tmp);

    // should be - but now action dependent rewards for now
    //tmp = Cudd_addApply(gbm,Cudd_addPlus,rewcounts[acttaken],ttmp);
    tmp = Cudd_addApply(gbm,Cudd_addPlus,rewcounts[0],ttmp);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(gbm,ttmp);
    //Cudd_RecursiveDeref(gbm,rewcounts[acttaken]);
    //rewcounts[acttaken][j] = tmp;
    Cudd_RecursiveDeref(gbm,rewcounts[0]);
    rewcounts[0] = tmp;
    
    // transform reward matrix ->> reward function
    normalizeFunction(rewcounts[0], &RewardD);
    //Cudd_RecursiveDeref(gbm,RewardD);
    //RewardD = rewcounts[0];
    //Cudd_Ref(RewardD);

    // update Q function for current state
    // and for k other states
    // s' is newvarvals
    // s is varvals
    updateQ(Vcurrent,acttaken,varvals);
   
    // update varvals
    for (j=0; j < numorigvars; j++)
      varvals[j] = newvarvals[j];

    // choose a = max_{a \in A} Q(s,a) 
    // but sometimes, choose randomly
    if (((double) rand())/((double) RAND_MAX+1.0) > (beta = beta*0.9999)) {
      acttaken = getMaxQ(Vcurrent, varvals);
      //acttaken = 0;
      fprintf(stderr,"%d *\n",acttaken);
    } else {
      acttaken = (int) floor(((double) rand())/((double) RAND_MAX+1.0)*numactions);
      fprintf(stderr,"%d -\n",acttaken);
    }
    numiters++;
  }
  // do a bit of approximation on the reward function
  Max_Error = 0.1;
  limit_type = LIMIT_ERROR;
  allPairsApprox(&RewardD);
  // extract policy from Q functions
  OptimalPolicy = extractQPolicy(Vcurrent,&OptimalValue);
  writeOptimalVPDot("learned");
  return 1;
}
void MDP::normalizeFunction(DdNode *counts, DdNode **f, bool primedvars)
{
  int i;
  DdNode *tmp, *ttmp; 
  if (*f != NULL)
    Cudd_RecursiveDeref(gbm,*f);
  // first, normalize the counts
  // build a cube of all variables
  int *varvals = new int[numorigvars];
  for (i=0; i<numorigvars; i++)
    varvals[i] = 0;
  tmp = buildCubeOrig(varvals,primedvars);

  // sum over all values of all variables
  ttmp = Cudd_addExistAbstract(gbm,counts,tmp);
  Cudd_Ref(ttmp);
  Cudd_RecursiveDeref(gbm,tmp);
  tmp = Cudd_addApply(gbm,Cudd_addDivide,counts,ttmp);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,ttmp);
  *f = tmp;
  delete [] varvals;
}
// normalizes function counts by dividing it by itself to yield (normalized) function f
// result is pre-reffed
void MDP::normalizeFunction(DdNode *counts, DdNode **f, int j)
{ 
  DdNode *ttmp; 
  if (*f != NULL)
    Cudd_RecursiveDeref(gbm,*f);
  // first, normalize the counts
  // build a cube of primed original variable j
  *f = buildOneCubeOrig(j,0,true);
  
  // abstract existentially (sums over all values of primed var j)
  ttmp = Cudd_addExistAbstract(gbm,counts,*f);
  Cudd_Ref(ttmp);
  Cudd_RecursiveDeref(gbm,*f);
  *f = Cudd_addApply(gbm,Cudd_addDivide,counts,ttmp);
  Cudd_Ref(*f);
  Cudd_RecursiveDeref(gbm,ttmp);
}
double MDP::simulater(int *state, int actionTaken, int *newstate, int *exclude)
{
  simulate(state,actionTaken,newstate,exclude);
  //compute reward achieved and return it
  Pair rval;
  // ********************************
  // THIS SHOULD BE 'state' here not 'newstate' for action-dependent rewards!!!
  getAV(gbm,RewardD, RewardD, rval, rval,  newstate,  vars, numvars, orig_vars, numorigvars);
  
  return (rval.get_max());
  
}

// simulates one step of the MDP given the optimal action from the policy
// assumes the optimal policy has been generated!
void MDP::simulate(int *state, int *newstate)
{
  simulate(state,consultOptimalPolicy(state),newstate);
}

// simulates one step of the MDP from 'state' given actionTaken
void MDP::simulate(int *state, int actionTaken, int *newstate, int *exclude) 
{
  int i;
  // build a cube of the original variables
  DdNode *tCube = buildCubeOrig(state);
  // restrict each newPrime diagram for actionTaken to be tCube
  DdNode *newVarDist, *temp, *pCube;
  for (i=0; i<numorigvars; i++) {
    if (exclude == NULL || !exclude[i])
      newstate[i] = -1;
  }
  // this has to go backwards since newprime diagrams can include
  // primed variables from later in the order - these must also
  // be restricted to their simulated values so they must be done first
  for (i=numorigvars-1; i>=0; i--) {
    if (exclude != NULL && exclude[i]) {
      // this used to be commented out, right?
      //newstate[i] = state[i];
    } else {
      newVarDist= Cudd_addRestrict(gbm,NewPrime[actionTaken][i],tCube);
      Cudd_Ref(newVarDist);
      
      // also need to restrict to the primed variables that we've already selected
      // in the order
      pCube = buildCubeOrig(newstate,true);
      temp = Cudd_addRestrict(gbm,newVarDist,pCube);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,pCube);
      Cudd_RecursiveDeref(gbm,newVarDist);
      newVarDist = temp;
      
      //selects one branch probabilistically according
      //to the distribution over all possible states (which should sum to 1)
      //and returns the index of that branch in the original variable index ordering
      newstate[i] = selectBranchProb(newVarDist,i);
      
      Cudd_RecursiveDeref(gbm,newVarDist);
    }
  }
}

// builds a CPT count cube - This is an ADD which has a 1 terminal for
// the values of primed variables in opvarvals and unprimed variables in ovarvals
// upon which NewPrime[acttaken][opvar] depends
DdNode *MDP::buildCPTCount(int opvar, int acttaken, int *opvarvals, int *ovarvals)
{
  int j;
  int *tvvp = new int[numorigvars];
  int *tvv = new int[numorigvars];
  DdNode *tmp, *ttmp, *count;
  for (j=0; j<numorigvars; j++) {
    tvv[j] = -1;
    tvvp[j] = -1;
    // if NewPrime[acttaken][opvar] depends on unprimed variable j
    if (ovarsupport[acttaken][opvar][j] == 1 || ovarsupport[acttaken][opvar][j] == 3) 
      tvv[j] = ovarvals[j];
    // if NewPrime[acttaken][opvar] depends on primed variable j
    if (ovarsupport[acttaken][opvar][j] > 1)
      tvvp[j] = opvarvals[j];
  }
  // build cube of unprimed conditioning variables
  tmp = buildCubeOrig(tvv,false);
  //ttmp = buildOneCubeOrig(opvar,opvarval,true);
  // build cube of primed conditioning variables
  ttmp = buildCubeOrig(tvvp,true);

  // mutliply together
  count = Cudd_addApply(gbm,Cudd_addTimes,tmp,ttmp);
  Cudd_Ref(count);
  Cudd_RecursiveDeref(gbm,ttmp);
  Cudd_RecursiveDeref(gbm,tmp);
  delete [] tvv;
  delete [] tvvp;
  return count;
}


// classifies the support of all NewPrime diagrams
// puts the result in ovarsupport
void MDP::classifySupport()
{
  int i,j,k;
  DdNode *oCube,*tf,*tg,*tcomm;

  // ovarsupport[i][j][k] is an index which says
  // whether or not the NewPrime[i][j] diagram depends on variable k
  // if ovarsupport[i][j][k] = 0 then no dependence
  // if                      = 1 then depends on unprimed variable k
  // if                      = 2 then depends on primed variable k
  // if                      = 3 then depends on both primed and unprimed variables k
  ovarsupport = new int**[numactions];
  for (i=0; i<numactions; i++) {
    ovarsupport[i] = new int*[numorigvars];
    for (j=0; j<numorigvars; j++) {
      ovarsupport[i][j] = new int[numorigvars];
      for (k=0; k<numorigvars; k++) 
	ovarsupport[i][j][k] = 0;
    }
  }

  for (i=0; i<numactions; i++) 
    for (j=0; j<numorigvars; j++) 
      for (k=0; k<numorigvars; k++) {

	// classify support in unprimed variables first
	// classify intersection of support between NewPrime and a 
	// cube of unprimemd orig variable k
	oCube = buildOneCubeOrig(k,0,false);
	
	if (!Cudd_ClassifySupport(gbm,NewPrime[i][j],oCube,&tcomm,&tf,&tg)) 
	  fprintf(stderr,"classify support failed\n");
	Cudd_Ref(tcomm); Cudd_Ref(tf); Cudd_Ref(tg);
	// if tcomm is non-zero, then we've got overlap and k is part of NewPrime[0][i]
	if (tcomm != Zero && tcomm != One) 
	  ovarsupport[i][j][k] = 1;
	Cudd_RecursiveDeref(gbm,tcomm);
	Cudd_RecursiveDeref(gbm,tf);
	Cudd_RecursiveDeref(gbm,tg);
	Cudd_RecursiveDeref(gbm,oCube);


	// now classify support in primed variables
	oCube = buildOneCubeOrig(k,0,true);
	
	if (!Cudd_ClassifySupport(gbm,NewPrime[i][j],oCube,&tcomm,&tf,&tg)) 
	  fprintf(stderr,"classify support failed\n");
	Cudd_Ref(tcomm); Cudd_Ref(tf); Cudd_Ref(tg);
	// if tcomm is non-zero, then we've got overlap and k is part of NewPrime[0][i]
	if (tcomm != Zero && tcomm != One) 
	  ovarsupport[i][j][k] += 2;
	Cudd_RecursiveDeref(gbm,tcomm);
	Cudd_RecursiveDeref(gbm,tf);
	Cudd_RecursiveDeref(gbm,tg);
	Cudd_RecursiveDeref(gbm,oCube);
      }
}

int MDP::selectBranchProb(DdNode *newVarDist, int ovar) 
{
  int i;
  DdNode *oCube,*tempM;
  // iterate over all possible values of the original variable ovar
  // build cube for each one over prime_vars 
  // multiply by newVarDist
  // could also use Cudd_Eval here
  // but have to be careful of primed-unprimed variables -
  // this is an easier way
  Pair dval;
  double *vvals = new double[orig_vars[ovar].nvals];
  for (i=0; i<orig_vars[ovar].nvals; i++) {
    oCube = buildOneCubeOrig(ovar,i,true);
    tempM = Cudd_addRestrict(gbm,newVarDist,oCube);  
    Cudd_Ref(tempM);
    dval = *(Cudd_V(tempM));
    Cudd_RecursiveDeref(gbm,tempM);
    Cudd_RecursiveDeref(gbm,oCube);
    vvals[i] = dval.get_min();
  }
  // now vvals[i] has the distribution - so choose one value
  // pick random number for 0-1
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  i=0;
  double sum = 0;
  while (i<orig_vars[ovar].nvals && sum < rnd) {
    sum += vvals[i];
    i++;
  }
  i--;
  return i;
}

// counts the number of *states* spanned by the add != 0
// over the novars original variables in listovars
double MDP::countSpan(DdNode *add, int *listovars, int novars, double spanStates) {
  int j;
  double cnt(0), missingcnt(1);
  DdNode *oCube, *tcomm, *tf, *tg;
  // check for constant first
  if (Cudd_IsConstant(add)) 
    if (add != Zero) {
      return spanStates;
    } else {
      return 0;
    }

  DdNode *radd = add;
  Cudd_Ref(add);
  double newSpanStates;

  // check if orig_vars[listovars[0]] is in support of add
  oCube = buildOneCubeOrig(listovars[0],0,false);
  if (!Cudd_ClassifySupport(gbm,add,oCube,&tcomm,&tf,&tg)) 
    fprintf(stderr,"classify support failed\n");
  Cudd_Ref(tcomm); Cudd_Ref(tf); Cudd_Ref(tg);
  // if tcomm is not constant, then we've got no overlap and listovars[i] is not part of add
  if (tcomm == Zero || tcomm == One) {
    missingcnt *= orig_vars[listovars[0]].nvals;
  } else {
    // compute number of states spanned by remainder
    newSpanStates = spanStates/(orig_vars[listovars[0]].nvals);
    for (j=0; j<orig_vars[listovars[0]].nvals; j++) {
      Cudd_RecursiveDeref(gbm,radd);
      radd = restrictVal(gbm,add,vars,numvars,orig_vars,numorigvars,listovars[0],j);
      Cudd_Ref(radd);
      cnt += countSpan(radd,listovars+1,novars-1,newSpanStates);
    }
  }
  return cnt*missingcnt;
}
// builds a 'same dd' for orig variable ovar - this is 
// the dd which is 1 for ovar = oval && ovar' = oval \forall ovals
// result is pre-reffed
DdNode *MDP::buildSameDD(int ovar) 
{
  DdNode *temp,*temp1, *temp2, *samedd;
  int j;
  samedd = Zero;
  Cudd_Ref(samedd);
  for (j=0; j < orig_vars[ovar].nvals; j++) {
    temp1 = buildOneCubeOrig(ovar,j,false);
    Cudd_Ref(temp1);
    temp2 = buildOneCubeOrig(ovar,j,true);
    Cudd_Ref(temp2);

    temp = Cudd_addApply(gbm,Cudd_addTimes,temp1,temp2);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp1);
    Cudd_RecursiveDeref(gbm,temp2);
    
    temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp,samedd);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,samedd);
    Cudd_RecursiveDeref(gbm,temp);
    
    samedd = temp1;
  }
  return samedd;
}
// builds a cube for the ovarval value of original variable ovar
// for either primed or unprimed variables (as dictated by primedvars)
// result is reffed
DdNode * MDP::buildOneCubeOrig(int ovar, int ovarval, bool primedvars)
{
  int j,nbv,nbvals,tmp;
  // figure out the variable values for ovar=ovarval
  nbv = orig_vars[ovar].nbvars;
  nbvals = int(pow(2.0,nbv));
  tmp = nbvals-ovarval-1;
  int *varass = new int[nbv];
  DdNode **varss = new DdNode*[nbv];
  for (j=nbv-1; j>=0; j--) {
    if (primedvars)
      varss[j] = prime_vars[orig_vars[ovar].var1index+j].add_var;
    else
      varss[j] = vars[orig_vars[ovar].var1index+j].add_var;
    Cudd_Ref(varss[j]);
    varass[j] = tmp%2;
    tmp = tmp/2;
  }
  // get cube for that ovar ovarval
  DdNode *oCube = Cudd_addComputeCube(gbm,varss,varass,nbv);
  Cudd_Ref(oCube);
  delete [] varass;
  for (j=0; j<nbv; j++)
    Cudd_RecursiveDeref(gbm,varss[j]);
  delete [] varss;
  return oCube;
}
// builds a cube for the values of the original variables in varvals
// vvars can be either vars (for unprimed) or prime_vars (for primed)
// if varval[i] < 0, it is not included
// result is pre-reffed
DdNode* MDP::buildCubeOrig(int *varvals, bool primedvars)
{
  int i;
  // cube for each original variable
  DdNode *oCube;
  // cube for all original variables
  DdNode *tCube, *tempM;
  tCube = One;
  Cudd_Ref(tCube);
  for (i=0; i<numorigvars; i++) {
    if (varvals[i] >= 0) {
      // build one cube
      oCube = buildOneCubeOrig(i,varvals[i],primedvars);
      
      // and add to the running product
      tempM = Cudd_addApply(gbm,Cudd_addTimes,tCube,oCube);
      Cudd_Ref(tempM);
      Cudd_RecursiveDeref(gbm,tCube);
      Cudd_RecursiveDeref(gbm,oCube);
      tCube = tempM;
    }
  }
  return tCube;
}
  
void MDP::getVal(DdManager *dd, DdNode *val, Pair & rval, int *varvals, int *varvalsp)
{
  int *varass = new int[2*numvars];
  getVarAss(varass, varvals, vars, numvars, orig_vars, numorigvars);
  getVarAssNoZero(varass, varvalsp, prime_vars, numvars, orig_vars, numorigvars);
  rval = *(Cudd_V(Cudd_Eval(dd,val,varass)));
  delete [] varass;
}
void MDP::getVal(DdManager *dd, DdNode *val, Pair & rval, int *varvals, bool primedvars)
{
  int *varass = new int[2*numvars];
  if (primedvars) {
    getVarAss(varass, varvals, prime_vars, numvars, orig_vars, numorigvars);
  } else {
    getVarAss(varass, varvals, vars, numvars, orig_vars, numorigvars);
  }
  rval = *(Cudd_V(Cudd_Eval(dd,val,varass)));
  delete [] varass;
}
void getVarAss(int *varass, int *varvals, rnum *v, int nvars, onum *ov, int novars) 
{
  int i;
  
  for (i=0; i<nvars*2; i++)
    varass[i] = 0;
  getVarAssNoZero(varass,varvals,v,nvars,ov,novars);
}
void getVarAssNoZero(int *varass, int *varvals, rnum *v, int nvars, onum *ov, int novars) 
{
  int i,j,nbv,tmp,nbvals;

  for (i=0; i<novars; i++) {
    nbv = ov[i].nbvars;
    nbvals = int(pow(2.0,nbv));
    tmp = nbvals-varvals[i]-1;
    for (j=nbv-1; j>=0; j--) {
      varass[Cudd_NodeReadIndex(v[ov[i].var1index+j].add_var)] = tmp%2;
      tmp = tmp/2;
    }
  }
}
void getAV(DdManager *dd, DdNode *val, DdNode *act, Pair & dval, Pair & aval, int *varvals,
	   rnum *v, int nvars, onum *ov, int novars) {
  //fprintf(stderr,"------------ nvars is %d\n",nvars);
  // query the  policy action and value with the current values of varvals
  // build an array varass of 1s and 0s over the variables v
  // with assignments corresponding to the varvals (values for the original variables  ov)
  // then call *(Cudd_V(Cudd_Eval(dd,act, varass)))
  int *varass = new int[2*nvars];
  int i,j,nbv,tmp,nbvals;
  for (i=0; i<nvars*2; i++)
    varass[i] = 0;

  for (i=0; i<novars; i++) {
    nbv = ov[i].nbvars;
    nbvals = int(pow(2.0,nbv));
    tmp = nbvals-varvals[i]-1;
    for (j=nbv-1; j>=0; j--) {
      varass[Cudd_NodeReadIndex(v[ov[i].var1index+j].add_var)] = tmp%2;
      tmp = tmp/2;
    }
  }
  /*  fprintf(stderr,"varass is ");
  for (i=0; i<nvars*2; i++)
    fprintf(stderr,"%d ",varass[i]);
  fprintf(stderr,"\n");
  */
  dval = *(Cudd_V(Cudd_Eval(dd,val,varass)));
  aval = *(Cudd_V(Cudd_Eval(dd,act,varass)));
  
  delete [] varass;
}
bool member(int *array, int arraysize, int index) {
  int i;
  for (i=0; i<arraysize; i++) {
    if (index == array[i])
      return true;
  }
  return false;
}
// picks an action out of those represented by lval
// picks the first action available from lval
int pickAction(double lval) {
  int i,j(0),ilval = ((int) lval);
  while (ilval >= 1) {
    i = ilval%2;
    if (i)
      return j;
    ilval = ilval >> 1;
    j++;
  }
  return -1;
}
int binaryAction(int lval)
{
  double foff = pow(2.0,lval);
  return ((int) foff);
}
void pdd(DdNode *b) 
{
  Cudd_PrintDebug(gbm,b,1,2);
}

void MDP::setOptValue(DdNode *opt) {
  if (OptimalValue) {
    Cudd_RecursiveDeref(gbm, OptimalValue);
  }
  OptimalValue = opt;
  Cudd_Ref(OptimalValue);
}
