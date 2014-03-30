#include "MDP.h"
int cntError;
MixGauss avgError;
double tolerance;
MixGauss tempOne,tempTwo;
MixGauss *pZero, *pOne;
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
  int i,j,k;
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
  toleranceP = new MixGauss(*(mdp->toleranceP));
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
    fprintf(stderr,"MDP not allocated properly\n");
    exit(0);
  }
}
void MDP::init() {
  // check if gbm exists already
  if (gbm == NULL)
    gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,MAXMEM);

  pZero = new MixGauss(0.0);
  pOne = new MixGauss(1.0);

  for (int i=0; i<MAXVARS; i++) {
    orig_vars[i].nvals = 0;
    orig_vars[i].nbvars = 0;
  }

  One = Cudd_addConst(gbm,pOne);
  Cudd_Ref(One);

  Zero = Cudd_addConst(gbm,pZero);
  Cudd_Ref(Zero);
  Cudd_AutodynDisable(gbm);

  MixGauss half;
  half.set(0.5);
  Half = Cudd_addConst(gbm,&half);
  Cudd_Ref(Half);

  MixGauss verysmall;
  verysmall.set(-500000000.0);
  VerySmall = Cudd_addConst(gbm,&verysmall);
  Cudd_Ref(VerySmall);

  tolerance = 0.1;
  bigadd = BIGADD;
  
  OptimalPolicy=NULL;
  OptimalValue = NULL;
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
  for(i=0;i<numactions;i++)
    {
      Cudd_RecursiveDeref(gbm,Allprime[i][0]);  
      if (mvinput) {
	for(j=0;j<numorigvars;j++)
	  {
	    Cudd_RecursiveDeref(gbm,NewPrime[i][j]);  
	  }
      } else {
	for(j=0;j<numvars;j++) {
	    Cudd_RecursiveDeref(gbm,Actions[i][j]);
	    Cudd_RecursiveDeref(gbm,NewPrime[i][j]);  
	}
      }
      
    }
  Cudd_RecursiveDeref(gbm,RewardD);
  switch (approxMeth)
    {
    case APPROX_NONE:
      Cudd_RecursiveDeref(gbm,OptimalValue);
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      //Cudd_RecursiveDeref(gbm,ApproximateValue);
      //Cudd_RecursiveDeref(gbm,ApproximatePolicy);
      break;
    };
  Cudd_RecursiveDeref(gbm,Half);
  Cudd_RecursiveDeref(gbm,VerySmall);
  Cudd_RecursiveDeref(gbm,discount);
  Cudd_RecursiveDeref(gbm,discount);
  Cudd_RecursiveDeref(gbm,Zero);
  Cudd_RecursiveDeref(gbm,One);

  
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
    for(i=0;i<numactions;i++)
      free(NewPrime[i]);
    free(NewPrime);
  }
  Cudd_Quit(gbm);
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
rnum *MDP::getPrimeVars() {
  return prime_vars;
}
char *MDP::getActionName(int i) {
  return actionnames[i];
}
DdNode *MDP::getNewPrime(int act, int var) {
  return NewPrime[act][var];
}
DdNode *MDP::getReward() {
  return RewardD;
}
DdNode *MDP::getRewardPrimed() {
  DdNode *temp = Cudd_addSwapVariables(gbm,RewardD,Array1,Array2,numvars);
  Cudd_Ref(temp);
  return temp;
}

void MDP::setNewPrime(int act, int var, DdNode *newNewPrime) {
  Cudd_RecursiveDeref(gbm,NewPrime[act][var]);
  NewPrime[act][var] = newNewPrime;
  Cudd_Ref(NewPrime[act][var]);
}
void decodevvals(int *varvals, onum *orig_vars, int numorigvars, int state)
{
  int i,k;
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
  int i,j,ii,jj,k,l,a;
  DdNode *temp, *temp1, *temp2;
  double pp;
  MixGauss rval;
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
  return !yyparse();
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

  setbuf(stdout,0); 
  

  // check for boolean input file
  mvinput=true;
  fscanf(yyin,"%c",&firstchar);
  if (firstchar == 'v') 
    mvinput = false; 

  // reset the input file
  rewind(yyin);

  if (!parseInput()) 
    return 0;
  

  /* convert the pruning percentage to error by using the span and discount */
  extR = get_extent(RewardD);

  // this is useless here because prune is not set in the input file, 
  // but as an input argument
  switch (toleranceMode) 
    {
    case TOLERANCE_FIXED:
      Max_Error = prune*discount_factor*extR/(1-discount_factor);
      break;
    case TOLERANCE_SLIDING:
      Max_Error = prune*extR;
      break;
    }

  
  /* CREATING THE CONSTANTS FOR THE DISCOUNTING AND THE PRIMING OF THE ACTION */
  fprintf(stderr,"THE DISCOUNT FACTOR IS\n");
  pdd(discount);

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

  toleranceP = new MixGauss(tolerance);

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
  MixGauss maxOPTP, maxPAIR;
  double maxApp,minApp;
  int i;

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
  toleranceP = new MixGauss(tolerance);

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
  DdNode **Vcurrent,*temp3;
  DdNode *Vpast,*temp,*temp1, *temp2;
  DdNode *Merged_Add;
  DdNode *tempA,*tempB;
  
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

int MDP::consultOptimalPolicy(DdNode *pol, DdNode *val, int *varvals)
{
  MixGauss dval, aval;
  getAV(gbm,pol,val, aval,dval, varvals, vars, numvars, orig_vars, numorigvars);
  int bestact;
  if ((bestact = pickAction((int) (aval.get_min()))) < 0) {
    fprintf(stderr,"error in action\n");
    exit(0);
  }
  return bestact;
  
}
int MDP::consultOptimalPolicy(int *varvals) {
  MixGauss dval, aval;
  getAV(gbm,OptimalPolicy, OptimalValue, aval,dval, varvals, vars, numvars, orig_vars, numorigvars);
  int bestact;
  if ((bestact = pickAction((int) (aval.get_min()))) < 0) {
    fprintf(stderr,"error in action\n");
    exit(0);
  }
  return bestact;
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
    case REORDER_EXACT:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_EXACT, 0))
	fprintf(stderr,"*****ERROR: Reorder exact failed\n");
      break;
    case REORDER_NONE:
    default:
      break;
    }


}  
void MDP::allocateMemory() {
  int i,j,k;
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
      allPlist[i][k] = Cudd_NodeReadIndex(prime_vars[orig_vars[j-1].var1index].add_var)+2*orig_vars[j-1].nbvars;
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
  int i,j,k,tempo,tempn,tempfoo,index,nodeCount;
  double temps0,temps1;
  long int usememory=0;
  DdNode *actionDD, *MidPoint_Add;
  DdNode *temp1,*temp2,*tempround;
  DdNode **Vcurrent,*Vpast;
  //DdNode *Reward;

  /* Structure to get the time of execution*/
  struct tms _t;
  long clk_tck = sysconf(_SC_CLK_TCK);
  
  FILE *value,*action,*both,*stats,*log,*log2,*reward,*actions;

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
  
  if (Reward == NULL) {
    Reward = RewardD;
    Cudd_Ref(Reward);
  }
  Merged_Add = Reward;
  Cudd_Ref(Merged_Add);
    
  
  actionDD = Zero;	/* Stupid ref to maintain loop invariant */
  Cudd_Ref(actionDD);

  MidPoint_Add = VerySmall;
  Cudd_Ref(MidPoint_Add);
  
  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  
  lastiteration=0;
  if (horizon == 1)
    lastiteration = 1;
  /* We get the starting up time */
  times(&_t);
  temps0 = 1.0*_t.tms_utime/clk_tck;

  /* VALUE ITERATION LOOP COMMENCES */
  while ((horizon == -1 && (test==0 || lastiteration == 1)) || (horizon >= 0 && counter < horizon)){ 
    

    /* KEEP TRACK OF THE NUMBER OF ITERATIONS */
    counter = counter + 1;
    fprintf(stderr,"\niteration %d\n",counter);
    
    /* update past Value function for convergence check */
    Cudd_RecursiveDeref(gbm,Vpast); 
    Vpast = Merged_Add; 
    Cudd_Ref(Vpast); 
    
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
      
#ifdef DEBUGP
      FILE *fileh = fopen("shit.dot","w");
      DumpDot_p(gbm, temp2, vars, prime_vars, numvars, orig_vars, numorigvars, NULL, fileh);
      fclose(fileh);
#endif

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
	fprintf(stderr,"Memory required exceeds availabe %ld bytes. \nTry a smaller bigadd limit constant\n",MAXMEMHARD);
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
  // success
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
  
  char basepath[256];
  char templog[MAXLINE];
  int  dagCount,treeCount,treeLeafCount,dagLeafCount;
  
  fprintf(stats,"\n\n-------------------- Spudd stats  --------------------\n\n");
  fprintf(stats,"Discount factor is : %f \n Tolerance is: %f \n horizon is: %f \n The BIGADD limit is set to: %d \n Target maximum memory use: %d bytes \n Hard Limit on memory use %d bytes\n",discount_factor, tolerance, horizon, bigadd,MAXMEM, MAXMEMHARD);
  
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
  
  int j;
  DdNode *temp, *temp1, *cube;
  int first, last;

  temp1 = Cudd_addApply(gbm,Cudd_addTimes,prime,reward);
  Cudd_Ref(temp1);

  if (mvinput) {
    temp = sumOutPrime(temp1,ovar,prime_vars);
  } else {
    cube = Cudd_addIte(gbm,varsum,One,Zero);
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
  DdNode *temp, *temp1, *cube;
  int first, last;

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
      cube = Cudd_addIte(gbm,prime_vars[j].add_var,One,Zero);
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
  int i,j,count,test,dagCount,dagLeafCount,treeCount,treeLeafCount;
  DdNode *temp,*temp1,*tempvar,*prodT, *prodE, *tempT, *tempE, *totT, *totE, *tot, *R;
  DdNode *node;
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
  int i,j,k;
  DdNode *temp,*temp2,*temp3, *X_4, *X_2, *offset;
  MixGauss offMixGauss;
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
  offMixGauss.set(poff);
  offset = Cudd_addConst(gbm,&offMixGauss);  
  Cudd_Ref(offset);
#else
  offMixGauss.set(off+1);
  offset = Cudd_addConst(gbm,&offMixGauss);
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
  DdNode *tempres, *tempres2;
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

int MDP::writeDualOptimal(char *binfilename, DdNode *action, DdNode *value, rnum *varsT,rnum *prime_varsT, onum *orig_varsT)
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
  
  Cudd_RecursiveDeref(gbm,List[0]);
  Cudd_RecursiveDeref(gbm,List[1]);

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
//This function reads from a binary file with the rnum structures and the Dual ADDs
//It needs to do error checking for the case where it compares the variables with what is in the current gbm
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
      Dddmp_cuddHeaderLoad(&ddtype,&nvars,&nsuppvars,&nRoots,&orderedVarNames,&suppVarNames,&varIDs,binfilename,fp);
      
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

// writes the CPTS out in 'flat' format as a matrix
void MDP::writeCPTRaw(char *fileprefix) {
  int i,j,k,l,jj,kk;
  DdNode *pCube, *oCube;
  DdNode *temp, *temp1;

  // size of the matrix
  int msze  = 1;
  for (j=0; j<numorigvars; j++) 
    msze *= orig_vars[j].nvals;
  double **mcpt = new double*[msze];
  for (j=0; j<msze; j++) 
    mcpt[j] = new double[msze];
  // variable values
  int *varvals = new int[numorigvars];
  // primed variable values
  int *pvarvals = new int[numorigvars];
  MixGauss theval;
  for(i=0;i<numactions;i++) {
    for (j=0; j<numorigvars; j++) 
      varvals[j] = 0;
    j=0; 
    while (j<msze) {
      for (k=0; k<numorigvars; k++) 
	pvarvals[k] = 0;
      k=0;
      // the unprimed cube
      //oCube = buildCubeOrig(varvals[k],false);
      while (k<msze) {
	// get j,k th matrix entry
	mcpt[j][k] = 1.0;
	for (kk=0; kk<numorigvars; kk++) {
	  // restrict to primed cube for this variable
	  //pCube = buildOneCubeOrig(kk,pvarvals[k],true);
	  getVal(gbm,NewPrime[i][kk],theval,varvals,pvarvals);
	  // now, temp should be a constant - keep a running product
	  mcpt[j][k] *= theval.get_val();
	}
	// increase pvarvals by 1
	kk=numorigvars-1;
	while (kk>= 0 && (++pvarvals[kk] >= orig_vars[kk].nvals)) {
	  pvarvals[kk] = 0;
	  kk--;
	}
	k++;
      }
      // increase varvals by 1
      jj=numorigvars-1;
      while (jj>= 0 && (++varvals[jj] >= orig_vars[jj].nvals)) {
	varvals[jj] = 0;
	jj--;
      }
      j++;
    }
    fprintf(stderr,"T: %s\n",actionnames[i]);
    for (j=0; j<msze; j++) {
      for (k=0; k<msze; k++) 
	fprintf(stderr,"%f ",mcpt[j][k]);
      fprintf(stderr,"\n");
    }
  }
  for (j=0; j<msze; j++) 
    delete [] mcpt[j];
  delete [] mcpt;
      
}
void MDP::writeMDPDot(char *fileprefix) {
  int i,j;
  FILE *fileh;
  char fname[256];

  // reward function
  sprintf(fname,"%sreward.dot",fileprefix);
  fileh = fopen(fname,"w");
  DumpDot(gbm, RewardD, vars, numvars, orig_vars, numorigvars, NULL, fileh);
  fclose(fileh);

  //Dual action diagrams
  for(i=0;i<numactions;i++) {
    for (j=0; j<numorigvars; j++) {
      sprintf(fname,"%sCPT%d-%d.dot",fileprefix,i,j);
      fileh = fopen(fname,"w");
      DumpDot_p(gbm, NewPrime[i][j], vars, prime_vars, numvars, orig_vars, numorigvars, NULL, fileh);
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
  DumpDot(gbm, valADD, vars, numvars, orig_vars, numorigvars, NULL, value);
  fclose(value);
  
  sprintf(fname,"%spolicy.dot",fileprefix);
  
  action = fopen(fname,"w");
  
  DumpDot(gbm, polADD, vars, numvars, orig_vars, numorigvars, actionnames,action);
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
          MixGauss newMixGauss;
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

          newMixGauss.set(mn);
          DdNode *res = Cudd_addConst(dd,&newMixGauss);
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
DdNode *replaceMixGausss(DdManager *dd, DdNode **f, DdNode **newC)
{
  DdNode *tmp, *res, *F, *ROE;
  double min_span,min_temp;
  
  F = *f; ROE = *newC;
  //Look for a constant that equals 
  //MixGaussOne and MixGaussTwo are the 2 constants that need to be replaced by the newC constant node
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

DdNode *  My_cuddAddApplyRecur(DdManager *dd, DdNode * (*op)(DdManager *, DdNode **, DdNode **), 
			       DdNode * f, DdNode * g)
{
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

    T = My_cuddAddApplyRecur(dd,op,fv,gv);
    if (T == NULL) 
      return(NULL);
    Cudd_Ref(T);

    E = My_cuddAddApplyRecur(dd,op,fvn,gvn);
    if (E == NULL) {
      Cudd_RecursiveDeref(dd,T);
      return(NULL);
    }
    Cudd_Ref(E);

    res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
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
DdNode * My_addApply(DdManager *dd,  DdNode * (*op)(DdManager *, DdNode **, DdNode **),DdNode * f,  DdNode * g)
{
    DdNode *res;
    res = My_cuddAddApplyRecur(dd,op,f,g);
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
MixGauss  get_span(DdNode *x) {
  MixGauss *res;
  res = (MixGauss *)malloc(sizeof(MixGauss));
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
  MixGauss error;
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
  DdNode *res;
  DdNode *F, *G;
  double gval,fval;
  double value;
  MixGauss pVal;
  
  
  F = *f; G = *g;
  return F;
}
DdNode * getAction(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *G;
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
  DdNode *F, *G, *res;
  double gval, fval;
  MixGauss pVal;
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
  DdNode *F, *G, *res;
  double gval, fval;
  int ifval;
  MixGauss pVal;
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
  DdNode *res;
  DdNode *F, *G;
  double gval,fval;
  double value;
  MixGauss pVal;
  
  
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
    double value;
    MixGauss pVal;


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
  DdNode *pickAct, *temp, *temp1, *temp2, *sum;
  MixGauss iMixGauss;

  pickAct = Cudd_addConst(gbm,pOne);
  Cudd_Ref(pickAct);
  
  sum = Cudd_addConst(gbm,pZero);
  Cudd_Ref(sum);
  
  temp1 = Cudd_addConst(gbm,pOne);
  Cudd_Ref(temp1);

  for (i=0; i<numact; i++) {
    
    //iMixGauss.set((double) i+1);
    iMixGauss.set((double) i);
    Cudd_RecursiveDeref(gbm,pickAct);
    pickAct = Cudd_addConst(gbm,&iMixGauss);
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
  DdNode *tempT, *tempE, *temp, *cube;
  int i,j,index;
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
DdNode *MDP::sumOutOnePrime(DdNode *dd, int j)
{
  return sumOutPrime(dd,orig_vars+j,prime_vars);
}
// sums out primed variable ovar
DdNode *sumOutPrime(DdNode *dd, onum *ovar, rnum *prime_vars) 
{
  int first, last;
  first = ovar->var1index;
  last = ovar->var1index+ovar->nbvars;      
  return newSumSubtrees(dd,first,last,prime_vars);
}
DdNode *MDP::sumOutAllPrimedOrigVars(DdNode *dd)
{
  int j;
  DdNode *temp, *res;
  res = dd;
  Cudd_Ref(res);
  for (j=0; j<numorigvars; j++) {
    temp = sumOutPrime(res,orig_vars+j,prime_vars);
    Cudd_RecursiveDeref(gbm,res);
    res = temp;
  }
  return res;
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
/*
Sums over the prime_vars in dd from index first to index last-1
for multi-valued diagrams
*/
DdNode *newSumSubtrees(DdNode *dd, int first, int last, rnum *prime_vars)
{
  DdNode *tempT, *tempE, *temp, *cube;
  int i,j,index;
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
  DdNode *tempT, *tempE, *temp, *cube;
  int i,j,index;
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
  DdNode *temp, *sum;
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
  
  MixGauss * pTol = new MixGauss(tol);

  int res = Cudd_EqualSupNorm(dd,val1,val2,pTol,0);
  Cudd_RecursiveDeref(dd,val1);
  Cudd_RecursiveDeref(dd,val2);
  return res;
}
// returns the restricted ADD val to the value of the original variable ovar = ovarval
// doesn't change val
// doesn't ref the return result
DdNode *restrictVal(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		    int ovar, int ovarval) {
  int i,j,nbv,nbvals,tmp;
  // figure out the variable values for ovar=ovarval
  nbv = ov[ovar].nbvars;
  nbvals = int(pow(2.0,nbv));
  tmp = nbvals-ovarval-1;
  int *varass = new int[nbv];
  DdNode **varss = new DdNode*[nbv];
  for (j=nbv-1; j>=0; j--) {
    varss[j] = v[ov[ovar].var1index+j].add_var;
    Cudd_Ref(varss[j]);
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
  DdNode *tmp;
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
    fprintf(fp,"horizon %d\n",horizon);
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
  DdNode *temp1, *temp2;
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
  int i,j;

  MixGauss aval;
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
  int i;
  DdNode *tmp, *ttmp; 
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
  MixGauss rval;
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
  int i,j;
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
  DdNode *snp,*tsnp,*oCube,*tf,*tg,*tcomm;
  int *indexarray = new int[numvars];

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
  MixGauss dval;
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
  double cnt(0), missingcnt(1);
  int i,j,k,l;
  DdNode *oCube, *tcomm, *tf, *tg;
  // check for constant first
  if (Cudd_IsConstant(add)) 
    if (add != Zero) 
      return spanStates;
    else
      return 0;


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
  int i,j;
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
  int i,j,nbv,nbvals,tmp;
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
  int i,j,nbv,tmp;
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
  
void MDP::getVal(DdManager *dd, DdNode *val, MixGauss & rval, int *varvals, int *varvalsp)
{
  int *varass = new int[2*numvars];
  getVarAss(varass, varvals, vars, numvars, orig_vars, numorigvars);
  getVarAssNoZero(varass, varvalsp, prime_vars, numvars, orig_vars, numorigvars);
  rval = *(Cudd_V(Cudd_Eval(dd,val,varass)));
  delete [] varass;
}
void MDP::getVal(DdManager *dd, DdNode *val, MixGauss & rval, int *varvals, bool primedvars)
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
void getAV(DdManager *dd, DdNode *val, DdNode *act, MixGauss & dval, MixGauss & aval, int *varvals,
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
