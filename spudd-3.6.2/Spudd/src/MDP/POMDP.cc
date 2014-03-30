 #include "POMDP.h"
// global for sampling
double samplesum;
int gotsample;
// for renormalization
int extraR;
// constructor
POMDP::POMDP(char *infile, double badd) :
    MDP(infile,badd) {
  // construct new obsfun
  totalObsFun = (DdNode **) malloc(numactions*sizeof(DdNode *));
  DdNode *temp;
  for (int a=0; a<numactions; a++) {
    totalObsFun[a] = One;
    Cudd_Ref(totalObsFun[a]);
    for (int o=0; o<numorigobs; o++) {
      temp = Cudd_addApply(gbm,Cudd_addTimes,totalObsFun[a],ObsFun[a][o]);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,totalObsFun[a]);
      totalObsFun[a] = temp;
    }
    //fprintf(stderr,"total obs fun %d is \n",a);
    //pdd(totalObsFun[a]);
  }
  
}

void POMDP::getMostLikelyBelief(int *varvals)
{
  int i,j;
  DdNode *temp;
  double bv, maxbv;
  for (i=0; i<numorigvars; i++) {
    maxbv = 0.0;
    for (j=0; j<orig_vars[i].nvals; j++) {
      temp = restrictVal(gbm,beliefState[i],prime_vars,numvars,orig_vars,numorigvars,i,j);
      Cudd_Ref(temp);
      bv = (*(Cudd_V(temp))).get_max();
      Cudd_RecursiveDeref(gbm,temp);
      if (bv > maxbv) {
	maxbv = bv;
	varvals[i] = j;
      }
    }
  }
}
// computes the total number of possible observations
void POMDP::computeTotalNumObs()
{
  int o;
  totalNumObs  = 1;
  for (o=0; o<numorigobs; o++) 
    totalNumObs *= orig_obs[o].nvals;
}
// constructs a new ADD observation variable for the novar original multi-valued variable
void POMDP::newADDObs() {
  char tmp[MAXLEN];
  
  // figure out how many new variables to add
  // I realise this could be done by taking logs and ceil, but this is safer
  int i,bvars(1),bnvals(2);
  //fprintf(stderr,"adding new ADD observation %d with %d values\n",numorigobs,orig_obs[numorigobs].nvals);
  while (bnvals < orig_obs[numorigobs].nvals) {
    bvars++;
    bnvals = bnvals*2;
  }
  orig_obs[numorigobs].nbvars = bvars;
  orig_obs[numorigobs].var1index = numobs;
  for (i=0; i<bvars; i++) {
    obs[numobs].orig_num = numorigobs;
    obs[numobs].number = numobs;

    sprintf(tmp,"%s%d",orig_obs[numorigobs].name,i);

    obs[numobs].name = strdup(tmp);
    obs[numobs].add_var = Cudd_addNewVar(gbm);
    Cudd_Ref(obs[numobs].add_var);
    numobs++;
  }
}
DdNode * POMDP::getBelief(int ovar, double *ovarvals, bool primedvars)
{
  int i;
  DdNode *tmp,*dv, *tmp2, *tmp3;
  Pair pv;
  tmp2 = Zero;
  Cudd_Ref(tmp2);
  for (i=0; i<orig_vars[ovar].nvals; i++) {
    tmp = buildOneCubeOrig(ovar,i,primedvars);
    pv.set(ovarvals[i]);
    dv = Cudd_addConst(gbm,&pv);
    Cudd_Ref(dv);
    tmp3 = Cudd_addApply(gbm,Cudd_addTimes,tmp,dv);
    Cudd_Ref(tmp3);
    Cudd_RecursiveDeref(gbm,tmp);
    tmp = Cudd_addApply(gbm,Cudd_addPlus,tmp2,tmp3);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(gbm,tmp2);
    Cudd_RecursiveDeref(gbm,tmp3);
    tmp2 = tmp;
  }
  return tmp2;
}
int POMDP::generatePolicy(int rMeth, int rAMeth, double tol, int badd, float prun, int max_size, 
			  int aMeth, int ltype, int tMode, int sFlag, bool dFlag) {
  return generatePolicy(false);
}
int POMDP::generatePolicy(bool factored_beliefs) {
  DdNode ***fbeliefs;
  DdNode **beliefs;
  DdNode *tmp,*tmp2,*tmp3;
  int numbeliefs(1);
  // figure out horizon
  horizon = computeHorizon(tolerance);
  int *policy;
  if (factored_beliefs) {
#define BDEBUG
    //#define TIGER2
#ifndef BDEBUG
    fbeliefs = (DdNode ***) malloc(numbeliefs*sizeof(DdNode **));
    fbeliefs[0] = beliefState;
    numbeliefs = expandBeliefs(numbeliefs,&fbeliefs,6,0.01,false);
#else
#ifndef TIGER2
    Pair pv;
    DdNode *dv;
    // shoehorn in some fake beliefs for testing
    numbeliefs = 7;
    fbeliefs = (DdNode ***) malloc(numbeliefs*sizeof(DdNode **));
    for (int i=0; i<numbeliefs; i++) {
      fbeliefs[i] = (DdNode **) malloc(1*sizeof(DdNode *));
      fbeliefs[i][0] = NULL;
    }
    fbeliefs[0] = beliefState;
    double tmpbeliefs[2];
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[1][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[2][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.97; tmpbeliefs[1] = 0.03;
    fbeliefs[3][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.03; tmpbeliefs[1] = 0.97;
    fbeliefs[4][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.995; tmpbeliefs[1] = 0.005;
    fbeliefs[5][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.005; tmpbeliefs[1] = 0.995;
    fbeliefs[6][0] = getBelief(0,tmpbeliefs,true);
#else
    numbeliefs = 17;
    fbeliefs = (DdNode ***) malloc(numbeliefs*sizeof(DdNode **));
    for (int i=0; i<numbeliefs; i++) {
      fbeliefs[i] = (DdNode **) malloc(1*sizeof(DdNode *));
      fbeliefs[i][0] = NULL;
      fbeliefs[i][1] = NULL;
    }
    fbeliefs[0] = beliefState;
    
    double tmpbeliefs[2];
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[1][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[2][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[3][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[4][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[5][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[6][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[7][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[8][0] = getBelief(0,tmpbeliefs,true);
    
    tmpbeliefs[0] = 1.0; tmpbeliefs[1] = 0.0;
    fbeliefs[9][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.0; tmpbeliefs[1] = 1.0;
    fbeliefs[10][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 1.0; tmpbeliefs[1] = 0.0;
    fbeliefs[11][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.0; tmpbeliefs[1] = 1.0;
    fbeliefs[12][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 1.0; tmpbeliefs[1] = 0.0;
    fbeliefs[13][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.0; tmpbeliefs[1] = 1.0;
    fbeliefs[14][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 1.0; tmpbeliefs[1] = 0.0;
    fbeliefs[15][0] = getBelief(0,tmpbeliefs,true);
    tmpbeliefs[0] = 0.0; tmpbeliefs[1] = 1.0;
    fbeliefs[16][0] = getBelief(0,tmpbeliefs,true);

    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[1][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[2][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[3][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[4][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.97; tmpbeliefs[1] = 0.03;
    fbeliefs[5][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.97; tmpbeliefs[1] = 0.03;
    fbeliefs[6][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.03; tmpbeliefs[1] = 0.97;
    fbeliefs[7][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.03; tmpbeliefs[1] = 0.97;
    fbeliefs[8][1] = getBelief(1,tmpbeliefs,true);
    
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[9][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.85; tmpbeliefs[1] = 0.15;
    fbeliefs[10][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[11][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.15; tmpbeliefs[1] = 0.85;
    fbeliefs[12][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.97; tmpbeliefs[1] = 0.03;
    fbeliefs[13][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.97; tmpbeliefs[1] = 0.03;
    fbeliefs[14][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.03; tmpbeliefs[1] = 0.97;
    fbeliefs[15][1] = getBelief(1,tmpbeliefs,true);
    tmpbeliefs[0] = 0.03; tmpbeliefs[1] = 0.97;
    fbeliefs[16][1] = getBelief(1,tmpbeliefs,true);
#endif
    
#endif
    policy = valueIteration((int) horizon, numbeliefs,fbeliefs,factored_beliefs);
    printPolicy(policy,numbeliefs,fbeliefs);
  } else {
    beliefs = (DdNode **) malloc(sizeof(DdNode *));
    beliefs[0] = NULL;
    normalizeFunction(initBelief,beliefs,true);

    numbeliefs = expandBeliefs(numbeliefs,&beliefs,10,0.01,false);
    policy = valueIteration((int) horizon, numbeliefs,&beliefs,factored_beliefs);
    printPolicy(policy,numbeliefs,beliefs);
  }   
}

void POMDP::simulate(double tol, int stages) {
  DdNode *tmp, *ttmp;
  fprintf(stderr,"before normalization \n");
  Cudd_PrintDebug(gbm,initBelief,1,2);

  //normalizeReduceFunction(initBelief,beliefs,2,true);
  tmp = NULL;
  normalizeFunction(initBelief,&tmp,true);
  
  fprintf(stderr,"after normalization %lf\n",  checkNormalization(tmp));

  Cudd_PrintDebug(gbm,tmp,1,2);


  simulate(tol, tmp, stages);
}
// write out the alpha vectors and belief states
void POMDP::printAlphas(DdNode **alphas, int numalphas)
{
  FILE *fp;
  fp = fopen("alphas.dat","w");
  printDdNode(gbm,alphas,numalphas,vars,prime_vars,numvars,orig_vars,numorigvars,fp);
  fclose(fp);
}
void POMDP::printBeliefs(DdNode **beliefSet, int numbeliefs)
{
  FILE *fp;
  fp = fopen("beliefs.dat","w");
  printDdNode(gbm,beliefSet,numbeliefs,vars,prime_vars,numvars,orig_vars,numorigvars,fp);
  fclose(fp);
}
void POMDP::getInitialAlphas(DdNode **alphas, int numalphas, int typ) 
{
  int i,j,a;
  double minval;
  Pair minvalp; 
  DdNode *tmp;

  switch (typ) {
  case 0:
    // using smallest possible
    minval = (*Cudd_V(Cudd_addFindMin(gbm,RewardD))).get_min();
    minval *= 1.0/(1.0-discount_factor);
    minvalp.set(minval);
    fprintf(stderr,"starting from %f\n",minval);
    for (i=0; i<numalphas; i++) {
      alphas[i] = Cudd_addConst(gbm,&minvalp);
      Cudd_Ref(alphas[i]);
    }
    break;
  case 1:
    // using just the reward function
    for (i=0; i<numalphas; i++) {
      alphas[i] = RewardD;
      Cudd_Ref(alphas[i]);
    }
    break;
  case 2:
    // using reward + cost
    numalphas = 0;
    for (a=0; a<numactions; a++) {
      tmp = Cudd_addApply(gbm,Cudd_addPlus,RewardD,actionCost[a]);
      Cudd_Ref(tmp);
      j=0;
      while (j<numalphas && tmp != alphas[j])
	j++;
      if (j>=numalphas) 
	alphas[numalphas++] = tmp;
      else 
	Cudd_RecursiveDeref(gbm,tmp);
    }
    break;
  };
}


// simulate the POMDP from belief state belief
// for stages actions
void POMDP::simulate(double tol, DdNode *belief, int stages) {
  int st,i,a,b,tmp;
  int *optacts = new int[numactions];
  int numoptacts;
  double thresh(2.0);
  // figure out horizon
  horizon = computeHorizon(tolerance);
  int numnewbeliefs, numbeliefs = 0;
  DdNode *newbelief, **beliefSet, **newbeliefSet, **tmpbeliefSet;
  beliefSet = NULL;
  bestactions;
  // initial alpha vectors
  numalphas = 1;
  alphas = (DdNode **) malloc(numalphas*sizeof(DdNode *));
  getInitialAlphas(alphas,numalphas,1);

  // generate the initial policy
  // just for consistencty this needs to be allocated 
  int *policy = new int[1];
  int theb, acttotake;
  for (st=0; st<stages; st++) {
    // if belief is not in  beliefSet
    // add it to the beliefSet and generate a new policy
    // this threshold can be set very highn
    if ((theb = member(belief,beliefSet,numbeliefs,thresh)) == -1) {
      tmpbeliefSet = (DdNode **) malloc(1*sizeof(DdNode *));
      tmpbeliefSet[0] = belief;
      Cudd_Ref(tmpbeliefSet[0]);
      // generate more policy from this new belief set
      // normally use horizon here, but ill just use stages
      numnewbeliefs = expandBeliefs(1,&tmpbeliefSet,10,0.08,false);
      // and add it to the old belief set
      newbeliefSet = (DdNode **) malloc((numnewbeliefs+numbeliefs)*sizeof(DdNode*));
      // copy over old beliefs
      // newbeliefSet takes beliefSet's refs
      for (i=0; i<numbeliefs; i++) 
	newbeliefSet[i] = beliefSet[i];
      // copy in the new beliefs
      // also takes the refs
      while (i<numnewbeliefs+numbeliefs) {
	newbeliefSet[i] = tmpbeliefSet[i-numbeliefs];
	i++;
      }
      numbeliefs = numnewbeliefs+numbeliefs;
      // delete tmp and old beliefs
      free(tmpbeliefSet);
      if (beliefSet != NULL)
	free(beliefSet);
      // beliefSet is re-assigned
      beliefSet = newbeliefSet;
      // the current belief will be at position numbeliefs (first one added)
      theb = numbeliefs-numnewbeliefs;
      //now, generate the new alpha vectors for this new belief set
      bestactions = fastPBVI((int) horizon, numbeliefs, &beliefSet, numalphas, &alphas,false);

      printBeliefs(beliefSet,numbeliefs);
      printAlphas(alphas,numalphas);
    }
    // get the policy for the current belief
    acttotake = getPolicy(numalphas,alphas,bestactions,belief);

    printBelief(belief);
    fprintf(stdout," action to take:\n");
    tmp = acttotake;
    a = 0;
    numoptacts = 0;
    while (tmp >= 1) {
      if (tmp%2) {
	fprintf(stdout," %s",actionlist[a].name);
	optacts[numoptacts] = a;
	numoptacts++;
      }
      a++;
      tmp = tmp/2;
    }
    fprintf(stdout,"\n\n");
    // select an action to take at random from the optimal choices
    acttotake = (int) floor(numoptacts*((double) rand())/((double) RAND_MAX+1.0));
    acttotake = optacts[acttotake];
    //simulate the action from belief state belief
    // normally this will be done by actually taking the action
    // and making an observation (instead of sampling), then
    // updating the belief state given the action + observation pair.
    newbelief = simulateAction(belief, acttotake);
    Cudd_RecursiveDeref(gbm,belief);
    belief = newbelief;
  }
}
// compute the horizon for value iteration as in Pineau & Thrun 
double POMDP::computeHorizon(double tolerance) {
  DdNode *temp1, *temp2, *tmp;
  int a;
  double temps0, temps1;
  if (horizon < 0) {
     for (a=0; a<numactions; a++) {
      tmp = Cudd_addApply(gbm,Cudd_addPlus,RewardD,actionCost[a]);
      Cudd_Ref(tmp);
      temp1 = Cudd_addFindMax(gbm, tmp);
      Cudd_Ref(temp1);
      temp2 = Cudd_addFindMin(gbm, tmp);
      Cudd_Ref(temp2);
      temps1 = ((*(Cudd_V(temp1))).get_max() - (*(Cudd_V(temp2))).get_min());
      Cudd_RecursiveDeref(gbm,temp1);
      Cudd_RecursiveDeref(gbm,temp2);
      Cudd_RecursiveDeref(gbm,tmp);
      if (a==0 || temps1 > temps0) 
	temps0 = temps1;
    }
    horizon = 0;
    while (temps0 >= tolerance) {
      temps0 *= discount_factor;
      horizon++;
    }
  }
  return horizon;
}  
int POMDP::getPolicy(int numalphas, DdNode **alphas, int *bestactions, DdNode *belief)
{
  int besta;
  double val = computeValue(numalphas, &besta, alphas, belief);
  return bestactions[besta];
}
int *POMDP::getPolicy(int numalphas, DdNode **alphas, int *bestactions, int numbeliefs, DdNode **beliefs)
{
  int *policy = new int[numbeliefs];
  for (int b=0; b<numbeliefs; b++) 
    policy[b] = getPolicy(numalphas,alphas,bestactions,beliefs[b]);
  return policy;
}
int POMDP::getPolicy(int numalphas, DdNode **alphas, int *bestactions, DdNode **belief)
{
  int besta;
  double val = computeValue(numalphas, &besta, alphas, belief);
  return bestactions[besta];
}
int *POMDP::getPolicy(int numalphas, DdNode **alphas, int *bestactions, int numbeliefs, DdNode ***beliefs)
{
  int *policy = new int[numbeliefs];
  for (int b=0; b<numbeliefs; b++) 
    policy[b] = getPolicy(numalphas,alphas,bestactions,beliefs[b]);
  return policy;
}
int *POMDP::valueIteration(int horizon, int & numbeliefs, DdNode ***beliefs, bool factbelief)
{
  int i,a,o,b;
  int *policy;
  // initial alpha vector
  numalphas = 1; //numactions;
  alphas = (DdNode **) malloc(numalphas*sizeof(DdNode *));
  getInitialAlphas(alphas,numalphas,1);
  bestactions = fastPBVI(horizon,numbeliefs,beliefs,numalphas,&alphas,factbelief);
  if (factbelief) {
    policy = getPolicy(numalphas, alphas, bestactions, numbeliefs,beliefs);
  } else {
    policy = getPolicy(numalphas, alphas, bestactions, numbeliefs,*beliefs);
  }
  //free(alphas);
  return policy;
}  
bool converged(double *Vb, double *oVb, int numb, int counter, int horizon, double tolerance)
{
  double maxvd(-1000000), vd;
  for (int b=0; b<numb; b++) {
    vd = fabs(Vb[b] - oVb[b]);
    if (vd > maxvd)
      maxvd = vd;
  }
  return (counter > horizon);
  //return (maxvd > 0.0 && maxvd < tolerance);
}

// this is Vlassis & Spaan 2003
// returns a list of the best actions to take for each alpha vector int the optimal value function
// if factbelief = true, then factored beliefs: beliefs[i][j] is the ith belief distribution over the jth original variable
int * POMDP::fastPBVI(int horizon, int & numbeliefs, DdNode ***beliefs, int & numalphas, DdNode ***palphas, bool factbelief)
{
  int i,a,o,b;
  DdNode **alphas = *palphas;
  // the temporary belief set
  DdNode ***tbeliefs;
  if (factbelief) {
    tbeliefs = (DdNode ***) malloc(numbeliefs*sizeof(DdNode **));
  } else {
    tbeliefs = (DdNode ***) malloc(sizeof(DdNode **));
    *tbeliefs = (DdNode **) malloc(numbeliefs*sizeof(DdNode*));
  }
  // original indices in temporary belief set (for value storage)
  int *tbindex = new int[numbeliefs];
  double *Vb = new double[numbeliefs];
  double *oldVb = new double[numbeliefs];
  double *tVb = new double[numbeliefs];
  int bestalpha;
  // value iteration
  int nb, sb, maxi, numnewalphas, oldnumalphas, sbestalpha;
  int newbestalpha;
  double maxgamma, sval, newVb;
  DdNode *sbelief;
  DdNode *** gamma_ao;
  DdNode **newalphas;
  int *bestactions = new int[numalphas];
  int *newbestactions;
  counter = 0;
  // compute value of all beliefs
  
  for (b=0; b<numbeliefs; b++) {
    if (factbelief) {
      Vb[b] = computeValue(numalphas, &bestalpha, alphas, beliefs[b]);
    } else {
      Vb[b] = computeValue(numalphas, &bestalpha, alphas, (*beliefs)[b]);
    }
    fprintf(stderr," %f",Vb[b]);
    oldVb[b] = -100000;
  }
  fprintf(stderr,"\n");

  // main loop commences
  int lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  while (lastiteration < 2){ 
    fprintf(stderr,"\niteration %d/%d numbeliefs %d numalphas %d\n",counter,horizon, numbeliefs, numalphas);
    //printAlphas(numalphas,alphas);
    nb = numbeliefs;
    // back up all alpha vectors
    gamma_ao = step1Backup(numalphas,alphas);

    // values of beliefs - and indices of best alphas
    // alphas are unprimed and beliefs are primed
    // reset the temporary belief set
    for (b=0; b<numbeliefs; b++) {
      if (factbelief) {
	tbeliefs[b] = beliefs[b];
      } else {
	tbeliefs[0][b] = (*beliefs)[b];
	Cudd_Ref(tbeliefs[0][b]);
      }
      tbindex[b] = b;
      tVb[b] = Vb[b];
      oldVb[b] = Vb[b];
      fprintf(stderr," %f",Vb[b]);
    }
    fprintf(stderr,"\n");
    // max numbeliefs new alphas
    newalphas = (DdNode **) malloc(numbeliefs*sizeof(DdNode *));
    newbestactions = new int[numbeliefs];
    numnewalphas = 0;
    while (nb > 0) {
      // sample one belief from tbeliefs
      sb = (int) floor(nb*((double) rand())/((double) RAND_MAX+1.0));
      // back it up - this gets the value in maxgamma
      if (factbelief) {
	newbestactions[numnewalphas] = step23Backup(numalphas, newalphas+numnewalphas, &maxgamma, gamma_ao, tbeliefs[sb]);
      } else {
	newbestactions[numnewalphas] = step23Backup(numalphas, newalphas+numnewalphas, &maxgamma, gamma_ao, tbeliefs[0][sb]);
      }	
      fprintf(stderr,"sampled %d which is at %d : old %f new %f  new alpha with action %d ",sb,tbindex[sb],Vb[tbindex[sb]],maxgamma,newbestactions[numnewalphas]);
      if (Vb[tbindex[sb]] > maxgamma)
	fprintf(stderr,"*********************************\n");
      else
	fprintf(stderr,"\n");
      Vb[tbindex[sb]] = maxgamma;
      //pdd(newalphas[numnewalphas]);
      
      numnewalphas++;
      //printAlphas(numnewalphas,newalphas);
      // remove this value from tbeliefs
      nb--;
      if (factbelief) {
	tbeliefs[sb] = tbeliefs[nb];
      } else {
	Cudd_RecursiveDeref(gbm,tbeliefs[0][sb]);
	tbeliefs[0][sb] = tbeliefs[0][nb];
      }
      // takes tbeliefs[nb]'s ref
      tVb[sb] = tVb[nb];
      tbindex[sb] = tbindex[nb];
      
      // check all remaining belief points
      // only if this is not the last iteration, in which case
      // we want to compute a backup for all belief points
      b=0;
      while (!lastiteration && b < nb) {
	if (factbelief) {
	  newVb = computeValue(numnewalphas, &bestalpha, newalphas, tbeliefs[b]);
	} else {
	  newVb = computeValue(numnewalphas, &bestalpha, newalphas, tbeliefs[0][b]);
	}
	// remove the belief point if its new value is greater
	if (newVb >= tVb[b]) {
	  nb--;
	  if (factbelief) {
	    tbeliefs[b] = tbeliefs[nb];
	  } else {	    
	    Cudd_RecursiveDeref(gbm,tbeliefs[0][b]);
	    tbeliefs[0][b] = tbeliefs[0][nb];
	  }
	  tVb[b] = tVb[nb];
	  Vb[tbindex[b]] = newVb;
	  tbindex[b] = tbindex[nb];
	} else {
	  b++;
	}
      }
    }
    // copy over all new alpha vectors
    for (i=0; i<numalphas; i++) 
      Cudd_RecursiveDeref(gbm,alphas[i]);
    free(alphas);
    free(bestactions);
    oldnumalphas = numalphas;
    numalphas = numnewalphas;

    alphas = (DdNode **) malloc(numalphas*sizeof(DdNode *));
    bestactions = new int[numalphas];
    for (i=0; i<numalphas; i++) {
      // takes the ref from newalphas
      // possibly want to *approximate* here
      alphas[i] = newalphas[i];
      bestactions[i] = newbestactions[i];
    }
    free(newalphas);
    delete [] newbestactions;
    
    counter = counter + 1;


    // delete the gamma_aos
    for (i=0; i<oldnumalphas; i++) {
      for (a=0; a<numactions; a++)
	Cudd_RecursiveDeref(gbm,gamma_ao[i][a]);
      free(gamma_ao[i]);
    }
    free(gamma_ao);
    if (lastiteration) {
      lastiteration++;
    } else {
      lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
    }
  }
  if (factbelief) {
    free(tbeliefs);
  } else {
    free(tbeliefs[0]);
    free(tbeliefs);
  }
  delete [] Vb;
  delete [] tVb;
  delete [] tbindex;
  *palphas = alphas;
  return bestactions;
}
// computes the value of a belief vector (max over alpha vectors)
// returns the value and maxi gets set to the index of the alpha vector that
// was maximal
double POMDP::computeValue(int numalphas, int * maxi, DdNode **alphas, DdNode *belief)
{
  int a,j,i;
  double maxgamma = -10000000.0;
  *maxi = -1;
  DdNode *temp, *temp1, *alphap;
  for (a=0; a<numalphas; a++) {
    //     dp = alpha . b
    // switch alpha to primed
    alphap = Cudd_addSwapVariables(gbm,alphas[a],Array1,Array2,numvars);
    Cudd_Ref(alphap);

    temp = Cudd_addApply(gbm,Cudd_addTimes,alphap,belief);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,alphap);
    for (j=0; j<numorigvars; j++) {
      // this actually sums over all variables, so temp should be a constant at the end
      temp1 = sumOutPrime(temp,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,temp);
      temp = temp1;
    }
    if (!Cudd_IsConstant(temp)) {
      fprintf(stderr,"something went wrong ..\n");
      exit(0);
    }
    if ((*Cudd_V(temp)).get_min() > maxgamma) {
      maxgamma = (*Cudd_V(temp)).get_min();
      *maxi = a;
    }
    Cudd_RecursiveDeref(gbm,temp);
  }
  return maxgamma;
}
// computes the value of a factored belief vector (max over alpha vectors)
// returns the value and maxi gets set to the index of the alpha vector that
// was maximal
double POMDP::computeValue(int numalphas, int * maxi, DdNode **alphas, DdNode **belief)
{
  int a,j,i;
  double tmg, maxgamma = -10000000.0;
  *maxi = -1;
  DdNode *temp, *temp1, *alphap;
  for (a=0; a<numalphas; a++) {
    alphap = computeValue(alphas[a],belief);
    if ((tmg = (*Cudd_V(alphap)).get_min()) > maxgamma) {
      maxgamma = tmg;
      *maxi = a;
    }
    Cudd_RecursiveDeref(gbm,alphap);
  }
  return maxgamma;
}
// compute value for a single alplha vector
// factored belief state
DdNode *POMDP::computeValue(DdNode *alpha, DdNode **belief) 
{
  int j;
  DdNode *alphap, *temp;
  //     dp = alpha . b
  // switch alpha to primed
  alphap = Cudd_addSwapVariables(gbm,alpha,Array1,Array2,numvars);
  Cudd_Ref(alphap);

  for (j=0; j<numorigvars; j++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,alphap,belief[j]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,alphap);
    // this actually sums over all variables, so temp should be a constant at the end
    alphap = sumOutPrime(temp,orig_vars+j,prime_vars);
    Cudd_RecursiveDeref(gbm,temp);
  }
  if (!Cudd_IsConstant(alphap)) {
    fprintf(stderr,"something went wrong ..\n");
    exit(0);
  }
  return alphap;
}
/*
// this version is Pineau 2003 PBVI
int *POMDP::valueIteration(int horizon, int & numbeliefs, DdNode ***beliefs)
{
  int i;
  //Number of iteration to convergence
  counter = 0;

  // get initial alphas
  int numalphas = numbeliefs;
  DdNode **alphas = (DdNode **) malloc(numalphas*sizeof(DdNode *));
  for (i=0; i<numalphas; i++) {
    alphas[i] = RewardD;
    Cudd_Ref(alphas[i]);
  }
  // VALUE ITERATION LOOP COMMENCES
  int *policy = new int[numbeliefs];
  while (counter < horizon){ 
    counter = counter + 1;
    fprintf(stderr,"\niteration %d numbeliefs %d numalphas %d\n",counter,numbeliefs, numalphas);
    delete [] policy;
    policy = new int[numbeliefs];
    numalphas = backupValue(policy,numalphas,&alphas, numbeliefs, *beliefs);
    if (counter < horizon && numbeliefs < MAXNUMBELIEFS) 
      numbeliefs = expandBeliefs(numbeliefs, beliefs,0.01,false);
  }
  printAlphas(numalphas,alphas);
  return policy;
}
*/
// backs up value functions in alpha vectors
// returns the number of generated alpha vectors
int POMDP::backupValue(int *policy, int numalphas, DdNode ***palpha, int numbeliefs, DdNode **pbeliefs)
{
  int a,o,b,i,maxi;
  Pair * alphaTol = new Pair(0.001);
  int dupalpha;
  DdNode *newalpha, **alpha;
  DdNode *** gamma_ao = step1Backup(numalphas,*palpha);
  alpha = *palpha;
  // delete old alphas
  for (i=0; i<numalphas; i++) 
    Cudd_RecursiveDeref(gbm,alpha[i]);
  free(alpha);
  // new alpha list  = empty
  // there will be at most numbeliefs new alpha vectors
  
  *palpha = (DdNode **) malloc(numbeliefs*sizeof(DdNode *));
  alpha = *palpha;
  
  // back up all the beliefs
  int newnumalphas=0;
  double maxgamma;
  for (b=0; b<numbeliefs; b++) {
    maxi = step23Backup(numalphas, &newalpha, &maxgamma, gamma_ao, pbeliefs[b]);

    // check if newalpha is already in the list of alphas
    dupalpha = 0;
    for (i=0; !dupalpha && i<newnumalphas; i++) 
      dupalpha = Cudd_EqualSupNorm(gbm,newalpha,alpha[i],alphaTol,0);

    //   add new alpha  to new alpha list
    if (!dupalpha) {
      alpha[newnumalphas] = newalpha;
      Cudd_Ref(alpha[newnumalphas]);
      newnumalphas++;
    }
    policy[b] = maxi;
  }
  // delete the gamma_aos
  for (i=0; i<numalphas; i++) {
    for (a=0; a<numactions; a++) {
      free(gamma_ao[i][a]);
    }
    free(gamma_ao[i]);
  }
  free(gamma_ao);
  
  return newnumalphas;
}
// performs step 1 of backup - constructs set \gamma^{a,o} expected value 
// for each alpha vector - returns a pointer to this set
// which is a numalphas*numactions*numorigobs
DdNode *** POMDP::step1Backup(int numalphas, DdNode **alpha)
{
  int i,a,o,ov;
  DdNode **alphap, ***gamma_ao;

  alphap = (DdNode **) malloc(numalphas*sizeof(DdNode *));
  gamma_ao = (DdNode ***) malloc(numalphas*sizeof(DdNode **));

  DdNode *temp, *temp1, *temp2;
  // first change all variables to primes in alpha(s) -> alpha(s')
  for (i=0; i<numalphas; i++) {
    alphap[i] = Cudd_addSwapVariables(gbm,alpha[i],Array1,Array2,numvars);
    Cudd_Ref(alphap[i]);
    gamma_ao[i] = (DdNode **) malloc(numactions*sizeof(DdNode *));
  }

  //************  construct alpha_i^{a,o}(s) (step 1)

  // for each element in alpha, alpha(s')
  for (i=0; i<numalphas; i++) {
    for (a=0; a<numactions; a++) {
      /*
      temp = alphap[i];
      Cudd_Ref(temp);
      for (o=0; o<numorigobs; o++) {
	//compute alpha_i^{a,o} <- discount*\sum_{s'} NewPrime(s,a,s')*obsFun(o,s',a)*alpha(s')
	temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp,ObsFun[a][o]);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp);
	temp = temp2;
      }
      */
      temp = Cudd_addApply(gbm,Cudd_addTimes,alphap[i],totalObsFun[a]);
      Cudd_Ref(temp);
      temp1 = multiplySumSet(temp,NewPrime[a]);
      Cudd_RecursiveDeref(gbm,temp);
      temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
      Cudd_Ref(temp2);
      Cudd_RecursiveDeref(gbm,temp1);
      gamma_ao[i][a] = temp2;
    }
  }
  for (i=0; i<numalphas; i++) 
    Cudd_RecursiveDeref(gbm,alphap[i]);
  free(alphap);
  
  return (gamma_ao);
  
}

// given the set of backed up alpha vectors for each action and observation, gamma_ao,
// generates a single new alpha vector by for belief pbelief
// this alpha is copied out in *palpha, and the best action(s) to take
// at this belief point are retured in binary (a 1 at pos i if action i is one of the best)
// maxgamma is the resulting value at the belief point
int POMDP::step23Backup(int numalphas, DdNode **palpha, double * maxgamma, DdNode ***gamma_ao, DdNode *pbelief)
{
  int i,j,a,o,ov,b,tmp,oo;
  Pair * alphaTol = new Pair(0.001);
  int dupalpha;
  // change all primed beliefs to unprimed
  DdNode *belief, *temp1,*temp, *temp2;
  belief = Cudd_addSwapVariables(gbm,pbelief,Array2,Array1,numvars);
  Cudd_Ref(belief);
  
  //*************  Construct gamma_ab (step 2)
  DdNode **gamma_ab = (DdNode **) malloc(numactions*sizeof(DdNode *));

  double maxalpha;
  int maxi;
  DdNode *tempgamma, *temp3;
  for (a=0; a<numactions; a++) {
    gamma_ab[a] = RewardD;
    Cudd_Ref(gamma_ab[a]);
    
    // add the action costs. 
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,gamma_ab[a],actionCost[a]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
    gamma_ab[a] = temp2;

    temp2 = Zero;
    Cudd_Ref(temp2);
    // this could be made more efficient by removing these
    // loops and the restrictVal, then maximizing over the result
    // after summing out variables - but also need to keep track
    // of the values of gamma_ao which contributed to the maximization
    // here we have to sum over all possible combinations of observations, 
    for (o=0; o<totalNumObs; o++) {
      maxalpha = -10000000.0;
      for (i=0; i<numalphas; i++) {
	// restrict gamma_ao[i][a] to this value of o
	tmp = o;
	tempgamma = gamma_ao[i][a];
	Cudd_Ref(tempgamma);
	for (oo=numorigobs-1; oo>=0; oo--) {
	  // value of variable oo for this o
	  ov = tmp%(orig_obs[oo].nvals);
	  tmp = tmp/(orig_obs[oo].nvals);
	  temp = restrictVal(gbm,tempgamma,obs,numobs,orig_obs,numorigobs,oo,ov);
	  Cudd_Ref(temp);
	  Cudd_RecursiveDeref(gbm,tempgamma);
	  tempgamma = temp;
	}
	//compute dotproduct alpha.b --> a scalar since this is \sum{s} alpha*b
	temp = Cudd_addApply(gbm,Cudd_addTimes,tempgamma,belief);
	Cudd_Ref(temp);
	for (j=0; j<numorigvars; j++) {
	  // this actually sums over all (unprimed) variables, so temp should only contain observations at the end
	  temp1 = sumOutPrime(temp,orig_vars+j,vars);
	  Cudd_RecursiveDeref(gbm,temp);
	  temp = temp1;
	}
	if (!Cudd_IsConstant(temp)) {
	  fprintf(stderr,"something went wrong ..\n");
	  exit(0);
	}
	if ((*Cudd_V(temp)).get_min() > maxalpha) {
	  maxalpha = (*Cudd_V(temp)).get_min();
	  temp3 = tempgamma;
	  Cudd_Ref(temp3);
	}
	Cudd_RecursiveDeref(gbm,tempgamma);
	Cudd_RecursiveDeref(gbm,temp);
      }
      temp = Cudd_addApply(gbm,Cudd_addPlus,temp3,temp2);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,temp2);
      Cudd_RecursiveDeref(gbm,temp3);
      temp2 = temp;
    }
    temp = Cudd_addApply(gbm,Cudd_addPlus,temp2,gamma_ab[a]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
    gamma_ab[a] = temp;
  }
  //*************  select best action for this belief point (step 3)
  // new alpha list  = empty
  numalphas = 0;
  *maxgamma = -10000000.0;
  maxi = 0;
  double gammaval;
  int maxfac = 1;
  *palpha = One;
  Cudd_Ref(*palpha);
  for (a=0; a<numactions; a++) {
    //     dp = gamma_b^a . b
    temp = Cudd_addApply(gbm,Cudd_addTimes,gamma_ab[a],belief);
    Cudd_Ref(temp);
    for (j=0; j<numorigvars; j++) {
      // this actually sums over all variables, so temp should be a constant at the end
      temp1 = sumOutPrime(temp,orig_vars+j,vars);
      Cudd_RecursiveDeref(gbm,temp);
      temp = temp1;
    }
    if (!Cudd_IsConstant(temp)) {
      fprintf(stderr,"something went wrong ..\n");
      exit(0);
    }
    gammaval = (*Cudd_V(temp)).get_min();
    if (gammaval >= *maxgamma) {
      if (gammaval == *maxgamma) {
	maxi = maxi+maxfac;
      } else {
	maxi = maxfac;
      }
      *maxgamma = gammaval;
      Cudd_RecursiveDeref(gbm,*palpha);
      *palpha = gamma_ab[a];
      Cudd_Ref(*palpha);
    }
    Cudd_RecursiveDeref(gbm,temp);
    maxfac *= 2;
  }
  for (a=0; a<numactions; a++)
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
  free(gamma_ab);
  
  return maxi;
}

// given the set of backed up alpha vectors for each action and observation, gamma_ao,
// generates a single new alpha vector by for belief pbelief
// this alpha is copied out in *palpha, and the best action(s) to take
// at this belief point are retured in binary (a 1 at pos i if action i is one of the best)
// maxgamma is the resulting value at the belief point
// this version with factored beliefs
int POMDP::step23Backup(int numalphas, DdNode **palpha, double * maxgamma, DdNode ***gamma_ao, DdNode **pbelief)
{
  int i,j,a,o,ov,b,tmp,oo;
  Pair * alphaTol = new Pair(0.001);
  int dupalpha;
  DdNode *temp1,*temp, *temp2;
  
  //*************  Construct gamma_ab (step 2)
  DdNode **gamma_ab = (DdNode **) malloc(numactions*sizeof(DdNode *));

  double maxalpha;
  int maxi;
  DdNode *tempgamma, *temp3;
  for (a=0; a<numactions; a++) {
    gamma_ab[a] = RewardD;
    Cudd_Ref(gamma_ab[a]);
    
    // add the action costs. 
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,gamma_ab[a],actionCost[a]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
    gamma_ab[a] = temp2;

    temp2 = Zero;
    Cudd_Ref(temp2);
    // this could be made more efficient by removing these
    // loops and the restrictVal, then maximizing over the result
    // after summing out variables - but also need to keep track
    // of the values of gamma_ao which contributed to the maximization
    // here we have to sum over all possible combinations of observations, 
    for (o=0; o<totalNumObs; o++) {
      maxalpha = -10000000.0;
      for (i=0; i<numalphas; i++) {
	// restrict gamma_ao[i][a] to this value of o
	tmp = o;
	tempgamma = gamma_ao[i][a];
	Cudd_Ref(tempgamma);
	for (oo=numorigobs-1; oo>=0; oo--) {
	  // value of variable oo for this o
	  ov = tmp%(orig_obs[oo].nvals);
	  tmp = tmp/(orig_obs[oo].nvals);
	  temp = restrictVal(gbm,tempgamma,obs,numobs,orig_obs,numorigobs,oo,ov);
	  Cudd_Ref(temp);
	  Cudd_RecursiveDeref(gbm,tempgamma);
	  tempgamma = temp;
	}
	//compute dotproduct alpha.b --> a scalar since this is \sum{s} alpha*b
	temp = computeValue(tempgamma,pbelief);
	if ((*Cudd_V(temp)).get_min() > maxalpha) {
	  maxalpha = (*Cudd_V(temp)).get_min();
	  temp3 = tempgamma;
	  Cudd_Ref(temp3);
	}
	Cudd_RecursiveDeref(gbm,tempgamma);
	Cudd_RecursiveDeref(gbm,temp);
      }
      temp = Cudd_addApply(gbm,Cudd_addPlus,temp3,temp2);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,temp2);
      Cudd_RecursiveDeref(gbm,temp3);
      temp2 = temp;
    }
    temp = Cudd_addApply(gbm,Cudd_addPlus,temp2,gamma_ab[a]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
    gamma_ab[a] = temp;
  }
  //*************  select best action for this belief point (step 3)
  // new alpha list  = empty
  numalphas = 0;
  *maxgamma = -10000000.0;
  maxi = 0;
  double gammaval;
  int maxfac = 1;
  *palpha = One;
  Cudd_Ref(*palpha);
  for (a=0; a<numactions; a++) {
    temp = computeValue(gamma_ab[a],pbelief);
    gammaval = (*Cudd_V(temp)).get_min();
    if (gammaval >= *maxgamma) {
      if (gammaval == *maxgamma) {
	maxi = maxi+maxfac;
      } else {
	maxi = maxfac;
      }
      *maxgamma = gammaval;
      Cudd_RecursiveDeref(gbm,*palpha);
      *palpha = gamma_ab[a];
      Cudd_Ref(*palpha);
    }
    Cudd_RecursiveDeref(gbm,temp);
    maxfac *= 2;
  }
  for (a=0; a<numactions; a++)
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
  free(gamma_ab);
  return maxi;
}

// expand beliefs over a horizon
// add beliefs only if  they differ (sup norm) by more than mthresh to some other belief
// if allbeliefs = false, only add the best belief for each action at each step
int POMDP::expandBeliefs(int numbeliefs, DdNode ***pbeliefs, double horizon, double mthresh, bool allbeliefs)
{
  int i,j;
  for (i=0; i<horizon; i++) {
    fprintf(stderr,"\niteration %d numbeliefs %d number of gcs so far %d\n",i,numbeliefs,Cudd_ReadGarbageCollections(gbm));
    numbeliefs = expandBeliefs(numbeliefs,pbeliefs,mthresh,allbeliefs);
    for (j=0; j<numbeliefs; j++) 
      fprintf(stderr,"%d ",Cudd_DagSize((*pbeliefs)[j]));
    fprintf(stderr,"\n");
  }
  return numbeliefs;
}

// expands the set of beliefs
// input pbeliefs are over primed variables
// for each action, sample a next belief state, then an observation
// from that, then update the belief based on that observation and action,
// but only take it if it is different from a belief already in the 
// set (sup. norm close by mthresh), if allbeliefs = false,
// only take the furthest belief over all the actions, otherwise, take
// all of them
int  POMDP::expandBeliefs(int numbeliefs, DdNode ***pbeliefs, double mthresh, bool allbeliefs)
{
  int i,a,b,c,o;
  double dist,*dists,maxdist;
  DdNode *temp = One;
  Cudd_Ref(temp);
  DdNode *maxb = One;
  Cudd_Ref(maxb);
  DdNode *bestb;
  DdNode **beliefs = *pbeliefs;
  //  new beliefs list <- empty
  // its at most twice the size of the old list
  
  int numnewbeliefs;
  if (allbeliefs)
    numnewbeliefs = numbeliefs*numactions;
  else
    numnewbeliefs = 2*numbeliefs;
  
  DdNode **newbeliefs = (DdNode **) malloc(numnewbeliefs*sizeof(DdNode *));
  numnewbeliefs = 0;
  // anydiff keeps track of if there are any new beliefs at all
  // diffdist checks an individual new belief for similarity with one 
  // that already exists
  int diffdist, anydiff;

  for (b=0; b<numbeliefs; b++) {
    newbeliefs[numnewbeliefs] = beliefs[b];
    // newbeliefs needs its own ref since we still need beliefs 
    // or do we? 
    Cudd_Ref(newbeliefs[numnewbeliefs]);
    numnewbeliefs++;
  }
  //  for each belief, b, in the input list, 
  for (b=0; b<numbeliefs; b++) {
    maxdist = 0.0;
    bestb = One;
    Cudd_Ref(bestb);
    for (a=0; a<numactions; a++) {
      //   b' <-- simulate one step forward from b on action a
      // newbeliefs is over primed variables
      Cudd_RecursiveDeref(gbm,temp);
      temp = simulateAction(beliefs[b],a);
      //fprintf(stderr,"reclaimed nodes %d\n",cuddGarbageCollect(gbm,1));
      // for each belief, c, in the new beliefs list (which is growing, 
      // but includes the old beliefs)
      diffdist = 1;
      for (c = 0; diffdist && c<numnewbeliefs; c++) {
	dist = supremumNorm(newbeliefs[c],temp);
	if (dist < mthresh) {
	  // already in the beliefs list, throw this one out
	  diffdist = 0;
	} else if (dist > maxdist) {
	  maxdist = dist;
	  Cudd_RecursiveDeref(gbm,maxb);
	  maxb = temp;
	  Cudd_Ref(maxb);
	  diffdist = 2;
	}
      }
      // we found a new one
      if (!allbeliefs && diffdist == 2) {
	Cudd_RecursiveDeref(gbm,bestb);
	bestb = maxb;
	Cudd_Ref(bestb);
      }
      // add it if its not already there
      if (allbeliefs && diffdist) {
	newbeliefs[numnewbeliefs] = temp;
	Cudd_Ref(newbeliefs[numnewbeliefs]);
	numnewbeliefs++;
      }
    }
    //put bestb in new beliefs list
    // if its been found
    if (!allbeliefs && bestb != One) {
      newbeliefs[numnewbeliefs] = bestb;
      Cudd_Ref(newbeliefs[numnewbeliefs]);
      numnewbeliefs++;
    }
    Cudd_RecursiveDeref(gbm,bestb);
  }
  // dereference and delete old beliefs list
  Cudd_RecursiveDeref(gbm,maxb);
  for (b=0; b<numbeliefs; b++) 
    Cudd_RecursiveDeref(gbm,beliefs[b]);
  free(beliefs);
  *pbeliefs = newbeliefs;
  // return new number of beliefs
  return numnewbeliefs;
}

// factored belief state versions
// expand beliefs over a horizon
// add beliefs only if  they differ (sup norm) by more than mthresh to some other belief
// if allbeliefs = false, only add the best belief for each action at each step
int POMDP::expandBeliefs(int numbeliefs, DdNode ****pbeliefs, double horizon, double mthresh, bool allbeliefs)
{
  int i,j;
  for (i=0; i<horizon; i++) {
    fprintf(stderr,"\niteration %d numbeliefs %d number of gcs so far %d\n",i,numbeliefs,Cudd_ReadGarbageCollections(gbm));
    numbeliefs = expandBeliefs(numbeliefs,pbeliefs,mthresh,allbeliefs);
  }
  return numbeliefs;
}

// expands the set of beliefs
// input pbeliefs are over primed variables
// for each action, sample a next belief state, then an observation
// from that, then update the belief based on that observation and action,
// but only take it if it is different from a belief already in the 
// set (sup. norm close by mthresh), if allbeliefs = false,
// only take the furthest belief over all the actions, otherwise, take
// all of them
int  POMDP::expandBeliefs(int numbeliefs, DdNode ****pbeliefs, double mthresh, bool allbeliefs)
{
  int i,j,a,b,c,o;
  double dist,*dists,maxdist;
  DdNode *temp;
  DdNode **maxb, **ftemp, **bestb;
  DdNode ***beliefs = *pbeliefs;
  //  new beliefs list <- empty
  // its at most twice the size of the old list
  
  int numnewbeliefs;
  if (allbeliefs)
    numnewbeliefs = numbeliefs*numactions;
  else
    numnewbeliefs = 2*numbeliefs;
  
  DdNode ***newbeliefs = (DdNode ***) malloc(numnewbeliefs*sizeof(DdNode **));
  numnewbeliefs = 0;
  // anydiff keeps track of if there are any new beliefs at all
  // diffdist checks an individual new belief for similarity with one 
  // that already exists
  int diffdist, anydiff;

  for (b=0; b<numbeliefs; b++) {
    newbeliefs[numnewbeliefs] = beliefs[b];
    numnewbeliefs++;
  }
  //  for each belief, b, in the input list, 
  for (b=0; b<numbeliefs; b++) {
    maxdist = 0.0;
    bestb = NULL;
    for (a=0; a<numactions; a++) {
      //   b' <-- simulate one step forward from b on action a
      // newbeliefs is over primed variables
      temp = simulateAction(beliefs[b],a);
      ftemp = beliefMarginal(temp);
      Cudd_RecursiveDeref(gbm,temp);
      //fprintf(stderr,"reclaimed nodes %d\n",cuddGarbageCollect(gbm,1));
      // for each belief, c, in the new beliefs list (which is growing, 
      // but includes the old beliefs)
      diffdist = 1;
      for (c = 0; diffdist && c<numnewbeliefs; c++) {
	dist = supremumNormMar(newbeliefs[c],ftemp);
	if (dist < mthresh) {
	  // already in the beliefs list, throw this one out
	  diffdist = 0;
	} else if (dist > maxdist) {
	  maxdist = dist;
	  maxb = ftemp;
	  diffdist = 2;
	}
      }
      // we found a new one
      if (!allbeliefs && diffdist == 2) {
	bestb = maxb;
      }
      // add it if its not already there
      if (allbeliefs && diffdist) {
	newbeliefs[numnewbeliefs] = ftemp;
	numnewbeliefs++;
      }
    }
    //put bestb in new beliefs list
    // if its been found
    if (!allbeliefs && bestb != NULL) {
      newbeliefs[numnewbeliefs] = bestb;
      numnewbeliefs++;
    }
  }
  *pbeliefs = (DdNode ***) malloc(numnewbeliefs*sizeof(DdNode **));
  for (c=0; c<numnewbeliefs; c++) 
    pbeliefs[0][c] = newbeliefs[c];
  free(newbeliefs);
  // return new number of beliefs
  return numnewbeliefs;
}

// computes the marginal distributions over b
// returns pointer to an array of numorigvars DdNode *
// one for each marginal belief
DdNode **POMDP::beliefMarginal(DdNode *b) 
{
  int i,j;
  DdNode *temp, *temp1;
  double tmpval;
  DdNode **bm = (DdNode **) malloc(numorigvars*sizeof(DdNode *));
  for (j=0; j<numorigvars; j++) {
    temp = b;
    Cudd_Ref(temp);
    // sum out all primes *except* j
    for (i=0; i<numorigvars; i++) {
      if (i != j) {
	temp1 = sumOutPrime(temp,orig_vars+i,prime_vars);
	Cudd_RecursiveDeref(gbm,temp);
	temp = temp1;
      }
    }
    bm[j] = temp;
  }
  return bm;
}
// simulate from b on action a observation function  - P(O|B) - passed in
DdNode * POMDP::simulateAction(DdNode *bb, int a, DdNode *obsprob)
{
  DdNode *newb, *stsamp, *sampo, *temp1, *temp2;
  int o,i,j;

  // get newb with unprimed variables
  newb = Cudd_addSwapVariables(gbm,bb,Array2,Array1,numvars);
  Cudd_Ref(newb);
  for (j=0; j<numorigvars; j++) {
    temp2 = Cudd_addApply(gbm,Cudd_addTimes,newb,NewPrime[a][j]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp2;
  }
  // **** newb COULD BLOW UP HERE (includes primed and unprimed variables!!! *****
  // sum out unprimed variables
  for (j=0; j<numorigvars; j++) {
    temp2 = sumOutPrime(newb,orig_vars+j,vars);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp2;
  }
  temp1 = Cudd_addApply(gbm,Cudd_addTimes,obsprob,newb);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,newb);
  newb = temp1;
  Cudd_Ref(newb);
  // sum out primed variables
  for (j=0; j<numorigvars; j++) {
    temp2 = sumOutPrime(temp1,orig_vars+j,prime_vars);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp2;
  }
  temp2 = Cudd_addApply(gbm,Cudd_addDivide,newb,temp1);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);
  Cudd_RecursiveDeref(gbm, newb);
  newb = temp2;
  return newb;
}

// simulate from b on action a
DdNode * POMDP::simulateAction(DdNode *bb, int a)
{
  DdNode *newb, *stsamp, *sampo, *obsprod, *temp1, *temp2;
  int o,i,j;
  // first, sample from the belief state
  // this gives a primed sample
  /*
  double cnorm = checkNormalization(bb);
  if (fabs(cnorm-1.0) > 1e-8) {
    fprintf(stderr,"problem with normalization of belief %lf\n",cnorm);
    pdd(bb);
  }
  */

  stsamp = sample_Belief(bb);
  // swap variables (want stsamp to be unprimed)
  temp1 = Cudd_addSwapVariables(gbm,stsamp,Array2,Array1,numvars);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,stsamp);
  stsamp = temp1;
  // now, sample from the transition function for that unprimed state
  newb = One;
  Cudd_Ref(newb);
  for (j=0; j<numorigvars; j++) {
    // restrict NewPrime[a][j] to stsamp
    temp2 = Cudd_addRestrict(gbm,NewPrime[a][j],stsamp);
    Cudd_Ref(temp2);
    // temp2 should be a distribution over original primed varible j
    // keep running product
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,newb,temp2);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp1;
  }
  // now newb is new belief state, b, now
  // sample a state from b - returns a cube with the state
  stsamp = sample_Belief(newb);
  //Cudd_PrintDebug(gbm,newb,1,2);
  //Cudd_PrintDebug(gbm,stsamp,1,2);
  // sample an observation from ObsFun[a] for each o at the state given by stsamp
  // and compute product of all ObsFun for the sampled values of observations
  obsprod = One;
  Cudd_Ref(obsprod);
  for (o=0; o<numorigobs; o++) {
    sampo = Cudd_addApply(gbm,Cudd_addTimes,stsamp,ObsFun[a][o]);
    Cudd_Ref(sampo);
    // but sampo still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampo,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampo);
      sampo = temp1;
    }
    // this should be a distribution over observation o only - sample from it
    // returns the cube representing the sampled observation
    temp1 = sample_Observation(sampo,o);
    Cudd_RecursiveDeref(gbm,sampo);
    sampo = temp1;
    
    // now we want the ObsFun at that value of observation: sampo
    temp2 = Cudd_addRestrict(gbm,ObsFun[a][o],sampo);
    Cudd_Ref(temp2);
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,temp2,obsprod);
    Cudd_Ref(temp1);
    //Cudd_PrintDebug(gbm,temp1,1,2);
    Cudd_RecursiveDeref(gbm,obsprod);
    Cudd_RecursiveDeref(gbm,sampo);
    Cudd_RecursiveDeref(gbm,temp2);
    obsprod = temp1;
  }
  // now, we compute the bayesian update to the belief function
  // given this sampled observation
  //multiply belief by transition function and sum out variables
  Cudd_RecursiveDeref(gbm,newb);
  // get newb with unprimed variables
  newb = Cudd_addSwapVariables(gbm,bb,Array2,Array1,numvars);
  Cudd_Ref(newb);
  for (j=0; j<numorigvars; j++) {
    temp2 = Cudd_addApply(gbm,Cudd_addTimes,newb,NewPrime[a][j]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp2;
  }
  // **** newb COULD BLOW UP HERE (includes primed and unprimed variables!!! *****
  // sum out unprimed variables
  for (j=0; j<numorigvars; j++) {
    temp2 = sumOutPrime(newb,orig_vars+j,vars);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp2;
  }
  temp1 = Cudd_addApply(gbm,Cudd_addTimes,obsprod,newb);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,newb);
  newb = NULL;
  // might actually have to normalize this, too
  //fprintf(stderr,"before normalization\n");
  //Cudd_PrintDebug(gbm,temp1,1,2);

  //normalizeReduceFunction(temp1,&newb,2,true);
  // now, here what I think we want to do is *approximate*,
  // but if we do, how on earth shall we sample from an approximate belief function?
  normalizeFunction(temp1,&newb,true);

  //fprintf(stderr,"after normalization %f\n",checkNormalization(newb));
  //Cudd_PrintDebug(gbm,newb,1,2);

  Cudd_RecursiveDeref(gbm,temp1);
  Cudd_RecursiveDeref(gbm,obsprod);
  Cudd_RecursiveDeref(gbm,stsamp);
  return newb;
}
// simulate from b on action a
// for factored belief state
DdNode * POMDP::simulateAction(DdNode **bb, int a)
{
  DdNode *newb, *stsamp, *sampo, *obsprod, *temp1, *temp2;
  int o,i,j;

  // sample from factored belief state
  stsamp = sample_Belief(bb);

  // swap variables (want stsamp to be unprimed)
  temp1 = Cudd_addSwapVariables(gbm,stsamp,Array2,Array1,numvars);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,stsamp);
  stsamp = temp1;

  // now, sample from the transition function for that unprimed state
  newb = One;
  Cudd_Ref(newb);
  for (j=0; j<numorigvars; j++) {
    // restrict NewPrime[a][j] to stsamp
    temp2 = Cudd_addRestrict(gbm,NewPrime[a][j],stsamp);
    Cudd_Ref(temp2);
    // temp2 should be a distribution over original primed varible j
    // keep running product
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,newb,temp2);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,newb);
    newb = temp1;
  }
  // now newb is new belief state, b, now
  // sample a state from b - returns a cube with the state
  stsamp = sample_Belief(newb);
  //Cudd_PrintDebug(gbm,newb,1,2);
  //Cudd_PrintDebug(gbm,stsamp,1,2);
  // sample an observation from ObsFun[a] for each o at the state given by stsamp
  // and compute product of all ObsFun for the sampled values of observations
  obsprod = One;
  Cudd_Ref(obsprod);
  for (o=0; o<numorigobs; o++) {
    sampo = Cudd_addApply(gbm,Cudd_addTimes,stsamp,ObsFun[a][o]);
    Cudd_Ref(sampo);
    // but sampo still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampo,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampo);
      sampo = temp1;
    }
    // this should be a distribution over observation o only - sample from it
    // returns the cube representing the sampled observation
    temp1 = sample_Observation(sampo,o);
    Cudd_RecursiveDeref(gbm,sampo);
    sampo = temp1;
    
    // now we want the ObsFun at that value of observation: sampo
    temp2 = Cudd_addRestrict(gbm,ObsFun[a][o],sampo);
    Cudd_Ref(temp2);
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,temp2,obsprod);
    Cudd_Ref(temp1);
    //Cudd_PrintDebug(gbm,temp1,1,2);
    Cudd_RecursiveDeref(gbm,obsprod);
    Cudd_RecursiveDeref(gbm,sampo);
    Cudd_RecursiveDeref(gbm,temp2);
    obsprod = temp1;
  }
  // now, we compute the bayesian update to the belief function
  // given this sampled observation
  //multiply belief by transition function and sum out variables
  Cudd_RecursiveDeref(gbm,newb);
  newb = bayesianUpdate(bb,a,obsprod);
  Cudd_RecursiveDeref(gbm,obsprod);
  Cudd_RecursiveDeref(gbm,stsamp);
  return newb;
}
double POMDP::checkNormalization(DdNode *b)
{
  int i;
  Pair *thePair;
  int *thecube;
  DdGen * theGen;
  theGen = Cudd_FirstCube(gbm,b,&thecube,&thePair);
  int gen=1;
  double sum(0.0);
  double factor = 1;
  while (gen) {
    // see how many don't care's there are
    factor = 1;
    for (i=0; i<numvars; i++)
      if (thecube[2*i] == 2) 
	factor *= 2;
    sum += factor*(thePair->get_max());
    gen = Cudd_NextCube(theGen,&thecube,&thePair);
  }
  gen = Cudd_GenFree(theGen);
  return sum;
}
void POMDP::flattenT(double ***T, int & ns, int & na) 
{
  MDP::flattenT(T,ns,na);
}

// sample from factored belief
DdNode * POMDP::sample_Belief(DdNode **b)
{
  int i,j;
  DdNode *theSample = One;
  Cudd_Ref(theSample);
  double tmpval, rnd, sum;
  DdNode *temp1, *temp2;
  temp2 = One;
  Cudd_Ref(temp2);

  for (j=0; j<numorigvars; j++) {
    rnd = ((double) rand())/((double) RAND_MAX+1.0);
    sum = 0;
    for (i=0; sum < rnd && i<orig_vars[j].nvals; i++) {
      Cudd_RecursiveDeref(gbm,temp2);
      temp2 = buildOneCubeOrig(j,i,true);
      temp1 = Cudd_addApply(gbm,Cudd_addTimes,temp2,b[j]);
      Cudd_Ref(temp1);
      sum += (*Cudd_V(Cudd_addFindMax(gbm,temp1))).get_max();
      Cudd_RecursiveDeref(gbm,temp1);
    }
    // temp2 is the sample for variable j
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,theSample,temp2);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,theSample);
    theSample = temp1;
  }
  return theSample;
}
//Bayesian update of belief b given action a -
// obsprod gives the observation function at the value of the
// actual observation
DdNode * POMDP::bayesianUpdate(DdNode **b, int a, DdNode *obsprod)
{
  int j;
  // first we have to compute the whole belief state
  DdNode *belief = One;
  Cudd_Ref(belief);
  DdNode *temp;
  for (j=0; j<numorigvars; j++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,b[j],belief);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,belief);
    belief = temp;
  }
  // think we have to do this
  temp = Cudd_addSwapVariables(gbm,belief,Array2,Array1,numvars);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,belief);
  belief = temp;

  for (j=0; j<numorigvars; j++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,belief,NewPrime[a][j]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,belief);
    belief = temp;
  }
  //****************  belief could blow up here**************

  for (j=0; j<numorigvars; j++) {
    temp = sumOutPrime(belief,orig_vars+j,vars);
    Cudd_RecursiveDeref(gbm,belief);
    belief = temp;
  }
  
  temp = Cudd_addApply(gbm,Cudd_addTimes,obsprod,belief);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,belief);
  belief = NULL;
  
  normalizeFunction(temp,&belief,true);

  return belief;
}
// new, correct, sampling method
// iterates through the cubes
DdNode * POMDP::sample_Belief(DdNode *b)
{
  int i;
  // the sample random value
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  Pair *thePair;
  int *thecube;
  int *thevals = new int[numvars];
  DdGen * theGen;
  theGen = Cudd_FirstCube(gbm,b,&thecube,&thePair);
  int gen=1;
  double sum(0.0);
  double factor = 1;
  while (gen) {
    // see how many don't care's there are
    factor = 1;
    for (i=0; i<numvars; i++)
      if (thecube[2*i] == 2) 
	factor *= 2;
    sum += factor*(thePair->get_max());
    if (sum >= rnd) {
      // copy out the cube values
      for (i=0; i<numvars; i++) 
	thevals[i] = thecube[2*i];
      // free the generator
      gen = Cudd_GenFree(theGen);
    } else {
      gen = Cudd_NextCube(theGen,&thecube,&thePair);
    }
  }
  // now build the cube from the current values in 'thecube'
  DdNode **thevars =(DdNode **) malloc(numvars*sizeof(DdNode *));
  for (i=0; i<numvars; i++) {
    if (thevals[i] == 2) {
      // select randomly for don't cares
      rnd = ((double) rand())/((double) RAND_MAX+1.0);
      if (rnd >= 0.5) {
	thevals[i] = 1;
      } else {
	thevals[i] = 0;
      }
    }
    thevars[i] = prime_vars[i].add_var;
    Cudd_Ref(thevars[i]);
  }
  DdNode *theres = Cudd_addComputeCube(gbm,thevars,thevals,numvars);
  Cudd_Ref(theres);
  for (i=0; i<numvars; i++) 
    Cudd_RecursiveDeref(gbm,thevars[i]);
  free(thevars);
  // ?????????????
  delete [] thevals;
  return theres;

}
DdNode * POMDP::sample_Observation(DdNode *b, int o)
{
  
  int i;
  // the sample random value
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  double sum = 0.0;
  DdNode *temp, *theres;
  // iterate over all values of observation o
  for (i=0; i<orig_obs[o].nvals; i++) {
    sum += evalOrigObs(b,o,i);
    if (sum >= rnd) 
      break;
  }
  //fprintf(stderr,"observation made: %s\n",orig_obs[o].valname[i]);
  return buildOneCubeOrigObs(o,i);
}
// samples from the belief function b
DdNode * POMDP::sampleBelief(DdNode *b)
{
  // this has to be a global variable
  samplesum = 0.0;
  gotsample = 0;
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  Pair prnd;
  prnd.set(rnd);
  DdNode *drnd = Cudd_addConst(gbm,&prnd);
  Cudd_Ref(drnd);
  DdNode *bsamp = Cudd_addApply(gbm,My_addSample,b,drnd);
  Cudd_Ref(bsamp);
  Cudd_RecursiveDeref(gbm,drnd);
  return bsamp;
}
//samples from the ofun over values of origobs o
// ofun should only contain the variables corresponding to origobs o
DdNode * POMDP::sampleObservation(DdNode *ofun, int o)
{
  
}

DdNode *My_addSample(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;
  
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    //add value of F to running sum and compare to G
    double fmx = (*Cudd_V(F)).get_max();
    double gmx = (*Cudd_V(G)).get_max();
    samplesum += fmx;
    if (!gotsample && samplesum > gmx) {
      gotsample = 1;
      return (One);
    } else {
      return (Zero);
    }
  } else {
    return NULL;
  }
}
// compute L2 norm between two belief functions over primed variables
double POMDP::L2dist(DdNode *dd1, DdNode *dd2)
{
  double res(0);
  // subtract the two dds
  DdNode *temp = Cudd_addApply(gbm,Cudd_addMinus,dd1,dd2);
  Cudd_Ref(temp);
  DdNode *temp1 = Cudd_addApply(gbm,Cudd_addTimes,temp,temp);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,temp);
  
  // sum over all leaves
  for (int j=0; j<numorigvars; j++) {
    temp = sumOutPrime(temp1,orig_vars+j,prime_vars);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;
  }
  if (Cudd_IsConstant(temp1))
    res = (*(Cudd_V(temp1))).get_min();
  else {
    fprintf(stderr,"whoops - something wrong!\n");
    return -1;
  }
  return res;
}
// prints out belief state in factored form 
// assumed b is over prime vars
void POMDP::printBelief(DdNode *b)
{
  int i,j;
  DdNode *temp, *temp1;
  double tmpval;
  for (j=0; j<numorigvars; j++) {
    temp = b;
    Cudd_Ref(temp);
    // sum out all primes *except* j
    for (i=0; i<numorigvars; i++) {
      if (i != j) {
	temp1 = sumOutPrime(temp,orig_vars+i,prime_vars);
	Cudd_RecursiveDeref(gbm,temp);
	temp = temp1;
      }
    }
    fprintf(stdout,"%s ",orig_vars[j].name);
    for (i=0; i<orig_vars[j].nvals; i++) {
      temp1 = restrictVal(gbm,temp,prime_vars,numvars,orig_vars,numorigvars,j,i);
      Cudd_Ref(temp1);
      tmpval = 0.01*((int) (100*((*(Cudd_V(temp1))).get_max())));
      fprintf(stdout," %s: %3.2f",orig_vars[j].valname[i],tmpval);
      Cudd_RecursiveDeref(gbm,temp1);
    }
    fprintf(stdout,"\n");
  }
  
}
// print belief wihtout names
void POMDP::printBeliefComparetoActual(FILE *fid, DdNode *b, int *varvals) {
  int i,j;
  DdNode *temp, *temp1;
  double tmpval, sum, maxval, prob_actual;
  int maxi;
  // only look at first three
  for (j=0; j<numorigvars-1; j++) {
    temp = b;
    Cudd_Ref(temp);
    // sum out all primes *except* j
    for (i=0; i<numorigvars; i++) {
      if (i != j) {
	temp1 = sumOutPrime(temp,orig_vars+i,prime_vars);
	Cudd_RecursiveDeref(gbm,temp);
	temp = temp1;
      }
    }
    prob_actual = 0;
    maxval = -1;
    for (i=0; i<orig_vars[j].nvals; i++) {
      temp1 = restrictVal(gbm,temp,prime_vars,numvars,orig_vars,numorigvars,j,i);
      Cudd_Ref(temp1);
      tmpval = 0.01*((int) (100*((*(Cudd_V(temp1))).get_max())));
      fprintf(fid,"%3.2f ",tmpval);
      if (i==varvals[j]) 
	prob_actual = tmpval;
      //if (j==2 && i==1-varvals[j])
      //prob_actual = tmpval;
      if (tmpval > maxval) {
	maxval = tmpval;
	maxi = i;
      }
      Cudd_RecursiveDeref(gbm,temp1);
    }
    fprintf(fid,"%f ",prob_actual);
    if (maxi == varvals[j])
      fprintf(fid,"1 ");
    else 
      fprintf(fid,"0 ");
  }
  fprintf(fid,"\n");
}
void POMDP::printBelief(DdNode **b) {
  int j,i;
  DdNode *temp, *temp1;
  double tmpval;
  for (j=0; j<numorigvars; j++) {
    fprintf(stdout,"%s ",orig_vars[j].name);
    for (i=0; i<orig_vars[j].nvals; i++) {
      temp1 = restrictVal(gbm,b[j],prime_vars,numvars,orig_vars,numorigvars,j,i);
      Cudd_Ref(temp1);
      tmpval = 0.01*((int) (100*((*(Cudd_V(temp1))).get_max())));
      fprintf(stdout," %s: %3.2f",orig_vars[j].valname[i],tmpval);
      Cudd_RecursiveDeref(gbm,temp1);
    }
    fprintf(stdout,"\n");
  }
}
void POMDP::printPolicy(int *policy, int numbeliefs, DdNode ***beliefs)
{
  int b,a;
  int tmp;
  for (b=0; b<numbeliefs; b++) {
    fprintf(stderr,"\nbest actions for belief %d are ",b);
    tmp = policy[b];
    a = 0;
    while (tmp >= 1) {
      if (tmp%2) 
	fprintf(stderr," %s",actionlist[a].name);
      a++;
      tmp = tmp/2;
    }
    fprintf(stderr,"\n");
    printBelief(beliefs[b]);
  }
}
void POMDP::printPolicy(int *policy,int numbeliefs, DdNode **beliefs)
{
  int b,a;
  int tmp;
  for (b=0; b<numbeliefs; b++) {
    fprintf(stderr,"\nbest actions for belief %d are ",b);
    tmp = policy[b];
    a = 0;
    while (tmp >= 1) {
      if (tmp%2) 
	fprintf(stderr," %s",actionlist[a].name);
      a++;
      tmp = tmp/2;
    }
    fprintf(stderr,"\n");
    printBelief(beliefs[b]);
    //Cudd_PrintDebug(gbm,beliefs[b],1,2);
  }
}
void POMDP::printAlphas(int numalphas,  DdNode **alpha)
{ 
  int a;
  for (a=0; a<numalphas; a++) {
    fprintf(stderr,"alpha %d\n",a);
    Cudd_PrintDebug(gbm,alpha[a],1,2);
  }
}
// function to subtract two adds taking absolute value same time
DdNode *My_addAbsMinus(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double fmx = (*Cudd_V(F)).get_max();
    double gmx = (*Cudd_V(G)).get_max();
    Pair res;
    res.set(fabs(fmx-gmx));
    DdNode *dres = Cudd_addConst(gbm,&res);
    return dres;
  }
  return (NULL);
} 
// compute supremum norm between two functions
double supremumNorm(DdNode *dd1, DdNode *dd2)
{
  // subtract the two dds taking absolute values at same time
  DdNode *temp = Cudd_addApply(gbm,My_addAbsMinus,dd1,dd2);
  Cudd_Ref(temp);

  // get maximum
  double maxval = (*Cudd_V(Cudd_addFindMax(gbm,temp))).get_max();
  return maxval;
}
// compute supremum norm between two factored belief functions
double POMDP::supremumNormMar(DdNode **dd1, DdNode **dd2)
{
  // subtract the two dds taking absolute values at same time
  DdNode *temp;
  int j;
  double maxval, val;
  for (j=0; j<numorigvars; j++) {
    val = supremumNorm(dd1[j],dd2[j]);
    if (j==0 || val > maxval) 
      maxval = val;
  }
  return maxval;
}
DdNode* POMDP::buildCubeOrigObs(int *obsvals)
{
  int i,j,nbv,tmp;
  // cube for each original observation
  DdNode *oCube;
  // cube for all original observations
  DdNode *tCube, *tempM;
  tCube = One;
  Cudd_Ref(tCube);
  for (i=0; i<numorigobs; i++) {
    if (obsvals[i] >= 0) {
      // build one cube
      oCube = buildOneCubeOrigObs(i,obsvals[i]);
      
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

DdNode * POMDP::buildOneCubeOrigObs(int ovar, int ovarval)
{
  int i,j,nbv,nbvals,tmp;
  // figure out the variable values for ovar=ovarval
  nbv = orig_obs[ovar].nbvars;
  nbvals = int(pow(2.0,nbv));
  tmp = nbvals-ovarval-1;
  int *varass = new int[nbv];
  DdNode **varss = new DdNode*[nbv];
  for (j=nbv-1; j>=0; j--) {
    varss[j] = obs[orig_obs[ovar].var1index+j].add_var;
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

double POMDP::evalOrigObs(DdNode *theadd, int ovar, int ovarval)
{

  DdNode *temp = restrictVal(gbm,theadd,obs,numobs,orig_obs,numorigobs,ovar,ovarval);
  Cudd_Ref(temp);
  if (!Cudd_IsConstant(temp)) {
    fprintf(stderr,"something wrong in evalOrigObs\n");
    exit(0);
  }
  double res = (*Cudd_V(temp)).get_max();
  Cudd_RecursiveDeref(gbm,temp);
  return res;
}
// this normalizes the function, but also reduces the precision to n decimal places
// the result after normalization
void POMDP::normalizeReduceFunction(DdNode *counts, DdNode **f, int n, bool primedvars)
{
  int j;
  DdNode *theSpan, *theNorm, *temp1, *temp2, *tmp;
  // first, renormalize the function
  normalizeFunction(counts,f,primedvars);
  
  // first mulitply by normalizing factor 10^n and floor
  // this actually reduces the precision
  double nf = pow(10.0,n);
  Pair nfp;
  nfp.set(nf);
  DdNode *normfact = Cudd_addConst(gbm,&nfp);
  Cudd_Ref(normfact);
  theNorm = Cudd_addApply(gbm,My_addReducePrecision,*f,normfact);
  Cudd_Ref(theNorm);
  Cudd_RecursiveDeref(gbm,*f);

  // theNorm is now the original belief function 
  // multiplied by 10^n and floored
  // now, we want to compute how much it sums to
  // we want this sum to be 10^n
  temp1 = theNorm;
  Cudd_Ref(temp1);
  // sum over all variables
  for (j=0; j<numorigvars; j++) {
    if (primedvars)
      temp2 = sumOutPrime(temp1,orig_vars+j,prime_vars);
    else
      temp2 = sumOutPrime(temp1,orig_vars+j,vars);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp2;
  }
  // temp1 should be a constant now - get its value
  // this is the extra amount that needs to be added on 
  // in order to make theNorm sum to 10^n
  extraR = (int) (nf-(*(Cudd_V(temp1))).get_max());

  // get the span 
  theSpan =  getSpan(theNorm,extraR,primedvars);
  
  // redistribute the extra   
  //temp2 = Cudd_addApply(gbm,redistributeExtra,theNorm,theSpan);
  temp2 = Cudd_addApply(gbm,Cudd_addPlus,theNorm,theSpan);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,theNorm);
  Cudd_RecursiveDeref(gbm,theSpan);
  theNorm = temp2;

  // and divide by 10^n
  temp2 =  Cudd_addApply(gbm,Cudd_addDivide,theNorm,normfact);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,theNorm);
  theNorm = temp2;

  // check the normalization
  double checkn = checkNormalization(theNorm);
  if (fabs(checkn-1.0) > 1e-12) 
    fprintf(stderr,"whoops! normalization came to %f\n",checkn);


  *f = theNorm;
}
// constructs a new DD with the same structure as f, but 
// with the leaves replaced by the 'span' of the leaves - 
// that is, how many actual states are covered by that leaf
// the states are computed based on primed variables if primedvars!= 0,
DdNode *POMDP::getSpan(DdNode *f, double extra, int primedvars)
{
  int i,j;
  Pair *thePair;
  Pair *fPair = new Pair();
  int *thecube;
  int *thevals = new int[numvars];
  DdGen * theGen;
  DdNode *oCube, *temp1, *factorD;
  DdNode *theSpan = Zero;
  Cudd_Ref(theSpan);
  if (primedvars != 0) 
    primedvars = 1;
  theGen = Cudd_FirstCube(gbm,f,&thecube,&thePair);
  int gen=1;
  double sum;
  double factor = 1;
  int *varass = new int[numvars];
  DdNode **varss = new DdNode*[numvars];
  while (gen && extra > 0) {
    // see how many don't care's there are
    factor = 1;
    j=0;
    for (i=0; i<numvars; i++) {
      if (thecube[2*i+(1-primedvars)] == 2) {
	factor *= 2;
      } else {
	varass[j] = thecube[2*i+(1-primedvars)];
	if (primedvars) 
	  varss[j] = prime_vars[i].add_var;
	else
	  varss[j] = vars[i].add_var;
	j++;
      }
    }
    if (factor > extra) 
      factor = extra;
    extra = extra-factor;
    fPair->set(factor);
    factorD = Cudd_addConst(gbm,fPair);
    Cudd_Ref(factorD);
    oCube = Cudd_addComputeCube(gbm,varss,varass,j);
    Cudd_Ref(oCube);
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,oCube,factorD);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,oCube);
    Cudd_RecursiveDeref(gbm,factorD);
    oCube = temp1;
    temp1 = Cudd_addApply(gbm,Cudd_addPlus,oCube,theSpan);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,oCube);
    Cudd_RecursiveDeref(gbm,theSpan);
    theSpan = temp1;
    
    // get next cube
    gen = Cudd_NextCube(theGen,&thecube,&thePair);
  }
  // free the generator
  Cudd_GenFree(theGen);
  return theSpan;
}
DdNode *My_addReducePrecision(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double fmx = (*Cudd_V(F)).get_max();
    double gmx = (*Cudd_V(G)).get_max();
    Pair res;
    //res.set(ceil(fmx*gmx)/gmx);
    res.set(floor(fmx*gmx));
    DdNode *dres = Cudd_addConst(gbm,&res);
    return dres;
  }
  return (NULL);
} 
DdNode *redistributeExtra(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;

  if (F == DD_MINUS_INFINITY(dd)) return(G);
  if (G == DD_MINUS_INFINITY(dd)) return(F);

  if (F==Zero)
    return (Zero);
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double fmx = (*Cudd_V(F)).get_max();
    double gmx = (*Cudd_V(G)).get_max();
    int igmx = ((int) gmx);
    if (extraR >= igmx) {
      extraR = extraR-igmx;
      // add to f
      Pair res;
      res.set(fmx+1);
      DdNode *dres = Cudd_addConst(gbm,&res);
      return (dres);
    } 
    return (F);
  }
  return (NULL);
} 
// returns the position in beliefSet of the dd closest to belief (sup norm)
// for all dds in beliefSet closer than thresh
// if none are closer than thresh, returns -1
int member(DdNode *belief, DdNode **beliefSet, int numbeliefs, double thresh) {
  int b;
  double dist, mindist;
  int minb(-1);
  mindist = 10.0;
  for (b=0; b<numbeliefs; b++) {
    dist = supremumNorm(beliefSet[b],belief);
    if (dist < thresh) {
      if (dist < mindist) {
	mindist = dist;
	minb = b;
      }
    } 
  }
  return minb;
}
// numdvars is the number of state variables that are printed
// numfvars is the total number of observations
// which would normally be 2*numbehaviors+2
int countLines(char *filename, int numdvars, int numfvars) {
  int j, numlines, act, fno;
  FILE *fid;
  numlines = 0;
  fid = fopen(filename,"r");
  while (EOF != fscanf(fid,"%d ",&fno)) {
    fscanf(fid,"%*d %*d %*d");
    for (j=0; j<numdvars; j++) 
      fscanf(fid,"%*d");
    for (j=0; j<numfvars; j++) 
      fscanf(fid,"%*lf");
    numlines++;
  }
  fclose(fid);
  return numlines;
}

double POMDP::simulateHandwashing(char *filename, char *ofile, int numb, char *lfile) 
{
  int i,j,k,l,t;
  DdNode *obsprob, *tobsprob, *temp, *temp1, *temp2;
  // read in actions,  P(O|B) for data set
  int numdat;
  
  int act;
  FILE *pobf, *bfid;
  int tmp,a,acttotake,numoptacts;
  int *optacts = new int[numactions];
  int hw, hs, wf;
  double *wfp, *pob;
  wfp = new double[2];
  pob = new double[numb];
  // count number of lines
  numdat = countLines(filename,3,numb*2+2);
  pobf = fopen(filename,"r");
  bfid = fopen(ofile,"a");
  // construct Dds for obsprob
  // water flow is ovar[2]  and behavior is ovar[3]
  Pair tmpv;
  int vvals[3] = {2,3,3};
  DdNode *belief = initBelief;
  Cudd_Ref(belief);
  double psum(0.0);
  int frameno, vitdecb, clustb;   
  int *varvals = new int[4];
  Pair rval;
  for (t=0; t<numdat; t++) {
    fscanf(pobf,"%d %d %d",&frameno,&vitdecb,&clustb);
    fscanf(pobf,"%d ",&act);
    // read in state
    fscanf(pobf,"%d %d %d\n",&hw,&hs,&wf);

    // possibly modify the belief here to be exactly this
    // for action verification only
    obsprob = One;
    Cudd_Ref(obsprob);

    for (j=0; j<3; j++) {
      tobsprob = Zero;
      Cudd_Ref(tobsprob);
      for (i=0; i<orig_vars[vvals[j]].nvals; i++) {
	fscanf(pobf,"%lf ",pob+i);
	//pob[i] = 1.0/7.0;
	if (j==0) {
	  pob[i] = 0.5;
	} 
	/*
	else {
	  pob[i] = 1.0/7.0;
	}
	*/
	// use j!=1 here to get normalized probabilities
	// use j!=2 here to get raw likelihoods  (which may be negative)
	// also change below
	if (j != 1) {
	  temp = buildOneCubeOrig(vvals[j],i,true);
	  tmpv.set(pob[i]);
	  temp1 = Cudd_addConst(gbm,&tmpv);
	  Cudd_Ref(temp1);
	  temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp,temp1);
	  Cudd_Ref(temp2);
	  Cudd_RecursiveDeref(gbm,temp);
	  Cudd_RecursiveDeref(gbm,temp1);
	  temp1 = Cudd_addApply(gbm,Cudd_addPlus,tobsprob,temp2);
	  Cudd_Ref(temp1);
	  Cudd_RecursiveDeref(gbm,tobsprob);
	  Cudd_RecursiveDeref(gbm,temp2);
	  tobsprob = temp1;
	}
      }
      if (j != 1) {
	temp1 = Cudd_addApply(gbm,Cudd_addTimes,obsprob,tobsprob);
	Cudd_Ref(temp1);
	Cudd_RecursiveDeref(gbm,obsprob);
	Cudd_RecursiveDeref(gbm,tobsprob);
	obsprob = temp1;
      }
    }
    // compare belief to actual
    // evaluate belief at cube
    varvals[0] = hw; varvals[1]=hs; varvals[2]=wf;
    varvals[3] = orig_vars[vvals[2]].nvals;
    varvals[3] = 0;
    psum = 0.0;
    // evaluate sum by summing over observations
    for (i=0; i<orig_vars[vvals[2]].nvals; i++) {
      getVal(gbm,belief,rval,varvals,true);
      psum  += rval.get_min();
      varvals[3]++;
    }
    fprintf(bfid,"%f ",psum);
    // compare to *old* belief here because
    // we're reading in the value of C (the parents of the behavior)
    // which are S_{t-1} - the previous state
    printBeliefComparetoActual(bfid,belief,varvals);

    // propagate beliefs over data set
    temp2 = simulateAction(belief,act,obsprob);
    Cudd_RecursiveDeref(gbm,belief);
    Cudd_RecursiveDeref(gbm,obsprob);
    belief = temp2;

    // get the policy for the current belief
    //acttotake = getPolicy(numalphas,alphas,bestactions,belief);
    printBelief(belief);
    /*
    fprintf(stdout," action to take:\n");
    tmp = acttotake;
    a = 0;
    numoptacts = 0;
    while (tmp >= 1) {
      if (tmp%2) {
	fprintf(stdout," %s",actionlist[a].name);
	optacts[numoptacts] = a;
	numoptacts++;
      }
      a++;
      tmp = tmp/2;
    }
    fprintf(stdout,"\n\n");
    */
  }
  fclose(pobf);
  // compare the final belief to the final state from the log file
  // this final state doesn't exist in the obsprob file because
  // the behavior is written out with its parents, whcih are S_{t-1}
  if (lfile != NULL) {
    pobf = fopen(lfile,"r");
    
    for (t=0; t<numdat; t++) {
      for (i=0; i<9; i++) {
	fscanf(pobf,"%*d");
      }
    }
    
    fscanf(pobf,"%d %*d %*d",&frameno);
    fscanf(pobf,"%d %d %d\n",&hw,&hs,&wf);
    fclose(pobf);
    
    varvals[0] = hw; varvals[1]=hs; varvals[2]=wf;
    varvals[3] = 0;
    psum = 0.0;
    // evaluate sum by summing over observations
    for (i=0; i<orig_vars[vvals[2]].nvals; i++) {
      getVal(gbm,belief,rval,varvals,true);
      psum  += rval.get_min();
      varvals[3]++;
    }
    fprintf(bfid,"%f ",psum);
    // compare to *old* belief here because
    // we're reading in the value of C (the parents of the behavior)
    // which are S_{t-1} - the previous state
    printBeliefComparetoActual(bfid,belief,varvals);
  }
  fclose(bfid);
  delete [] varvals;

  return psum/numdat;
}
