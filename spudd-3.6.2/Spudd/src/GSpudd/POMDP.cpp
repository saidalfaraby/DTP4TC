 #include "POMDP.h"
// global for sampling
double samplesum;
int gotsample;
// for renormalization
int extraR;
POMDP *__tempMDP;
image obsimage;
Sample **obssamples;
  int totnumsamples;
char ***edgeimage;
int num2Dobsfuns;
int *index2Dobsfuns;
int currBeliefObsImage;

POMDP *__thePOMDP;

//#define FUNKYTIGER
#define DOFASTPBVI 1
//#define SAMPLEOBS
//#define SEQUENTIAL
//#define NORMALIZEOBSFUNS 1
// constructor

POMDP::POMDP(char *infile, double badd) :
    MDP(infile,badd) {
  dosampling = false;
  numalphas = 0;
}
POMDP::POMDP(char *infile, bool ds, double badd) : 
    MDP(infile,badd) {
  dosampling = ds;
  numalphas = 0;
}
void POMDP::allocate2Dstructures(int nsamples) 
{
  int a,i,k;
  // allocate structures to keep track of 
  // regions for two dimensional functions
  // and for displaying regions 
  fprintf(stderr,"allocating edgeimage\n");
  obsimage.nx = 160;
  obsimage.ny = 120;
  obsimage.pixels = new int[obsimage.nx*obsimage.ny];
  obsimage.vals = new double[obsimage.nx*obsimage.ny];
  obsimage.action = new int[obsimage.nx*obsimage.ny];
  if (dosampling) {
    totnumsamples = nsamples;  // total for regular grid: obsimage.nx*obsimage.ny
  } else {
    totnumsamples = obsimage.nx*obsimage.ny;
  }
  /*
  fprintf(stderr,"using %d samples\n",totnumsamples);
  obssamples = new Sample*[totnumsamples];
  for (i=0; i<totnumsamples; i++){
    obssamples[i] = new Sample(2);
  }
  */
  // THIS IS problem dependent - get edge images for only those 
  // observations which are type 2 (2D Gaussians)
  edgeimage = new char**[numactions];
  num2Dobsfuns = 0;
  for (i=0; i<numorigobs; i++) 
    if (orig_obs[i].type == 2) 
      num2Dobsfuns++;
  if (num2Dobsfuns > 0) {
    index2Dobsfuns = new int[numorigobs];
    for (a=0; a<numactions; a++) {
      edgeimage[a] = new char*[num2Dobsfuns];
      k = 0;
      for (i=0; i<numorigobs; i++) 
	if (orig_obs[i].type == 2) {
	  //fprintf(stderr,"k is %d i is %d type is %d name is %s\n",k,i,orig_obs[i].type,orig_obs[i].name);
	  edgeimage[a][k] = new char[obsimage.nx*obsimage.ny];
	  getEdgeImage2D(a,i,edgeimage[a][k]);
	  index2Dobsfuns[i] = k;
	  k++;
	}
    }
    // write the edge images
    char fname[256], imname[256];
    for (a=0; a<numactions; a++) {
      // write edge images
      sprintf(fname,"/project/pomdp/edges%d.pgm",a);
      writeEdgeImagePGM(fname,edgeimage[a][index2Dobsfuns[1]]);
    }
    // write edge images superimposed on color image
    sprintf(imname,"/var/tmp/jhoey/papers/pomdp/handpos.PPM");
    for (a=0; a<numactions; a++) {
      sprintf(fname,"/project/pomdp/eim%d.ppm",a);
      writeEdgeImagePPM(fname,imname,edgeimage[a][index2Dobsfuns[1]]);
    }
    for (a=0; a<numactions; a++) {
      sprintf(fname,"/project/pomdp/ereg%d.ppm",a);
      printObsImage(fname,imname,1,a);
    }
  }
#ifdef NORMALIZEOBSFUNS
  // check normalization of observation functions
  MixGauss *region;
  DdNode *tmp1, *tmp2, *ddregion;
  for (i=0; i<totnumsamples; i++) {
    obssamples[i]->index = 0;
    obssamples[i]->z[0] = i%160;
    obssamples[i]->z[1] = i/160;
  }
  for (i=0; i<numorigobs; i++) {
    
    if (orig_obs[i].type == 2) {
      // all the same region
      region = new MixGauss(1);
    }  else if (orig_obs[i].type == 1) {
      region = new MixGauss(2,0);
      // integrate from -inf: +inf
      region->val = 3;
    }
    ddregion = Cudd_addConst(gbm,region);
    Cudd_Ref(ddregion);
    for (a=0; a<numactions; a++) {
      tmp2 = My_addApply(gbm,ddintegrate,totalObsFun[a][i],ddregion);
      Cudd_Ref(tmp2);
      fprintf(stderr,"normalization for observation function %d action %d is \n",i,a);
      pdd(tmp2);
      tmp1 = Cudd_addApply(gbm,Cudd_addDivide,totalObsFun[a][i],tmp2);
      Cudd_Ref(tmp1);
      Cudd_RecursiveDeref(gbm,tmp2);
      Cudd_RecursiveDeref(gbm,totalObsFun[a][i]);
      totalObsFun[a][i] = tmp1;

      // check again
      tmp2 = My_addApply(gbm,ddintegrate,totalObsFun[a][i],ddregion);
      Cudd_Ref(tmp2);
      fprintf(stderr,"AFTER normalizing, normalization for observation function %d action %d  is \n",i,a);
      pdd(tmp2);
    }
  }
#endif
}
DdNode* POMDP::getInitBelief()
{
  return initBelief;
}
DdNode* POMDP::getInitBeliefPrimed()
{
  DdNode *temp = Cudd_addSwapVariables(gbm,initBelief,Array1,Array2,numvars);
  Cudd_Ref(temp);
  return temp;
}

void POMDP::getMostLikelyBelief(int *varvals)
{
  int i,j;
  DdNode *temp;
  double bv, maxbv;
  for (i=0; i<numorigvars; i++) {
    maxbv = 0.0;
    for (j=0; j<orig_vars[i].nvals; j++) {
      temp = restrictVal(gbm,initBeliefState[i],prime_vars,numvars,orig_vars,numorigvars,i,j);
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
// add a constant to the reward function and costs
// to avoid the root-finding problem
void POMDP::addToReward(double addval)
{
  int a;
  MixGauss pv;
  pv.set(addval);
  DdNode *aconst = Cudd_addConst(gbm,&pv);
  Cudd_Ref(aconst);
  DdNode *temp = Cudd_addApply(gbm,Cudd_addPlus,RewardD,aconst);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,RewardD);
  RewardD = temp;
  for (a=0; a<numactions; a++) {
    temp = Cudd_addApply(gbm,Cudd_addPlus,actionCost[a],aconst);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,actionCost[a]);
    actionCost[a] = temp;
  }
  Cudd_RecursiveDeref(gbm,aconst);
}
DdNode * POMDP::getBelief(int ovar, double *ovarvals, bool primedvars)
{
  int i;
  DdNode *tmp,*dv, *tmp2, *tmp3;
  MixGauss pv;
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
int POMDP::generatePolicy(int nsamples, int minnbeliefs, char *outpath, char *beliefFile) {
  int b;
  DdNode *tmp,*tmp0,*tmp1,*tmp2,*tmp3;
  char fname[256];
  //pdd(totalObsFun[0][1]);
  /*
  int ns, na;
  double **T;
  flattenT(&T,ns,na);
  for (int i=0;i<na; i++) {
    fprintf(stderr,"action %d\n",i);
    for (int j=0;j<ns; j++) {
      for (int k=0; k<ns; k++) {
	fprintf(stderr,"%f ",T[i][j*ns+k]);
      }
      fprintf(stderr,"\n");
    }
  }
  */
  allocate2Dstructures(nsamples);

#ifdef USE_MATHEMATICA
  initialize_mlink();
#endif


  //analyze observation function - turn into discrete
  //discretizeObservationFunction(0,1,obsimage);
  // figure out horizon
  horizon = computeHorizon(tolerance);
  int *policy;

  if (beliefFile == NULL) {
    // generate beliefs
    numbeliefs = 1;
    beliefs = new DdNode **[numbeliefs]; 
    beliefs[0] = initBeliefState;
    expandBeliefs(minnbeliefs,0.1,false);
    //write out beliefs
    sprintf(fname,"%s/beliefs.mdp",outpath);
    printBeliefsFile(fname);
  } else {
    //read beliefs
    readBeliefsFile(beliefFile);
    if (dosampling) {
#ifdef SEQUENTIAL
    // double-up the beliefs using the last variable (obs_step) either 1 or 0
    // create 2 cubes for two values of obs_step
    int i,k;
    DdNode **os = new DdNode*[2];
    os[0] = buildOneCubeOrig(numorigvars-1,0,true);
    os[1] = buildOneCubeOrig(numorigvars-1,1,true);
    DdNode ***newbeliefs = new DdNode**[numbeliefs*2];
    for (b=0; b<numbeliefs; b++) {
      newbeliefs[2*b] = new DdNode*[numorigvars];
      newbeliefs[2*b+1] = new DdNode*[numorigvars];
      // copy everything but the obsstep one (which will be wrong)
      // first of the new beliefs takes old beliefs ref
      for (i=0; i<numorigvars-1; i++) {
	for (k=0; k<2; k++) {
	  // beliefs will all be multiplied by 2 since obs_step doesn't exist in belief file
	  
	  tmp = Cudd_addApply(gbm,Cudd_addTimes,beliefs[b][i],Half);
	  Cudd_Ref(tmp);
	  newbeliefs[2*b+k][i] = tmp;
	  //newbeliefs[2*b+k][i] = NULL;
	  //normalizeFunction(tmp,&newbeliefs[2*b+k][i],true);
	  //Cudd_RecursiveDeref(gbm,tmp);
	}
	Cudd_RecursiveDeref(gbm,beliefs[b][i]);
      }
      delete [] beliefs[b];
      // fudge in the obsstep one
      for (k=0; k<2; k++) {
	newbeliefs[2*b+k][numorigvars-1] = os[k];
	Cudd_Ref(newbeliefs[2*b+k][numorigvars-1]);
      }
    }
    Cudd_RecursiveDeref(gbm,os[0]);
    Cudd_RecursiveDeref(gbm,os[1]);
    delete [] os;
    // delete old beliefs
    delete [] beliefs;
    beliefs = newbeliefs;
    numbeliefs *= 2;
#endif
    }
  }

  for (b=0; b<numbeliefs; b++) {
    fprintf(stdout,"belief %d is\n",b);
    printBelief(beliefs[b]);
    fprintf(stdout,"\n");
  }
  double *beliefValues = new double[numbeliefs];
  policy = valueIteration((int) horizon, beliefValues, outpath);
  fprintf(stderr,"%d alpha vectors in total:\n",numalphas);
  printPolicy(policy, numbeliefs, beliefs, beliefValues);

  char htmlfile[256];
  sprintf(htmlfile,"%s/output.html",outpath);
  for (int o=0; o<numorigobs; o++) {
    if (orig_obs[o].type == 2) {
      printPolicyHTML(htmlfile, policy, numbeliefs, beliefs, beliefValues, (int) (horizon+1),o);
    }
  }
#ifdef USE_MATHEMATICA
  terminate_mlink();
#endif
}

// write out the alpha vectors and belief states
void POMDP::printAlphasFile(char *filename)
{
  FILE *fp;
  fp = fopen(filename,"w");
  fprintf(fp,"%d ",numalphas);
  for (int i=0; i<numalphas; i++) {
    fprintf(fp,"%d ",bestactions[i]);
  }
  fprintf(fp,"\n");
  printDdNode(gbm,alphas,numalphas,vars,prime_vars,numvars,orig_vars,numorigvars,fp,"OPTALPHA");
  fclose(fp);
}
void POMDP::readAlphasFile(char *filename)
{
  FILE *fp;
  fp = fopen(filename,"r");
  __thePOMDP = this;
  fscanf(fp,"%d ",&numalphas);
  bestactions = new int[numalphas];
  alphas = new DdNode*[numalphas];
  for (int i=0; i<numalphas; i++) 
    fscanf(fp,"%d ",bestactions+i);
  int ona = numalphas;
  numalphas = 0;
  yyin = fp;
  yyparse();
  fclose(fp);
  if (numalphas != ona) {
    fprintf(stderr,"something wrong - file %s specified %d alpha vectors but there were only %d\n",filename,ona,numalphas);
    exit(0);
  }
}
// write out the alpha vectors and belief states
void POMDP::printBeliefsFile(char *filename)
{
  FILE *fp;
  // first compute ufactored beliefs
  ufbeliefs = new DdNode*[numbeliefs];
  for (int i=0; i<numbeliefs; i++) 
    ufbeliefs[i] = beliefProduct(beliefs[i]);
  fp = fopen(filename,"w");
  fprintf(fp,"%d\n",numbeliefs);
  printDdNode(gbm,ufbeliefs,numbeliefs,vars,prime_vars,numvars,orig_vars,numorigvars,fp,"BELIEFSAMPLE");
  fclose(fp);
}
void POMDP::readBeliefsFile(char *filename)
{
  FILE *fp;
  fp = fopen(filename,"r");
  __thePOMDP = this;
  fscanf(fp,"%d ",&numbeliefs);
  // reads in ufactored beliefs
  ufbeliefs = new DdNode*[numbeliefs];
  yyin = fp;
  int onb = numbeliefs;
  numbeliefs = 0;
  yyparse();
  fclose(fp);
  if (numbeliefs != onb) {
    fprintf(stderr,"something wrong - file %s specified %d alpha vectors but there were only %d\n",filename,onb,numbeliefs);
    exit(0);
  }
  beliefs = new DdNode**[numbeliefs];
  // compute belief marginals
  for (int i=0; i<numbeliefs; i++) 
    beliefs[i] = beliefMarginal(ufbeliefs[i]);
  initBeliefState = beliefs[0];
}
void POMDP::printConditionalPlansForBeliefs(int *pb, int numpb, char *outpath)
{
  int b;
  DdNode *tmp,*tmp2,*tmp3;
  // argument here should not matter
  // done outside now
  //allocate2Dstructures(0);

  // figure out horizon
  horizon = computeHorizon(tolerance);
  int *policy;

  // generate the beliefs
  // read in beliefs
  double * beliefValues = new double[numbeliefs];
  policy = getPolicy(bestactions, numbeliefs, beliefs, beliefValues);

  for (int o=0; o<numorigobs; o++) {
    if (orig_obs[o].type == 2) {
      printConditionalPlanRegions2(outpath,policy,numpb,pb,beliefs,o);
    }
  }
}

void POMDP::printBeliefs(DdNode **beliefSet, int numbeliefs)
{
  FILE *fp;
  fp = fopen("beliefs.dat","w");
  printDdNode(gbm,beliefSet,numbeliefs,vars,prime_vars,numvars,orig_vars,numorigvars,fp);
  fclose(fp);
}
void POMDP::getInitialAlphas(int typ) 
{
  int i,j,a;
  double minval;
  MixGauss minvalp;
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

// compute the horizon for value iteration as in Pineau & Thrun 
double POMDP::computeHorizon(double tolerance) {
  DdNode *temp1, *temp2, *tmp;
  int a;
  double temps1, temps0;
  if (horizon < 0) {
    // Find the largest reward+cost over all actions
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
int *POMDP::getPolicy(int *bestactions, int numbeliefs, DdNode ***beliefs, double *beliefValues)
{
  int *policy = new int[numbeliefs];
  int besta;
  for (int b=0; b<numbeliefs; b++) 
    policy[b] = getPolicy(bestactions,beliefs[b],beliefValues[b]);
  return policy;
}

// for one belief
int POMDP::getPolicy(int *bestactions, DdNode **belief, double & beliefValue)
{
  int besta;
  beliefValue = computeValue(besta, belief); 
  return bestactions[besta];
}
// non-factored beliefs
int POMDP::getPolicy(int *bestactions, DdNode *belief, double & beliefValue)
{
  int besta;
  beliefValue = computeValue(besta, belief); 
  return bestactions[besta];
}
// for one belief
int POMDP::getPolicy(DdNode **belief, double & beliefValue)
{
  int besta;
  beliefValue = computeValue(besta, belief); 
  return besta;
}
// non-factored beliefs
int POMDP::getPolicy(DdNode *belief, double & beliefValue)
{
  int besta;
  beliefValue = computeValue(besta, belief); 
  return besta;
}

int *POMDP::valueIteration(int horizon, double *beliefValues, char *outpath)
{
  int i,a,o,b;
  int *policy;
  // initial alpha vectors (at most numactions of them)
  //numalphas = numactions;
  if (numalphas == 0) {
    numalphas = 1;
    alphas = new DdNode*[numalphas];
    // this can reduce the number of alpha vectors (if reward+cost for different actions are the same)
    getInitialAlphas(1);
  } else {
    // already have alphas
  }
  if (dosampling) {
    bestactions = fastPBVISampled2(horizon,outpath);
  } else {
    bestactions = fastPBVI(horizon, outpath);
  }
  policy = getPolicy(bestactions, numbeliefs, beliefs, beliefValues);
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
// factored beliefs: beliefs[i][j] is the ith belief distribution over the jth original variable
int * POMDP::fastPBVI(int horizon, char *outpath)
{
  int i,k,a,o,b;
  // the temporary belief set
  DdNode ***tbeliefs;
  char buf[256];
  tbeliefs = new DdNode **[numbeliefs]; 
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
  DdNode *sbelief, *temp;
  DdNode **newalphas;
  bestactions = new int[numalphas];
  int *newbestactions;
  counter = 0;
  // compute value of all beliefs
  for (b=0; b<numbeliefs; b++) {
    Vb[b] = computeValue(bestalpha, beliefs[b]);
    fprintf(stderr," %f",Vb[b]);
    oldVb[b] = -100000;
  }
  fprintf(stderr,"\n");
  
  // also compute the forward propagated belief states
  baz = new DdNode***[numactions];
  for (i=0; i<numactions; i++) {
    baz[i] = new DdNode**[numorigobs];
    for (o=0; o<numorigobs; o++) 
      baz[i][o] = new DdNode *[numbeliefs];
  }
  
  
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      temp = bayesianUpdate(beliefs[b],a);
      for (o=0; o<numorigobs; o++) {
	// multiply by observation function
	baz[a][o][b] = Cudd_addApply(gbm,Cudd_addTimes,totalObsFun[a][o],temp);
	Cudd_Ref(baz[a][o][b]);
      }
      Cudd_RecursiveDeref(gbm,temp);
    }
  }


  // main loop commences
  lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  while (lastiteration < 1){ 
    fprintf(stderr,"\niteration %d/%d numbeliefs %d numalphas %d\n",counter,horizon, numbeliefs, numalphas);
    //printAlphas();
    nb = numbeliefs;

    
    // values of beliefs - and indices of best alphas
    // alphas are unprimed and beliefs are primed
    // reset the temporary belief set
    for (b=0; b<numbeliefs; b++) {
      tbeliefs[b] = beliefs[b];
      tbindex[b] = b;
      tVb[b] = Vb[b];
      oldVb[b] = Vb[b];
      fprintf(stderr," %f",Vb[b]);
    }
    fprintf(stderr,"\n");
    // max numbeliefs new alphas
    newalphas = new DdNode*[numbeliefs];
    newbestactions = new int[numbeliefs];
    numnewalphas = 0;
    while (nb > 0) {
      // sample one belief from tbeliefs
      sb = (int) floor(nb*((double) rand())/((double) RAND_MAX+1.0));
      // for writing out obs image must know what belief we're doing
      currBeliefObsImage = tbindex[sb];
      // back it up - this gets the value in maxgamma
      newbestactions[numnewalphas] = backup(newalphas+numnewalphas, &maxgamma, tbindex[sb], tbeliefs[sb], outpath);
      fprintf(stderr,"sampled %d which is at %d : old %f new %f  new alpha with action %d ",sb,tbindex[sb],Vb[tbindex[sb]],maxgamma,newbestactions[numnewalphas]);
      if (Vb[tbindex[sb]] > maxgamma)
	fprintf(stderr,"******");
      fprintf(stderr,"  %d",nb);
      Vb[tbindex[sb]] = maxgamma;
      //pdd(newalphas[numnewalphas]);
      
      numnewalphas++;
 
      // remove this value from tbeliefs
      nb--;
      tbeliefs[sb] = tbeliefs[nb];
      // takes tbeliefs[nb]'s ref
      tVb[sb] = tVb[nb];
      tbindex[sb] = tbindex[nb];
      
      // check all remaining belief points
      // only if this is not the last iteration, in which case
      // we want to compute a backup for all belief points
      b=0;
      while (DOFASTPBVI && !lastiteration && b < nb) {
	newVb = computeValue(numnewalphas, bestalpha, newalphas, tbeliefs[b]);
	// remove the belief point if its new value is greater
	if (newVb >= tVb[b]) {
	  fprintf(stderr," -> %d (%d) %g ",nb,tbindex[b],newVb);
	  nb--;
	  tbeliefs[b] = tbeliefs[nb];
	  tVb[b] = tVb[nb];
	  Vb[tbindex[b]] = newVb;
	  tbindex[b] = tbindex[nb];
	} else {
	  b++;
	}
      }
      fprintf(stderr," -> %d\n",nb);
    }
    // copy over all new alpha vectors
    for (i=0; i<numalphas; i++) 
      Cudd_RecursiveDeref(gbm,alphas[i]);
    delete [] alphas;
    delete [] bestactions;
    bool *isdup = new bool[numnewalphas];
    // check for duplicate alpha vectors here
    for (i=0; i<numnewalphas; i++) 
      isdup[i] = false;
    int ndup = 0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) 
	for (k=i+1; k<numnewalphas; k++) 
	  if (!isdup[k] && newalphas[i] == newalphas[k]) {
	    isdup[k] = true;
	    ndup++;
	  }
    }
    oldnumalphas = numalphas;
    numalphas = numnewalphas-ndup;

    alphas = new DdNode*[numalphas];
    bestactions = new int[numalphas];
    k=0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) {
	// takes the ref from newalphas
	// possibly want to *approximate* here
	alphas[k] = newalphas[i];
	bestactions[k] = newbestactions[i];
	k++;
      }
    }
    delete [] newalphas;
    delete [] isdup;
    delete [] newbestactions;

    // write alphas to a file
    sprintf(buf,"%s/alphas%d.mdp",outpath,counter);
    printAlphasFile(buf);
    
    counter = counter + 1;

    lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  }
  delete [] tbeliefs; 
  delete [] Vb;
  delete [] tVb;
  delete [] tbindex;
  return bestactions;
}
// this is Vlassis & Spaan 2003
// returns a list of the best actions to take for each alpha vector int the optimal value function
// factored beliefs: beliefs[i][j] is the ith belief distribution over the jth original variable
int * POMDP::fastPBVISampled(int horizon, char *outpath)
{
  int i,k,a,o,b;
  // the temporary belief set
  DdNode ***tbeliefs;
  tbeliefs = new DdNode **[numbeliefs]; 
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
  DdNode *sbelief, *temp, *temp2;
  DdNode **newalphas;
  int *bestactions = new int[numalphas];
  int *newbestactions;
  counter = 0;
  // compute value of all beliefs
  for (b=0; b<numbeliefs; b++) {
    Vb[b] = computeValue(bestalpha, beliefs[b]);
    fprintf(stderr," %f",Vb[b]);
    oldVb[b] = -100000;
  }
  fprintf(stderr,"\n");
  
  // now, draw a set of samples
  //  for each belief+action combination
  
  obsSampleProb = new DdNode **[numactions];
  fprintf(stderr,"computing sample probabilities");
  FILE **sofd = new FILE*[numorigobs];
  char filename[256];
  for (o=0; o<numorigobs; o++) {
    sprintf(filename,"%s/samples%d.txt",outpath,o);
    sofd[o] = fopen(filename,"w");
  }
  for (a=0; a<numactions; a++) {
    obsSampleProb[a] = new DdNode *[totnumsamples];
    // draw a set of samples
    for (i=0; i<totnumsamples; i++) {
      fprintf(stderr,".");
      // draw a sample from each observation function
      obsSampleProb[a][i] = One;
      Cudd_Ref(obsSampleProb[a][i]);
      for (o=0; o<numorigobs; o++) {
	temp = My_addApply(gbm,drawSample,totalObsFun[a][o],One);
	Cudd_Ref(temp);
	// get the samples out for printing
	extractPrintSamples(sofd[o],temp);
	// swap to unprimed!
	temp2 = Cudd_addSwapVariables(gbm,temp,Array1,Array2,numvars);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp);
	temp = temp2;
	temp2 = Cudd_addApply(gbm,computePDF,totalObsFun[a][o],temp);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp);
	temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,obsSampleProb[a][i]);
	Cudd_Ref(temp);
	Cudd_RecursiveDeref(gbm,temp2);
	Cudd_RecursiveDeref(gbm,obsSampleProb[a][i]);
	obsSampleProb[a][i] = temp;
      }
    }
  }
  for (o=0; o<numorigobs; o++) 
    fclose(sofd[o]);
  fprintf(stderr,"computed sample probabilities\n");
  // compute the forward propagated belief states
  baz = new DdNode***[numactions];
  for (i=0; i<numactions; i++) {
    baz[i] = new DdNode**[totnumsamples];
    for (o=0; o<totnumsamples; o++) 
      baz[i][o] = new DdNode *[numbeliefs];
  }
  
  
  fprintf(stderr,"computing forward propagated beliefs",a,b);
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      temp = bayesianUpdate(beliefs[b],a);
      fprintf(stderr,"action %d belief %d\n",a,b);
      for (i=0; i<totnumsamples; i++) {
	fprintf(stderr,".");
	// multiply by sample PDF 
	temp2 = Cudd_addApply(gbm,Cudd_addTimes,obsSampleProb[a][i],temp);
	Cudd_Ref(temp2);
	baz[a][i][b] = NULL;
	normalizeFunction(temp2,&baz[a][i][b],true);
	Cudd_RecursiveDeref(gbm,temp2);
      }
      fprintf(stderr,"\n");
      Cudd_RecursiveDeref(gbm,temp);
    }
  }

  fprintf(stderr,"computing forward propagated beliefs");

  // main loop commences
  lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  while (lastiteration < 2){ 

    //print out info
    //Cudd_PrintInfo(gbm,stderr);

    fprintf(stderr,"\niteration %d/%d numbeliefs %d numalphas %d\n",counter,horizon, numbeliefs, numalphas);
    //printAlphas();
    nb = numbeliefs;
    // values of beliefs - and indices of best alphas
    // alphas are unprimed and beliefs are primed
    // reset the temporary belief set
    for (b=0; b<numbeliefs; b++) {
      tbeliefs[b] = beliefs[b];
      tbindex[b] = b;
      tVb[b] = Vb[b];
      oldVb[b] = Vb[b];
      fprintf(stderr," %f",Vb[b]);
    }
    fprintf(stderr,"\n");
    // max numbeliefs new alphas
    newalphas = new DdNode*[numbeliefs];
    newbestactions = new int[numbeliefs];
    numnewalphas = 0;
    fprintf(stderr,"collecting garbage for the %d time, %d nodes collected\n",Cudd_ReadGarbageCollections(gbm)+1,cuddGarbageCollect(gbm,1));
    while (nb > 0) {
      // sample one belief from tbeliefs
      sb = (int) floor(nb*((double) rand())/((double) RAND_MAX+1.0));
      // for writing out obs image must know what belief we're doing
      currBeliefObsImage = tbindex[sb];
      // back it up - this gets the value in maxgamma
      newbestactions[numnewalphas] = backup(newalphas+numnewalphas, &maxgamma, tbindex[sb], tbeliefs[sb], outpath);
      fprintf(stderr,"sampled %d which is at %d : old %f new %f  new alpha with action %d ",sb,tbindex[sb],Vb[tbindex[sb]],maxgamma,newbestactions[numnewalphas]);
      if (Vb[tbindex[sb]] > maxgamma)
	fprintf(stderr,"******");
      fprintf(stderr,"  %d",nb);
      Vb[tbindex[sb]] = maxgamma;
      //pdd(newalphas[numnewalphas]);
      
      numnewalphas++;
 
      // remove this value from tbeliefs
      nb--;
      tbeliefs[sb] = tbeliefs[nb];
      // takes tbeliefs[nb]'s ref
      tVb[sb] = tVb[nb];
      tbindex[sb] = tbindex[nb];
      
      // check all remaining belief points
      // only if this is not the last iteration, in which case
      // we want to compute a backup for all belief points
      b=0;
      while (!lastiteration && b < nb) {
	newVb = computeValue(numnewalphas, bestalpha, newalphas, tbeliefs[b]);
	// remove the belief point if its new value is greater
	if (newVb >= tVb[b]) {
	  fprintf(stderr," -> %d (%d) %g ",nb,tbindex[b],newVb);
	  nb--;
	  tbeliefs[b] = tbeliefs[nb];
	  tVb[b] = tVb[nb];
	  Vb[tbindex[b]] = newVb;
	  tbindex[b] = tbindex[nb];
	} else {
	  b++;
	}
      }
      fprintf(stderr," -> %d\n",nb);
    }
    // copy over all new alpha vectors
    for (i=0; i<numalphas; i++) 
      Cudd_RecursiveDeref(gbm,alphas[i]);
    delete [] alphas;
    delete [] bestactions;
    bool *isdup = new bool[numnewalphas];
    // check for duplicate alpha vectors here
    for (i=0; i<numnewalphas; i++) 
      isdup[i] = false;
    int ndup = 0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) 
	for (k=i+1; k<numnewalphas; k++) 
	  if (!isdup[k] && newalphas[i] == newalphas[k]) {
	    isdup[k] = true;
	    ndup++;
	  }
    }
    oldnumalphas = numalphas;
    numalphas = numnewalphas-ndup;

    alphas = new DdNode*[numalphas];
    bestactions = new int[numalphas];
    k=0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) {
	// takes the ref from newalphas
	// possibly want to *approximate* here
	alphas[k] = newalphas[i];
	bestactions[k] = newbestactions[i];
	k++;
      }
    }
    delete [] newalphas;
    delete [] isdup;
    delete [] newbestactions;
    
    counter = counter + 1;

    if (lastiteration) {
      lastiteration++;
    } else {
      lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
    }
  }
  delete [] tbeliefs; //free(tbeliefs);
  delete [] Vb;
  delete [] tVb;
  delete [] tbindex;
  return bestactions;
}
// this is Vlassis & Spaan 2003
// returns a list of the best actions to take for each alpha vector int the optimal value function
// factored beliefs: beliefs[i][j] is the ith belief distribution over the jth original variable
// second version using sampling from proposal
int * POMDP::fastPBVISampled2(int horizon, char *outpath)
{
  int i,k,a,o,b;
  // the temporary belief set
  DdNode ***tbeliefs;
  char buf[256];
  tbeliefs = new DdNode **[numbeliefs];
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
  DdNode *sbelief, *temp, *temp2;
  DdNode **newalphas;
  bestactions = new int[numalphas];
  int *newbestactions;
  counter = 0;
  // compute value of all beliefs
  for (b=0; b<numbeliefs; b++) {
    Vb[b] = computeValue(bestalpha, beliefs[b]);
    fprintf(stderr," %f",Vb[b]);
    oldVb[b] = -100000;
  }
  fprintf(stderr,"\n");
  
  // first , compute the b^a (belief states propagated forward without observations)
  // b^a_j = \sum_{s} T(s,a,s')b_j(s)
  ba = new DdNode **[numactions];
  for (a=0; a<numactions; a++) {
    ba[a] = new DdNode *[numbeliefs];
    for (b=0; b<numbeliefs; b++) {
      temp = bayesianUpdate(beliefs[b],a);
      ba[a][b] =  NULL;
      normalizeFunction(temp,&ba[a][b],true);
      Cudd_RecursiveDeref(gbm,temp);
    }
  }
  // now, draw a set of samples
  //  for each belief+action combination
  proposalsamples = new MixGauss****[numactions];
  for (a=0; a<numactions; a++) {
    proposalsamples[a] = new MixGauss***[numbeliefs];
    for (b=0; b<numbeliefs; b++) {
      proposalsamples[a][b] = new MixGauss**[totnumsamples];
      for (i=0; i<totnumsamples; i++) {
	proposalsamples[a][b][i] = new MixGauss*[numorigobs];
      }
    }
  }
  // draw the samples
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) {
	drawSampleFromBelief(ba[a][b],a,proposalsamples[a][b][i]);
      }
    }
  }
  // write samples to a file
  FILE *sofd;
  char filename[256];
  sprintf(filename,"%s/samples.txt",outpath);
  sofd = fopen(filename,"w");
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) {
	fprintf(sofd,"%d %d %d ",a,b,i);
	for (o=0; o<numorigobs; o++) 
	  for (k=0; k<orig_obs[o].type; k++) 
	    fprintf(sofd,"%f ",proposalsamples[a][b][i][o]->mixweights[k]);
	fprintf(sofd,"\n");
      }
    }
  }
  fclose(sofd);

  //compute p_o^{ij} (product of observation functions over obs variables for each sample)

  fprintf(stderr,"computing sample probabilities");
  pojk = new DdNode ***[numactions];
  MixGauss themg;
  for (a=0; a<numactions; a++) {
    pojk[a] = new DdNode **[numbeliefs];
    for (b=0; b<numbeliefs; b++) {
      pojk[a][b] = new DdNode *[totnumsamples];
      for (i=0; i<totnumsamples; i++) {
	pojk[a][b][i] = One;
	Cudd_Ref(pojk[a][b][i]);
	for (o=0; o<numorigobs; o++) {
	  temp = Cudd_addConst(gbm,proposalsamples[a][b][i][o]);
	  Cudd_Ref(temp);
	  temp2 = Cudd_addApply(gbm,computePDF,totalObsFun[a][o],temp);
	  Cudd_Ref(temp2);
	  Cudd_RecursiveDeref(gbm,temp);
	  temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,pojk[a][b][i]);
	  Cudd_Ref(temp);
	  Cudd_RecursiveDeref(gbm,temp2);
	  Cudd_RecursiveDeref(gbm,pojk[a][b][i]);
	  pojk[a][b][i] = temp;

	  if (pojk[a][b][i] == Zero) {
	    fprintf(stderr,"\npojk  is Zero at %d %d %d %d: %s\n",a,b,i,o,proposalsamples[a][b][i][o]->toString());
	    fprintf(stderr,"\n setting to uniform ..\n");
	    
	  }
	}
      }
    }
  }
  fprintf(stderr,"computed sample probabilities\n");

  // compute the forward propagated belief states
  baz = new DdNode***[numactions];
  for (i=0; i<numactions; i++) {
    baz[i] = new DdNode**[totnumsamples];
    for (o=0; o<totnumsamples; o++) {
      baz[i][o] = new DdNode*[numbeliefs];
    }
  }
  
  fprintf(stderr,"computing forward propagated beliefs",a,b);
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) {
	// multiply by sample PDF 
	temp2 = Cudd_addApply(gbm,Cudd_addTimes,pojk[a][b][i],ba[a][b]);
	Cudd_Ref(temp2);
	baz[a][i][b] = temp2;
	//baz[a][i][b] = NULL;
	//normalizeFunction(temp2,&baz[a][i][b],true);
	//Cudd_RecursiveDeref(gbm,temp2);
      }
    }
  }

  //compute najk - baz summed over s'
  // also divide out the number of samples
  MixGauss mg;
  mg.set(totnumsamples);
  DdNode *numsampdd = Cudd_addConst(gbm,&mg);
  Cudd_Ref(numsampdd);

  fajk = new DdNode ***[numactions];
  for (a=0; a<numactions; a++) {
    fajk[a] = new DdNode **[totnumsamples];
    for (i=0; i<totnumsamples; i++) 
      fajk[a][i] = new DdNode *[numbeliefs];
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) {
	temp = sumOutAllPrimedOrigVars(baz[a][i][b]);
	sval = (*Cudd_V(temp)).get_val();
	fajk[a][i][b] = Cudd_addApply(gbm,Cudd_addDivide,pojk[a][b][i],temp);
	Cudd_Ref(fajk[a][i][b]);
	Cudd_RecursiveDeref(gbm,temp);
      }
    }
  }
  Cudd_RecursiveDeref(gbm,numsampdd);
  
  /*
  // can delete ba, proposalsamples, and pojk now
  // delete proposal samples
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) {
	for (o=0; o<numorigobs; o++) 
	  delete proposalsamples[a][b][i][o];
	delete [] proposalsamples[a][b][i];
      }
      delete [] proposalsamples[a][b];
    }
    delete []  proposalsamples[a];
  }
  delete [] proposalsamples;
  */
  // delete ba
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) 
      Cudd_RecursiveDeref(gbm,ba[a][b]);
    delete [] ba[a];
  } 
  delete [] ba;
  //delete pojk
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      for (i=0; i<totnumsamples; i++) 
	Cudd_RecursiveDeref(gbm,pojk[a][b][i]);
      delete [] pojk[a][b];
    }
    delete [] pojk[a];
  }
  delete [] pojk;
  // garbage collect to try to clean up a bit
  fprintf(stderr,"collecting garbage, %d nodes collected\n",cuddGarbageCollect(gbm,1));
  
  
  // main loop commences
  lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  while (lastiteration < 1) {
    //print out info
    //Cudd_PrintInfo(gbm,stderr);

    fprintf(stderr,"\niteration %d/%d numbeliefs %d numalphas %d\n",counter,horizon, numbeliefs, numalphas);
    //printAlphas();
    nb = numbeliefs;
    // values of beliefs - and indices of best alphas
    // alphas are unprimed and beliefs are primed
    // reset the temporary belief set
    for (b=0; b<numbeliefs; b++) {
      tbeliefs[b] = beliefs[b];
      tbindex[b] = b;
      tVb[b] = Vb[b];
      oldVb[b] = Vb[b];
      fprintf(stderr," %f",Vb[b]);
    }
    fprintf(stderr,"\n");
    // max numbeliefs new alphas
    newalphas = new DdNode*[numbeliefs];
    newbestactions = new int[numbeliefs];
    numnewalphas = 0;
    //fprintf(stderr,"collecting garbage for the %d time, %d nodes collected\n",Cudd_ReadGarbageCollections(gbm)+1,cuddGarbageCollect(gbm,1));
    while (nb > 0) {
      // sample one belief from tbeliefs
      sb = (int) floor(nb*((double) rand())/((double) RAND_MAX+1.0));
      // for writing out obs image must know what belief we're doing
      currBeliefObsImage = tbindex[sb];
      // back it up - this gets the value in maxgamma
      newbestactions[numnewalphas] = backup(newalphas+numnewalphas, &maxgamma, tbindex[sb], tbeliefs[sb], outpath);
      fprintf(stderr,"sampled %d which is at %d : old %f new %f  new alpha with action %d ",sb,tbindex[sb],Vb[tbindex[sb]],maxgamma,newbestactions[numnewalphas]);
      if (Vb[tbindex[sb]] > maxgamma)
	fprintf(stderr,"******");
      fprintf(stderr,"  %d",nb);
      Vb[tbindex[sb]] = maxgamma;
      //pdd(newalphas[numnewalphas]);
      
      numnewalphas++;
 
      // remove this value from tbeliefs
      nb--;
      tbeliefs[sb] = tbeliefs[nb];
      // takes tbeliefs[nb]'s ref
      tVb[sb] = tVb[nb];
      tbindex[sb] = tbindex[nb];
      
      // check all remaining belief points
      // only if this is not the last iteration, in which case
      // we want to compute a backup for all belief points
      b=0;
      while (b< nb) {
	newVb = computeValue(numnewalphas, bestalpha, newalphas, tbeliefs[b]);
	// remove the belief point if its new value is greater
	if (newVb >= tVb[b]) {
	  fprintf(stderr," -> %d (%d) %g ",nb,tbindex[b],newVb);
	  nb--;
	  tbeliefs[b] = tbeliefs[nb];
	  tVb[b] = tVb[nb];
	  Vb[tbindex[b]] = newVb;
	  tbindex[b] = tbindex[nb];
	} else {
	  b++;
	}
      }
      fprintf(stderr," -> %d\n",nb);
    }
    // copy over all new alpha vectors
    for (i=0; i<numalphas; i++) 
      Cudd_RecursiveDeref(gbm,alphas[i]);
    delete [] alphas;
    delete [] bestactions;
    bool *isdup = new bool[numnewalphas];
    // check for duplicate alpha vectors here
    for (i=0; i<numnewalphas; i++) 
      isdup[i] = false;
    int ndup = 0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) 
	for (k=i+1; k<numnewalphas; k++) 
	  if (!isdup[k] && newalphas[i] == newalphas[k]) {
	    isdup[k] = true;
	    ndup++;
	  }
    }
    oldnumalphas = numalphas;
    numalphas = numnewalphas-ndup;

    alphas = new DdNode*[numalphas];
    bestactions = new int[numalphas];
    k=0;
    for (i=0; i<numnewalphas; i++) {
      if (!isdup[i]) {
	// takes the ref from newalphas
	// possibly want to *approximate* here
	alphas[k] = newalphas[i];
	bestactions[k] = newbestactions[i];
	k++;
      }
    }
    delete [] newalphas;
    delete [] isdup;
    delete [] newbestactions;
    // write alphas to a file
    sprintf(buf,"%s/alphas%d.mdp",outpath,counter);
    printAlphasFile(buf);
    
    counter = counter + 1;
    lastiteration = converged(Vb,oldVb,numbeliefs,counter,horizon,tolerance);
  }
  delete [] tbeliefs; //free(tbeliefs);
  delete [] Vb;
  delete [] tVb;
  delete [] tbindex;
  return bestactions;
}
// computes the value of a factored belief vector (max over alpha vectors)
// returns the value and maxi gets set to the index of the alpha vector that
// was maximal
double POMDP::computeValue(int & maxi, DdNode **belief)
{
  return computeValue(numalphas,maxi,alphas,belief);
}

double POMDP::computeValue(int & maxi, DdNode *belief)
{
  return computeValue(numalphas,maxi,alphas,belief);
}
// computes the value of a factored belief vector (max over alpha vectors)
// with alpha vectors passed in
// returns the value and maxi gets set to the index of the alpha vector that
// was maximal
double POMDP::computeValue(int numalphs, int & maxi, DdNode **alphs, DdNode **belief)
{
  int a,j,i;
  double tmg, maxgamma = -10000000.0;
  maxi = -1;
  DdNode *temp, *temp1, *alphap;
  for (a=0; a<numalphs; a++) {
    alphap = computeValue(alphs[a],belief);
    if ((tmg = (*Cudd_V(alphap)).get_min()) > maxgamma) {
      maxgamma = tmg;
      maxi = a;
    }
    Cudd_RecursiveDeref(gbm,alphap);
  }
  return maxgamma;
}
double POMDP::computeValue(int numalphs, int & maxi, DdNode **alphs, DdNode *belief)
{
  int a,j,i;
  double tmg, maxgamma = -10000000.0;
  maxi = -1;
  DdNode *temp, *temp1, *alphap;
  for (a=0; a<numalphs; a++) {
    alphap = computeValue(alphs[a],belief);
    if ((tmg = (*Cudd_V(alphap)).get_min()) > maxgamma) {
      maxgamma = tmg;
      maxi = a;
    }
    Cudd_RecursiveDeref(gbm,alphap);
  }
  return maxgamma;
}
// compute value for a single alpha vector at factored belief belief
DdNode *POMDP::computeValue(DdNode *alpha, DdNode **belief) 
{
  int j;
  DdNode *alphap, *temp;
  //     dp = alpha . b
  // switch alpha to primed
  alphap = Cudd_addSwapVariables(gbm,alpha,Array1,Array2,numvars);
  Cudd_Ref(alphap);
  
  
  temp = multiplySumSet(alphap,belief);
  Cudd_RecursiveDeref(gbm,alphap);
  alphap = temp;
  return alphap;
}
// compute value for a single alpha vector at unfactored belief belief
DdNode *POMDP::computeValue(DdNode *alpha, DdNode *belief) 
{
  int j;
  DdNode *alphap, *temp, *temp2;
  //     dp = alpha . bn
  // switch alpha to primed
  alphap = Cudd_addSwapVariables(gbm,alpha,Array1,Array2,numvars);
  Cudd_Ref(alphap);
  temp = Cudd_addApply(gbm,Cudd_addTimes,alphap,belief);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,alphap);
  alphap = temp;
  // this actually sums over all variables, so temp should be a constant at the end
  temp = sumOutAllPrimedOrigVars(alphap);
  Cudd_RecursiveDeref(gbm,alphap);
  alphap = temp;

  return alphap;
}

// given the set of backed up alpha vectors for each action and observation, gamma_ao,
// generates a single new alpha vector by for belief pbelief
// this alpha is copied out in *palpha, and the best action(s) to take
// at this belief point are retured in binary (a 1 at pos i if action i is one of the best)
// maxgamma is the resulting value at the belief point
// this version with factored beliefs
int POMDP::backup(DdNode **palpha, double * maxgamma, int pbelief_ind, DdNode **pbelief, char *outpath)
{
  int i,j,a,o,ov,b,tmp,oo;
  MixGauss * alphaTol = new MixGauss(0.001);
  int dupalpha;
  DdNode *temp1,*temp, *temp2;
  
  //*************  Construct gamma_ab (step 2)
  DdNode **gamma_ab = new DdNode*[numactions];

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

    if (dosampling) {
      // do the backup by computing the sample-based aggregate observation function
      temp2 = sampleBackup2(pbelief_ind,a, outpath); 
    } else {
      // do the backup by partitioning the observation space, summing over regions
      // and doing a Bellman backup
      temp2 = partitionIntegrateObsSpace(pbelief_ind,a);
    }
    temp = Cudd_addApply(gbm,Cudd_addPlus,temp2,gamma_ab[a]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,gamma_ab[a]);
    gamma_ab[a] = temp;
  }
  //*************  select best action for this belief point (step 3)
  // new alpha list  = empty
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
  delete [] gamma_ab; 
  return maxi;
}
DdNode *POMDP::sampleBackup(int pbelief_ind, int a)
{
  int i,j, dpind,k;
  DdNode *temp, *temp2, *intval, *maxdd;
  double dpval,sum,maxdpval;

  //DdNode **alphap;
  //alphap = new DdNode *[numalphas];
  DdNode **maxindexsum = new DdNode*[numalphas];
  DdNode **nmaxindex = new DdNode*[numalphas];
  for (k=0; k<numalphas; k++) {
    maxindexsum[k] = Zero;
    Cudd_Ref(maxindexsum[k]);
  }
  for (i=0; i<totnumsamples; i++) {
    for (j=0; j<numalphas; j++) {
      temp = computeValue(alphas[j],baz[a][i][pbelief_ind]);
      if (j==0) {
	// maxdd takes temp's ref
	maxdd = temp;
	nmaxindex[j] = One;
	Cudd_Ref(nmaxindex[j]);
      } else {
        temp2 = Cudd_addApply(gbm,Cudd_addMaximum,maxdd,temp);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,maxdd);
	maxdd = temp2;
	temp2 = Cudd_addApply(gbm,Cudd_addMinus,temp,maxdd);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp);
	temp = Cudd_addBddThreshold(gbm,temp2,pZero);
	Cudd_Ref(temp);
	Cudd_RecursiveDeref(gbm,temp2);
	nmaxindex[j] = Cudd_BddToAdd(gbm,temp);
	Cudd_Ref(nmaxindex[j]);
	Cudd_RecursiveDeref(gbm,temp);
	temp2 = Cudd_addApply(gbm,Cudd_addMinus,One,nmaxindex[j]);
	Cudd_Ref(temp2);
	for (k=0; k<j; k++) {
	  temp = Cudd_addApply(gbm,Cudd_addTimes,nmaxindex[k],temp2);
	  Cudd_Ref(temp);
	  Cudd_RecursiveDeref(gbm,nmaxindex[k]);
	  nmaxindex[k] = temp;
	}
	Cudd_RecursiveDeref(gbm,temp2);
      }
    }
    Cudd_RecursiveDeref(gbm,maxdd);
    for (k=0; k<numalphas; k++) {
      temp = Cudd_addApply(gbm,Cudd_addPlus,maxindexsum[k],nmaxindex[k]);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,nmaxindex[k]);
      Cudd_RecursiveDeref(gbm,maxindexsum[k]);
      maxindexsum[k] = temp;
    }
  }
  intval = Zero;
  Cudd_Ref(intval);
  MixGauss mg;
  mg.set(totnumsamples);
  DdNode *numsampdd = Cudd_addConst(gbm,&mg);
  Cudd_Ref(numsampdd);
  
  for (i=0; i<numalphas; i++) {
    // compute observation probability
    temp2 = Cudd_addApply(gbm,Cudd_addDivide,maxindexsum[i],numsampdd);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,maxindexsum[i]);
    // multiply 
    temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,alphas[i]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);

    // swap to primed variables
    temp2 = Cudd_addSwapVariables(gbm,temp,Array2,Array1,numvars);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    temp = temp2;
    
    // then by the transition function and sum out primed variable
    temp2 = multiplySumSet(temp,NewPrime[a]);
    Cudd_RecursiveDeref(gbm,temp);

    // finally, multiply by the discount
    temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,discount);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);

    // add to sum
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp,intval);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,intval);
    intval = temp2;
  }
  Cudd_RecursiveDeref(gbm,numsampdd);
  delete [] maxindexsum;
  delete [] nmaxindex;
  return intval;
}
DdNode *POMDP::sampleBackup2(int pbelief_ind, int a, char *outpath)
{
  int i,j, o, dpind,k, index;
  DdNode *temp, *temp2, *intval, *maxdd;
  double dpval,sum,maxdpval;

  DdNode **maxindexsum = new DdNode*[numalphas];
  DdNode **nmaxindex = new DdNode*[numalphas];
  for (k=0; k<numalphas; k++) {
    maxindexsum[k] = Zero;
    Cudd_Ref(maxindexsum[k]);
  }
  int samplesinimage=0;
  for (i=0; i<totnumsamples; i++) {
    // find the best alpha vector
    for (j=0; j<numalphas; j++) {
      temp = computeValue(alphas[j],baz[a][i][pbelief_ind]);
      if (!Cudd_IsConstant(temp)) {
	fprintf(stderr,"Where is Yayo? \n");
	exit(0);
      }
      dpval = (*Cudd_V(temp)).get_val();
      if (j==0 || dpval > maxdpval) {
	maxdpval = dpval;
	dpind = j;
      }
      Cudd_RecursiveDeref(gbm,temp);
    }
    // dpind is best alpha and maxdpval is its value for sample proposalsample[a][pblief_ind][i]
    temp = Cudd_addApply(gbm,Cudd_addPlus,maxindexsum[dpind],fajk[a][i][pbelief_ind]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,maxindexsum[dpind]);
    maxindexsum[dpind] = temp;
  }
  // renormalize maxindexsum
  temp = Zero;
  Cudd_Ref(temp);
  for (j=0; j<numalphas; j++) {
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,maxindexsum[j],temp);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    temp = temp2;
  }
  for (j=0; j<numalphas; j++) {
    temp2 = Cudd_addApply(gbm,Cudd_addDivide,maxindexsum[j],temp);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,maxindexsum[j]);
    maxindexsum[j] = temp2;
  }
  intval = Zero;
  Cudd_Ref(intval);

  // back up the alpha vectors
  for (i=0; i<numalphas; i++) {
    // multiply  by alpha
    temp = Cudd_addApply(gbm,Cudd_addTimes,maxindexsum[i],alphas[i]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,maxindexsum[i]);

    // swap to primed variables
    temp2 = Cudd_addSwapVariables(gbm,temp,Array2,Array1,numvars);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    temp = temp2;
    
    // then by the transition function and sum out primed variable
    temp2 = multiplySumSet(temp,NewPrime[a]);
    Cudd_RecursiveDeref(gbm,temp);

    // finally, multiply by the discount
    temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,discount);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);

    // add to sum
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp,intval);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,intval);
    intval = temp2;
  }
  delete [] maxindexsum;
  delete [] nmaxindex;
  return intval;
}

// factored belief state versions
// expand beliefs over a horizon
// add beliefs only if  they differ (sup norm) by more than mthresh to some other belief
// if allbeliefs = false, only add the best belief for each action at each step
void POMDP::expandBeliefs(double horizon, double mthresh, bool allbeliefs)
{
  int i,j;
  for (i=0; i<horizon; i++) {
    fprintf(stderr,"\niteration %d numbeliefs %d number of gcs so far %d\n",i,numbeliefs,Cudd_ReadGarbageCollections(gbm));
    expandBeliefs(mthresh,allbeliefs);
  }
}

void POMDP::expandBeliefs(int minnumbeliefs, double mthresh, bool allbeliefs)
{
  int i,j;
  i=0;
  int onb;
  while (numbeliefs < minnumbeliefs) {
    onb = numbeliefs;
    fprintf(stderr,"\niteration %d numbeliefs %d number of gcs so far %d\n",i,numbeliefs,Cudd_ReadGarbageCollections(gbm));
    expandBeliefs(mthresh,allbeliefs);
    // found all beliefs at this mthresh, so divide mthresh
    if (numbeliefs == onb) {
      mthresh = mthresh/2.0;
      fprintf(stderr,"reducing mthresh to %f\n",mthresh);
    }
  }
}

// expands the set of beliefs
// input pbeliefs are over primed variables
// for each action, sample a next belief state, then an observation
// from that, then update the belief based on that observation and action,
// but only take it if it is different from a belief already in the 
// set (sup. norm close by mthresh), if allbeliefs = false,
// only take the furthest belief over all the actions, otherwise, take
// all of them
void POMDP::expandBeliefs(double mthresh, bool allbeliefs)
{
  int i,j,a,b,c,o;
  double dist,*dists,maxdist;
  DdNode *temp;
  DdNode **maxb, **ftemp, **bestb;
  //  new beliefs list <- empty
  // its at most twice the size of the old list
  
  int numnewbeliefs;
  if (allbeliefs)
    numnewbeliefs = numbeliefs*numactions;
  else
    numnewbeliefs = 2*numbeliefs;
  
  DdNode ***newbeliefs = new DdNode**[numnewbeliefs];
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
  double r;
  for (b=0; b<numbeliefs; b++) {
    maxdist = 0.0;
    bestb = NULL;
    for (a=0; a<numactions; a++) {
      //   b' <-- simulate one step forward from b on action a
      // newbeliefs is over primed variables
      // have to check if its zero
      do {
	temp = simulateAction(beliefs[b],a);
      } while (temp == Zero);
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
  delete [] beliefs;
  beliefs = newbeliefs;
  numbeliefs = numnewbeliefs;
}
// computes the product of a factored belief state
DdNode *POMDP::beliefProduct(DdNode **b) 
{
  DdNode *temp;
  DdNode *belf = One;
  Cudd_Ref(belf);
  for (int i=0; i<numorigvars; i++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,belf,b[i]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,belf);
    belf = temp;
  }
  return belf;
}

// computes the marginal distributions over b
// returns pointer to an array of numorigvars DdNode *
// one for each marginal belief
DdNode **POMDP::beliefMarginal(DdNode *b) 
{
  int i,j;
  DdNode *temp, *temp1;
  double tmpval;
  DdNode **bm = new DdNode*[numorigvars]; 
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
// simluator for  POMDP
// starts at initb belief 
// from which the actual initial state, astate,  is drawn
// initb -> belf (the belief state)
// for stages iterations
//   consult the policy at belf -> a
//   take action a, update astate -> astatep (returns reward and observation)
//                  update belf -> belfp using observation
//   belfp -> belf;
//   astatep -> astate
// end
// assumes  alpha vectors and bestactions have already been computed
double POMDP::simulator(int stages, char *outpath, bool verbose, bool funky)
{
  return simulator(initBeliefState, stages, outpath, verbose, funky);
}
double POMDP::simulator(DdNode **initb, int stages, char *outpath, bool verbose, bool funky)
{
  int o,t,i,j,a, ap, tmp, acttotake;
  double bv;
  MixGauss rew;
  char buf[128];
  DdNode *astate, *belf, *astatep, *belfp, **belff, *obsprob, *temp;
  DdNode **rewcost;
  astate = sample_Belief(initb);

  // swap astate to unprimed
  temp = Cudd_addSwapVariables(gbm,astate,Array2,Array1,numvars);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,astate);
  astate = temp;

  belf = beliefProduct(initb);
  double cumrew(0.0);
  rewcost = new DdNode*[numactions];
  for (a=0; a<numactions; a++) {
    rewcost[a] = Cudd_addApply(gbm,Cudd_addPlus,RewardD,actionCost[a]);
    Cudd_Ref(rewcost[a]);
  }
  double df =  (*(Cudd_V(discount))).get_max();
  double currdf = df;

  for (t=0; t<stages; t++) {
    // get best action to take
    tmp = getPolicy(bestactions, belf, bv);
    
    // choose one action
    if (verbose) {
      fprintf(stderr,"belief is \n");
      printBelief(belf,stderr);
      fprintf(stderr,"actions to take are "); 
    }
    ap = 0;
    while (tmp >= 1) {
      if (tmp%2) {
	if (verbose) 
	  fprintf(stderr," %s",actionlist[ap].name);
	a = ap;
      }
      ap++;
      tmp = tmp/2;
    }
    acttotake= a;
    if (verbose)
      fprintf(stderr,"\n ******************************* taking action %s\n",actionlist[a].name);

    // get reward at that state
    rew.copyFrom(evaluateDdAtCube(rewcost[a],astate,false));
    cumrew += currdf*rew.val;
    currdf *= df;
    // print this out
    if (verbose) 
      fprintf(stderr,"reward was %f    total  after %d stages is %f\n",rew.val,t,cumrew);
    
    // update the state with that action
    astatep = simulateTransition(astate, a);


    if (verbose) {
      fprintf(stderr,"\nactual state is \n");
      printBelief(astatep,stderr);
    }
    
    // sample an observation from astatep and compute prob. of observation based on it
    if (funky) {
      obsprob = computeSampleObsProdFunky(astatep,a,verbose);
    } else {
      obsprob = computeSampleObsProd(astatep,a,verbose);
    }
    // swap astatep to unprimed and replace astate
    temp = Cudd_addSwapVariables(gbm,astatep,Array2,Array1,numvars);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,astatep);
    Cudd_RecursiveDeref(gbm,astate);
    astate = temp;

    if (verbose) {
      //fprintf(stdout,"print conditional plan (y/n)? ");
      //fscanf(stdin,"%c",&buf);
      //buf[0] = 'y';
      //if (buf[0] == 'y'){
      if (1) {
	sprintf(buf,"%s/sim%d",outpath,t);
	for (o=0; o<numorigobs; o++) {
	  if (orig_obs[o].type == 2) {
	    printConditionalPlanRegions3(buf,acttotake,belf,o);
	  }
	}
      }
    }
    // update belf with this observation
    belfp = bayesianUpdate(belf,a,obsprob);
    Cudd_RecursiveDeref(gbm,belf);
    Cudd_RecursiveDeref(gbm,obsprob);
    belf = belfp;
  }
  for (a=0; a<numactions; a++) 
    Cudd_RecursiveDeref(gbm,rewcost[a]);
  delete [] rewcost;
  Cudd_RecursiveDeref(gbm,belf);
  Cudd_RecursiveDeref(gbm,astate);
  return cumrew;
}
// evaluates thedd at thecube  - returns a constant
MixGauss POMDP::evaluateDdAtCube(DdNode *thedd, DdNode *thecube, bool primed)
{
  DdNode *temp, *temp2;
  temp = Cudd_addApply(gbm,Cudd_addTimes,thedd,thecube);
  Cudd_Ref(temp);
  if (primed) {
    temp2 = sumOutAllPrimedOrigVars(temp);
  } else {
    temp2 = sumOutAllOrigVars(temp);
  }
  Cudd_RecursiveDeref(gbm,temp);
  MixGauss res = (*Cudd_V(temp2));
  Cudd_RecursiveDeref(gbm,temp2);
  return res;
}
// samples a transition from state stsamp on action a
// returns the cube representing the new state
DdNode *POMDP::simulateTransition(DdNode *stsamp, int a)
{
  int j;
  DdNode *newb, *temp1, *temp2;
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
  temp1 = sample_Belief(newb);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,newb);
  return temp1;
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

  temp1 = simulateTransition(stsamp,a);
  Cudd_RecursiveDeref(gbm,stsamp);
  stsamp = temp1;

  obsprod = computeSampleObsProd(stsamp,a);
  // now, we compute the bayesian update to the belief function
  // given this sampled observation
  //multiply belief by transition function and sum out variables
  newb = bayesianUpdate(bb,a,obsprod);
  Cudd_RecursiveDeref(gbm,obsprod);
  Cudd_RecursiveDeref(gbm,stsamp);
  return newb;
}
 
// draws a sample observation from the state stsamp
// and returns the product of observation probabilities of that sample for action a
DdNode *POMDP::computeSampleObsProd(DdNode *stsamp, int a, bool verbose)
{ 
  DdNode *obsprod, *sampo, *temp1, *temp2;
  int o,j;
  // sample an observation from ObsFun[a] at the state given by stsamp
  // sample from each observation separately
  obsprod = One;
  Cudd_Ref(obsprod);
  for (o=0; o<numorigobs; o++) {
    sampo = My_addApply(gbm,drawSample,totalObsFun[a][o],stsamp);
    Cudd_Ref(sampo);

    // but sampo still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampo,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampo);
      sampo = temp1;
    }
    if (verbose) 
      pdd(sampo);
    // now, sampo is a function over the observation variable
    // now we want the ObsFun at that value of observation: sampo
    if (orig_obs[o].type == 0) {
      temp2 = Cudd_addApply(gbm,computePDFType0,totalObsFun[a][o],sampo);
    } else {
      temp2 = Cudd_addApply(gbm,computePDF,totalObsFun[a][o],sampo);
    }
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,sampo);

    // keep a running product
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,obsprod,temp2);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,obsprod);
    obsprod = temp1;
  }
  return obsprod;
}
// funky version where the first obs function is the discrete version we want to 
// evaluate, but the second and third are the continuous ones we want to draw a sample from
// for tiger problem, its simpler
DdNode *POMDP::computeSampleObsProdFunky(DdNode *stsamp, int a, bool verbose)
{ 
  DdNode *obsprod, **sampof, *temp1, *temp2;
  int o,j, index;
  // sample an observation from ObsFun[a] at the state given by stsamp
  // sample from each observation separately
  sampof = new DdNode*[2];
  DdNode *sampo;
  obsprod = One;
  Cudd_Ref(obsprod);
  MixGauss tmp;
  bool noobs = false;
  double maxval;
  int *vv = new int[numorigvars];
  int *vindex = new int[numorigobs];
  int *vardep = new int[numorigobs];

#ifdef FUNKYTIGER
  // dependence of observation function on variables
  vardep[0] = 0;
  for (o=0; o<numorigobs-1; o++) {
    sampo = My_addApply(gbm,drawSample,totalObsFun[a][o+1],stsamp);
    Cudd_Ref(sampo);
    // but sampo still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampo,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampo);
      sampo = temp1;
    }
    if (verbose) 
      pdd(sampo);
    
    temp2 = Cudd_addApply(gbm,computePDF,totalObsFun[a][o+1],sampo);
    Cudd_Ref(temp2);
    sampo = temp2;

    for (j=0; j<numorigvars; j++) 
      vv[j] = 0;
    for (j=0; j<orig_vars[vardep[o]].nvals; j++) {
      vv[vardep[o]] = j;
      temp2 = buildCubeOrig(vv,true);
      tmp = evaluateDdAtCube(sampo,temp2,true);
      Cudd_RecursiveDeref(gbm,temp2);
      if (j==0 || tmp.val > maxval) {
	vindex[o] = j;
	maxval = tmp.val;
      }
    }
  }
  // finally, build index from vindexes
  // need to know how many observations there are here in general
  index = vindex[0];
  if (verbose) 
    fprintf(stderr,"index is %d\n",index);
#else
  for (o=0; o<numorigobs-1; o++) {
    sampof[o] = My_addApply(gbm,drawSample,totalObsFun[a][o+1],stsamp);
    Cudd_Ref(sampof[o]);
    // but sampof still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampof[o],orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampof[o]);
      sampof[o] = temp1;
    }

    if (verbose) 
      pdd(sampof[o]);
    

    // check to see if its  the no-obs observation
    if (o==1) {
      tmp.copyFrom((*Cudd_V(sampof[o])));
      if (tmp.mixweights[0] < 0 && tmp.mixweights[1] < 0) {
	noobs = true;
      }
    }
    
    // now, sampof[o] is an observation o
    // we want to classify it based on totalObsFun to get what discrete observation it corresponds to
    // so get its probability given each state 
    temp2 = Cudd_addApply(gbm,computePDF,totalObsFun[a][o+1],sampof[o]);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,sampof[o]);
    sampof[o] = temp2;
  }
  // get values of variables water flow (observation 1 - variable 2 (2 observations))
  int maxwindex, maxhpindex;
  for (j=0; j<numorigvars; j++)
    vv[j] = 0;

  for (j=0; j<orig_vars[2].nvals; j++) {
    vv[2] = j;
    temp2 = buildCubeOrig(vv,true);
    tmp = evaluateDdAtCube(sampof[0],temp2,true);
    Cudd_RecursiveDeref(gbm,temp2);
    if (j==0 || tmp.val > maxval) {
      maxwindex = j;
      maxval = tmp.val;
    }
  }
  // for  hand position (observation 2 - variable 3 (6 observations))
  if (noobs) {
    maxhpindex = 0;
  } else {
    vv[2] = 0;
    for (j=0; j<orig_vars[3].nvals; j++) {
      vv[3] = j;
      temp2 = buildCubeOrig(vv,true);
      tmp = evaluateDdAtCube(sampof[1],temp2,true);
      Cudd_RecursiveDeref(gbm,temp2);
      if (j==0 || tmp.val > maxval) {
	maxhpindex = 5-j;
	maxval = tmp.val;
      }
    }
  }
  // now reconstruct the discrete observation from maxhpindex and maxwindex
  index = maxwindex*6+maxhpindex;
  if (verbose) 
    fprintf(stderr,"water index %d hp index %d comb index %d\n",maxwindex,maxhpindex,index);
#endif
  // make a new MixGauss of of it
  tmp.set(index);
  temp1 = Cudd_addConst(gbm,&tmp);
  Cudd_Ref(temp1);
  // now we want the ObsFun at that value of observation: sampof
  temp2 = Cudd_addApply(gbm,computePDFType0,totalObsFun[a][0],temp1);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);
  return temp2;
}
// draws a  sample from newb
void POMDP::drawSampleFromBelief(DdNode *newb, int a, MixGauss **v) 
{
  int o,j,k;
  DdNode *stsamp, *sampo, *temp1;
  MixGauss *mg;
  // now newb is new belief state, b, now
  // sample a state from b - returns a cube with the state
  stsamp = sample_Belief(newb);

  // sample an observation from ObsFun[a] at the state given by stsamp
  // sample from each observation separately
  for (o=0; o<numorigobs; o++) {
    sampo = My_addApply(gbm,drawSample,totalObsFun[a][o],stsamp);
    Cudd_Ref(sampo);

    // but sampo still also includes the primed variables in stsamp, so
    // we need to sum these out
    for (j=0; j<numorigvars; j++) {
      temp1 = sumOutPrime(sampo,orig_vars+j,prime_vars);
      Cudd_RecursiveDeref(gbm,sampo);
      sampo = temp1;
    }
    // get value of sampo
    if (!Cudd_IsConstant(sampo)) {
      fprintf(stderr,"I got a funny feeling about this ...\n");
      exit(0);
    }
    v[o] = new MixGauss(*Cudd_V(sampo));
    //for (k=0; k<orig_obs[o].type; k++) 
    //v[o][k] = mg->mixweights[k];
    Cudd_RecursiveDeref(gbm,sampo);
  }
  Cudd_RecursiveDeref(gbm,stsamp);
}
double POMDP::checkNormalization(DdNode *b)
{
  int i;
  MixGauss *theMixGauss;
  int *thecube;
  DdGen * theGen;
  theGen = Cudd_FirstCube(gbm,b,&thecube,&theMixGauss);
  int gen=1;
  double sum(0.0);
  double factor = 1;
  while (gen) {
    // see how many don't care's there are
    factor = 1;
    for (i=0; i<numvars; i++)
      if (thecube[2*i] == 2) 
	factor *= 2;
    sum += factor*(theMixGauss->get_max());
    gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
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
// bayesian update with unfactored beliefs
DdNode * POMDP::bayesianUpdate(DdNode *belief, int a)
{
  int j;
  DdNode *temp, *temp1;
  temp = Cudd_addSwapVariables(gbm,belief,Array2,Array1,numvars);
  Cudd_Ref(temp);

  for (j=0; j<numorigvars; j++) {
    temp1 = Cudd_addApply(gbm,Cudd_addTimes,temp,NewPrime[a][j]);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp);
    temp = temp1;
  }
  //****************  belief could blow up here**************

  temp1 = sumOutAllOrigVars(temp);
  Cudd_RecursiveDeref(gbm,temp);
  temp = temp1;

  return temp;
}
// updates the belief on action a by multiplying by the transition function and summing
DdNode * POMDP::bayesianUpdate(DdNode **b, int a)
{
  int j;
  // first we have to compute the whole belief state
  DdNode *temp;
  DdNode *belief;
  belief = beliefProduct(b);
  temp = bayesianUpdate(belief,a);
  return temp;
}
// multiplies the blief by the obsprod
DdNode *POMDP::multBeliefObsProd(DdNode *b, DdNode *obsprod)
{
  DdNode *temp;
  temp = Cudd_addApply(gbm,Cudd_addTimes,obsprod,b);
  Cudd_Ref(temp);
  DdNode *belief = NULL;
  // if the observation is so unlikely that its probability is 0, we 
  // assume that it didn't even happen so the belief stays the same
  if (temp != Zero){
    normalizeFunction(temp,&belief,true);
  } else {
    belief = b;
    Cudd_Ref(belief);
  }
  Cudd_RecursiveDeref(gbm,temp);
  return belief;
}
//Bayesian update of belief b given action a -
// obsprod gives the observation function at the value of the
// actual observation
DdNode * POMDP::bayesianUpdate(DdNode **b, int a, DdNode *obsprod)
{
  int j;
  DdNode *temp, *belief;
  // first we have to compute the whole belief state
  temp = bayesianUpdate(b,a);
  belief = multBeliefObsProd(temp,obsprod);
  Cudd_RecursiveDeref(gbm,temp);
  return belief;
}
DdNode * POMDP::bayesianUpdate(DdNode *b, int a, DdNode *obsprod)
{
  int j;
  DdNode *temp, *belief;
  // first we have to compute the whole belief state
  temp = bayesianUpdate(b,a);
  belief = multBeliefObsProd(temp,obsprod);
  Cudd_RecursiveDeref(gbm,temp);
  return belief;
}
// new, correct, sampling method
// iterates through the cubes
DdNode * POMDP::sample_Belief(DdNode *b)
{
  int i;
  // the sample random value
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  MixGauss *theMixGauss;
  int *thecube;
  int *thevals = new int[numvars];
  DdGen * theGen;
  theGen = Cudd_FirstCube(gbm,b,&thecube,&theMixGauss);
  int gen=1;
  double sum(0.0);
  double factor = 1;
  while (gen) {
    // see how many don't care's there are
    factor = 1;
    for (i=0; i<numvars; i++)
      if (thecube[2*i] == 2) 
	factor *= 2;
    sum += factor*(theMixGauss->get_max());
    if (sum >= rnd) {
      // copy out the cube values
      for (i=0; i<numvars; i++) 
	thevals[i] = thecube[2*i];
      // free the generator
      gen = Cudd_GenFree(theGen);
    } else {
      gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
    }
  }
  // now build the cube from the current values in 'thecube'
  DdNode **thevars =new DdNode*[numvars]; //(DdNode **) malloc(numvars*sizeof(DdNode *));
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
  delete [] thevars; //free(thevars);
  delete [] thevals;
  return theres;

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
// merges two sorted arrays of doubles 
// replaces szc with the new, sorted array
// replaces snr with the total number of elements in szc (the new sorted array)
// removes any duplicates
void merge(const double *zc, const int nr, double **szc, int & snr)
{
  double *newzc = new double[nr+snr];
  int k(0),l(0),j(0);
  while (k < nr || l < snr) {
    if ((k >= nr && l < snr) || ((l<snr) && ((*szc)[l] < zc[k]))) 
      newzc[j++] = (*szc)[l++];
    if ((l>= snr && k < nr) || ((k < nr) && (zc[k] < (*szc)[l]))) 
      newzc[j++] = zc[k++];
    if (k < nr && l < snr && zc[k] == (*szc)[l]) {
      newzc[j++] = zc[k++];
      l++;
    } 
  }
  snr = j;
  k=0;
  // remove duplicates
  double tmp = newzc[k];
  j=0;
  while (k < snr) {
    k++;
    while(k<snr && fabs(newzc[k]-tmp) < 0.0001) 
      k++;
    j++;
    if (k < snr) {
      tmp = newzc[k];
      newzc[j] = tmp;
    }
  }
  snr = j;
  delete [] *szc;
  *szc = new double[snr];
  for (k=0; k<snr; k++) 
    (*szc)[k] = newzc[k];
  delete [] newzc;
}
// partitions the observation space according to the alpha vectors at pbelief
// then integrates the alpha vectors over these partitions
DdNode *POMDP::partitionIntegrateObsSpace(int pbelief, int a) 
{
  // First,  compute an array of DdNodes over observation: the gamma_ao vectors at pbelief
  int i,k,l;
  double maxval, midp, val;
  DdNode **intObsFuns = new DdNode *[numalphas];
  DdNode *intval,*temp,*tempM;

  // compute Regions in which each mgs is maximal 
  // and integrate the observation functions over the regions for each alpha
  // giving the array intObsFuns
  computeIntegrateRegions(totalObsFun[a], intObsFuns, pbelief, a);
  
  // backup alpha vectors using aggregated observation functions
  // gives the new alpha vector
  tempM = backupAlphas(intObsFuns, NewPrime[a]);

  for (i=0; i<numalphas; i++) 
    Cudd_RecursiveDeref(gbm,intObsFuns[i]);
  delete [] intObsFuns;
  return tempM;
}
// computes the regions over which each mg[0..numalphas-1] is maximal
// this has to happen over each of the branches of obs variable
// integrates the observation functions over those regions for each alpha
// 
void POMDP::computeIntegrateRegions(DdNode **obsfun, DdNode **intObsFuns, int b, int a)
{
  int i,o;
  MixGauss **mgs = new MixGauss*[numalphas];
  MixGauss **regions = new MixGauss*[numalphas];
  char fname[256];

  DdNode *tmp1,*tmp2;
  for (i=0; i<numalphas; i++) {
    intObsFuns[i] = One;
    Cudd_Ref(intObsFuns[i]);
  }
  
  // go over each branch of 'obs'
  for (o=0; o<numorigobs; o++) {
    // compute value of each alpha vector for this observation variable at the forward propagated belief states
    for (i=0; i<numalphas; i++) {
      tmp1 = computeValue(alphas[i], baz[a][o][b]);
      mgs[i] = new MixGauss(*(Cudd_V(tmp1)));
      Cudd_RecursiveDeref(gbm,tmp1);
    }
    // call a computeregions function on all the mgs
    partitionObservationSpace(orig_obs[o].type,mgs,regions);

    // integrate observation function obsfun over those regions
    // and keep a running product over numorigobs for each alpha
    for (i=0; i<numalphas; i++) {
      tmp1 = Cudd_addConst(gbm,regions[i]);
      Cudd_Ref(tmp1);
      tmp2 = My_addApply(gbm,ddintegrate,obsfun[o],tmp1);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp1);
      // this should be saved as intObsFun[o][i] for each observation+alpha pair
      // this product is the problem
      tmp1 = Cudd_addApply(gbm,Cudd_addTimes,tmp2,intObsFuns[i]);
      Cudd_Ref(tmp1);
      Cudd_RecursiveDeref(gbm,tmp2);
      Cudd_RecursiveDeref(gbm,intObsFuns[i]);
      intObsFuns[i] = tmp1;
    }
    for (i=0; i<numalphas; i++) {
      delete mgs[i];
      delete regions[i];
    }
    // possibly write out the obsimage for each 2D observation function
    /*
    if (lastiteration && orig_obs[o].type == 2) {
      sprintf(fname,"/h/23/jhoey/project/pomdp/obsreg%d_%d_%d_%d.dat",counter,currBeliefObsImage,a,o); 
      writeObsImageRaw(fname);
      sprintf(fname,"/h/23/jhoey/project/pomdp/obsreg%d_%d_%d_%d.mgs",counter,currBeliefObsImage,a,o); 
      writeObsImageMgs(fname);
    }
    */
  }
  delete [] mgs;
  delete [] regions;
}

// backs up the alpha vectors using the obsfuns for each one that have
// been aggregated over the regions
DdNode *POMDP::backupAlphas(DdNode **intObsFun, DdNode **transition)
{
  int i;
  DdNode *temp, *temp2, *intval, *ofintval;

  DdNode **alphap;
  alphap = new DdNode *[numalphas];

  // first change all variables to primes in alpha(s) -> alpha(s')
  for (i=0; i<numalphas; i++) {
    alphap[i] = Cudd_addSwapVariables(gbm,alphas[i],Array1,Array2,numvars);
    Cudd_Ref(alphap[i]);
  }

  intval = Zero;
  Cudd_Ref(intval);
  for (i=0; i<numalphas; i++) {
    // multiply integrated Observation Functionby alpha vector
    temp = Cudd_addApply(gbm,Cudd_addTimes,intObsFun[i],alphap[i]);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,alphap[i]);

    // then by the transition funciton and sum out primed variable
    temp2 = multiplySumSet(temp,transition);
    Cudd_RecursiveDeref(gbm,temp);

    // finally, multiply by the discount
    temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,discount);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp2);

    // add to sum
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp,intval);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,intval);
    intval = temp2;
  }
  delete [] alphap;
  return intval;
}
DdNode * POMDP::regionIntegrateBackup(int b, int a)
{
  int i,o,j;
  char fname[256];
  DdNode *tmp1,*tmp2;
  DdNode *temp, *temp2, *intval, *alphasum;
  DdNode **alphap = new DdNode *[numalphas];
  MixGauss **mgs = new MixGauss*[numalphas];
  MixGauss ***regions = new MixGauss**[numorigobs];
  for (o=0; o<numorigobs; o++) 
    regions[o] = new MixGauss*[numalphas];

  // first change all variables to primes in alpha(s) -> alpha(s')
  for (i=0; i<numalphas; i++) {
    alphap[i] = Cudd_addSwapVariables(gbm,alphas[i],Array1,Array2,numvars);
    Cudd_Ref(alphap[i]);
  }
    // compute value of each alpha vector for each  observation variable at the forward propagated belief states
  for (o=0; o<numorigobs; o++) {
    for (i=0; i<numalphas; i++) {
      tmp1 = computeValue(alphas[i], baz[a][o][b]);
      mgs[i] = new MixGauss(*(Cudd_V(tmp1)));
      Cudd_RecursiveDeref(gbm,tmp1);
    }
    // call a computeregions function on all the mgs
    partitionObservationSpace(orig_obs[o].type,mgs,regions[o]);
    for (i=0; i<numalphas; i++) 
      delete mgs[i];
  }


  intval = Zero;
  Cudd_Ref(intval);
  for (i=0; i<numalphas; i++) {
    alphasum = Zero;
    Cudd_Ref(alphasum);
    // for each primed variable separately
    for (j=0; j<numorigvars; j++) {
      // multiply observation function for this variable
      // here, we need to get the observation o generated from primed variable j
      // for now just assume its the same
      o = j;
      temp = Cudd_addApply(gbm,Cudd_addTimes,totalObsFun[a][o],alphap[i]);
      Cudd_Ref(temp);

      // then by the transition function and sum out the primed variable
      temp2 = multiplySum(temp,NewPrime[a][j],prime_vars[j].add_var,orig_vars+j);
      Cudd_RecursiveDeref(gbm,temp);

      // integrateover the regions for this alpha
      tmp1 = Cudd_addConst(gbm,regions[o][i]);
      Cudd_Ref(tmp1);
      tmp2 = My_addApply(gbm,ddintegrate,temp2,tmp1);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp1);
      Cudd_RecursiveDeref(gbm,temp2);

      // and keep a running sum for each alpha
      tmp1 = Cudd_addApply(gbm,Cudd_addPlus,tmp2,alphasum);
      Cudd_Ref(tmp1);
      Cudd_RecursiveDeref(gbm,tmp2);
      Cudd_RecursiveDeref(gbm,alphasum);
      alphasum = tmp1;
    }
    // add to sum
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,alphasum,intval);
    Cudd_Ref(temp2);
    Cudd_RecursiveDeref(gbm,alphasum);
    Cudd_RecursiveDeref(gbm,intval);
    intval = temp2;
  }
  // finally, multiply by the discount
  temp = Cudd_addApply(gbm,Cudd_addTimes,intval,discount);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,intval);
  intval = temp;

    // possibly write out the obsimage for each 2D observation function
  if (lastiteration && orig_obs[o].type == 2) {
    sprintf(fname,"/h/23/jhoey/project/pomdp/obsreg%d_%d_%d_%d.dat",counter,currBeliefObsImage,a,o); 
    writeObsImageRaw(fname);
    sprintf(fname,"/h/23/jhoey/project/pomdp/obsreg%d_%d_%d_%d.mgs",counter,currBeliefObsImage,a,o); 
    writeObsImageMgs(fname);
  }
  for (o=0; o<numorigobs; o++) {
    for (i=0; i<numalphas; i++) 
      delete regions[o][i];
    delete [] regions[o];
  }
  delete [] regions;
  delete [] mgs;
  delete [] alphap;
  return intval;
}
// partitions the observation space into regions where a single mgs is maximal
// the result is returned in an array of regions MixGauss'. 
// for 0-D
// regions[i]'s mixweight[j]
// give whether or not observation j is best at this mg (its either 1 or 0)
// for 1-D
// the result is returned in an array of regions MixGauss'. regions[i]'s mixweights
// give the regions over which mgs[i] is maximal (in pairs)
// the value gives the integration type
// if value = 0 each pair is normal
// if value = 1 the first pair is from -inf to second element of pair
// if value = 2 the last pair is from first element of pair to +inf
// if value = 3 the first pair is from -inf to inf
// for 2-D
// the result is returned in an array of regions MixGauss'. regions[i]'s mixweights
// give the pixels indices for which mgs[i] is maximal (in pairs)
// now gives the index only - use the global array obsimage to store everything else
void POMDP::partitionObservationSpace(int type, MixGauss **mgs, MixGauss **regions)
{
  int j,i,o,k,l;
  int nobs = mgs[0]->nmix;
  double maxval,val;
  int snr,numr,maxmg;
  double *bounds, *maxmgs, *regs;
  if (type == 0) {
    if (nobs == 0 && mgs[0]->val == 0.0) {
      // this is just 0, but its really for all observations
      for (i=0; i<numalphas; i++) 
	regions[i] = new MixGauss(1.0);
    } else {
      maxmgs = new double[nobs];
      for (o=0; o<nobs; o++) {
	for (i=0; i<numalphas; i++) {
	  if (i==0 || mgs[i]->pdf(o) > maxval) {
	    maxval = mgs[i]->pdf(o);
	    maxmgs[o] = i;
	  }
	}
      }
      double *v = new double[nobs];
      for (i=0; i<numalphas; i++) {
	regions[i] = new MixGauss(0,nobs,v);
	for (o=0; o <nobs; o++) 
	  regions[i]->mixweights[o] = (maxmgs[o] == i);
      }
      delete [] v;
    }
  } else if (type == 1) {
    partitionObservationSpace(mgs,snr,bounds,maxmgs);
    for (i=0; i<numalphas; i++) {
      // count number of regions
      numr = 0;
      for (k=0; k<snr-1; k++) {
	if (maxmgs[k] == i)
	  numr++;
      }
      if (numr > 0) {
	regs = new double[2*numr];
	l=0;
	for (k=0; k<snr-1; k++) {
	  if (maxmgs[k] == i) {
	    regs[l++] = bounds[k];
	    regs[l++] = bounds[k+1];
	  }
	}
	regions[i] = new MixGauss(0,numr*2,regs);
	// set the integration type for these regions
	// val = 0 means each region is normal 
	// val = 1 means first is really from -inf to regs[1]
	// val = 2 means last is really from regs[numr*2-2] to +inf
	// val = 3 means first is really from -inf to +inf
	regions[i]->val = 0;
	// if first region belongs to this alpha - then its from -inf
	if (maxmgs[0] == i) {
	  regions[i]->val = 1;
	} 
	// if last region belongs to this alpha - then its to +inf
	if (maxmgs[snr-2] == i) {
	  regions[i]->val = 2;
	}
	// if first and last regions belong to this alpha then its -inf and +inf 
	if (maxmgs[0] == i && maxmgs[snr-2] == i) {
	  regions[i]->val = 3;
	}
	delete [] regs;
      } else {
	regions[i] = new MixGauss(0.0);
      }
    }
    delete [] bounds;
    delete [] maxmgs;
  } else if (type == 2) {
    j=0;
    if (!dosampling) {
      // a regular grid (every pixel
      for (i=0; i<totnumsamples; i++) {
	// but for now, we just let it be
	obssamples[j]->z[0] = j%obsimage.nx;
	obssamples[j]->z[1] = j/obsimage.nx;
	j++;
      }
    } else {
     
      /*
    for (l=0; l<numalphas; l++) {
      numsamples[l] = MIN(numsamplespermg,(totnumsamples-tmp));
      tmp += numsamples[l];
      for (i=0; i<numsamples[l]; i++) {
	// normally, this would be like this:
	if (mgs[l].drawSample(obssamples[j]->z)) {
	  obssamples[j]->z[0] = floor(obssamples[j]->z[0]);
	  obssamples[j]->z[1] = floor(obssamples[j]->z[1]);
	  if (obssamples[j]->z[0] >= obsimage.nx)
	    obssamples[j]->z[0] = obsimage.nx-1;
	  if (obssamples[j]->z[0] < 0)
	    obssamples[j]->z[0] = 0;
	  if (obssamples[j]->z[1] >= obsimage.ny)
	    obssamples[j]->z[1] = obsimage.ny-1;
	  if (obssamples[j]->z[1] < 0)
	    obssamples[j]->z[1] = 0;
	} else {
	  // this means that the mgs is Zero 
	  // so we need to sample uniformly on a grid
	  obssamples[j]->z[0] = ((int) floor(j*samplefac))%obsimage.nx;
	  obssamples[j]->z[1] = ((int) floor(j*samplefac))/obsimage.nx;
	}
	j++;
      }
    }
    */
    }
    //fprintf(stderr,"drew %d samples\n",j);
    for (j=0; j<totnumsamples; j++) {
      //fprintf(stderr,"%f %f\n",obssamples[j]->z[0],obssamples[j]->z[1]);
      for (l=0; l<numalphas; l++) {
	val = mgs[l]->pdf(obssamples[j]->z);
	if (l==0 || val > maxval) {
	  maxval = val;
	  maxmg = l;
	}
      }
      obssamples[j]->index = maxmg;
      obssamples[j]->val = maxval;
    }
    for (i=0; i<numalphas; i++) 
      regions[i] = new MixGauss(i+1);
  }
}
// partitions the observation space into regions where a single mgs is maximal
// the result is returned in a double array in bounds, and the maximal mg between each bound is returned in maxmgs;
// thus, maxmgs[k] gives the index of the best mg in between bounds[k] and bounds[k+1]
// snr is the total number of bounds found (including -inf and +inf)
// this is for 1D
void POMDP::partitionObservationSpace(MixGauss **mgs, int & snr, double * & bounds, double * &maxmgs)
{
  int i,k,l,nr;
  double maxval, midp, val;
  double *zc;// = new double[256];
  double *tszc = new double[1024];
  snr = 0;
  
  for (k=0; k<numalphas; k++) {
    // compare kth to all others 
    for (l=k+1; l<numalphas; l++) {
      nr = mgs[k]->intersections(*(mgs[l]), &zc);
      merge(zc,nr,&tszc,snr);
    }
  }
  /*
  for (k=0; k<numalphas; k++) {
    fprintf(stderr,"%s\n",mgs[k]->toString());
  }
  fprintf(stderr,"%d intersections are :",snr);
  for (k=0;k<snr; k++) 
    fprintf(stderr,"%g ",tszc[k]);
  fprintf(stderr,"\n");

  // do it again using old version
  snr = 0;
  for (k=0; k<numalphas; k++) {
    // compare kth to all others 
    for (l=k+1; l<numalphas; l++) {
      nr = mgs[k]->oldintersections(*(mgs[l]), &zc);
      merge(zc,nr,&tszc,snr);
    }
  }
  fprintf(stderr,"%d intersections are :",snr);
  for (k=0;k<snr; k++) 
    fprintf(stderr,"%g ",tszc[k]);
  fprintf(stderr,"\n");
  */
  

  //The regions do not contain the
  // end points -inf & +inf. These must be added now
  bounds = new double[snr+2];
  maxmgs = new double[snr+1];
  if (snr > 0) {
    bounds[0] = tszc[0]-1.0;
  } else {
    bounds[0] = 0.0;
  }
  for (k=1; k<=snr; k++) 
    bounds[k] = tszc[k-1];
  if (snr > 0) {
    bounds[snr+1] = tszc[snr-1]+1.0;
  } else {
    bounds[snr+1] = 0;
  }
  snr += 2;
  for (k=0; k<snr-1; k++) {
    // look in region bounds[k],bounds[k+1]
    // figure out the midpoint,
    midp = (bounds[k+1]-bounds[k])/2.0+bounds[k];

    //find maximal mg
    for (l=0; l<numalphas; l++) {
      val = mgs[l]->pdf(&midp);
      if (l==0 || val > maxval) {
	maxval = val;
	maxmgs[k] = l;
      }
    }
  }
  if (numalphas >= 2) 
    delete [] zc;
  delete [] tszc;
}
// to be used with Cudd_addApply. The function f will have its
// leaves integrated over (and hence replaced by MixGauss' with just
// a val in them 
// g gives the range of integration in the first two weights 
// (its a MixGauss with nmix = 2 and mixweights[0] and mixweights[1] 
// give the integration range. 
// the typ of the integration (0 = normal, 1 = from -inf, 2 = to +inf, 3 = -inf:+inf)
// is given by the val of the G MixGauss
DdNode *ddintegrate(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *G;
  
  F = *f; G = *g;
  
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    if (F == One) {
      return One;
    }
    if (G == Zero) {
      return Zero;
    }
    double sum = integrate(*(Cudd_V(F)), *(Cudd_V(G)));
    MixGauss newMG(sum);
    DdNode *result = Cudd_addConst(gbm,&newMG);
    Cudd_Ref(result);
    return result;
  } else {
    return NULL;
  }
}

// integrates this mg over the region given by region
// works for 0, 1 or 2D functions
//
// for 0-D
// regions[i]'s mixweight[j]
// give whether or not observation j is best at this mg (its either 1 or 0)
// for 1-D
// the result is returned in an array of regions MixGauss'. regions[i]'s mixweights
// give the regions over which mgs[i] is maximal (in pairs)
// the value gives the integral type as described in POMDP.cpp
// for 2-D
// the result is returned in an array of regions MixGauss'. regions[i]'s mixweights
// give the pixels indices for which mgs[i] is maximal (in pairs)
// now its just the index of the mg
double integrate(MixGauss & integrand, MixGauss & region)
{
  int i;
  double sum(0);
  int type = integrand.fvdim;
  if (type == 0) {
    if (region.nmix == 0) {
      sum = 1;
    } else {
      for (i=0; i<integrand.nmix; i++) {
	if (region.mixweights[i] != 0) 
	  sum += integrand.mixweights[i];
      }
    }
  } else if (type == 1) {
    i = 0;
    double lo, hi;
    // the type of integration  (0 = normal, 1 = from -inf, 2 = to +inf, 3 = -inf:+inf)
    int integraltype, ityp;
    integraltype = (int) (region.val);
    i=0;
    while (i < region.nmix-1) {
      lo = region.mixweights[i];
      hi = region.mixweights[i+1];
      if (integraltype == 0) {
	//vanilla - all regions are good
	sum += integrand.integrate(integraltype, lo, hi);
      } else if (integraltype == 1 || integraltype == 3) {
	if (i==0 && region.nmix == 2) {
	  sum += integrand.integrate(integraltype,lo,hi);
	} else if (i==0) {
	  sum += integrand.integrate(1,lo,hi);
	} else {
	  sum += integrand.integrate(0,lo,hi);
	}
      } else if (integraltype == 2 || integraltype == 3) {
	if (i < region.nmix-2) {
	  sum += integrand.integrate(0,lo,hi);
	} else {
	  sum += integrand.integrate(2,lo,hi);
	}
      }
      i+=2;
    }
  } else if (type == 2) {
    for (i=0; i<totnumsamples; i++) {
      if (obssamples[i]->index == (region.val-1)) {
	sum += integrand.pdf(obssamples[i]->z);
      }
    }
  } 
  return sum;
}
DdNode *computePDF(DdManager *dd, DdNode **f, DdNode **g)
{
    DdNode *F, *G;
    DdNode *result;
  F = *f; G = *g;
  double res;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    if (G == Zero || F == Zero) {
      return Zero;
    }
    if (G == One || F == One) {
      return One;
    }
    res = (*Cudd_V(F)).pdf(*Cudd_V(G));
    if (res < 1e-10) {
      res = 1e-10;
    }
    MixGauss newMG(res);
    result = Cudd_addConst(gbm,&newMG);
    Cudd_Ref(result);
    return result;
  } else {
    return NULL;
  }
}
// only for fvdim=0 mixtures
DdNode *computePDFType0(DdManager *dd, DdNode **f, DdNode **g)
{
    DdNode *F, *G;
    DdNode *result;
  F = *f; G = *g;
  double res;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    if (F == Zero) {
      return Zero;
    }
    res = (*Cudd_V(F)).pdf(*Cudd_V(G));
    MixGauss newMG(res);
    result = Cudd_addConst(gbm,&newMG);
    Cudd_Ref(result);
    return result;
  } else {
    return NULL;
  }
}
// to be used with addApply
// draws a sample from the observation branch if g==1 (g is a cube)
DdNode *drawSample(DdManager *dd, DdNode **f, DdNode **g)
{
  DdNode *F, *G;
  F = *f; G = *g;
  DdNode *temp, *result;
  int oo;
  if (Cudd_IsConstant(G)) {
    if (G==Zero) {
      return Zero; 
    } else if (Cudd_IsConstant(F)) {
      if (F==One) {
	return One;
      }
      if (F==Zero) {
	return Zero;
      }
      // now, F is a constant over some observation, so call its sample function. 
      MixGauss *mg = Cudd_V(F);
      MixGauss newMG;
      if (mg->fvdim == 0) {
	//newMG = new MixGauss(mg->drawMixtureSample());
	newMG.set(mg->drawMixtureSample());
      } else {
	//newMG = new MixGauss(0, mg->fvdim, mg->drawSample());
	newMG.set(0,mg->fvdim,mg->drawSample());
      }
      DdNode *result = Cudd_addConst(gbm,&newMG);
      Cudd_Ref(result);
      return result;
    }
  }
  return NULL;
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
  temp = sumOutAllPrimedOrigVars(temp1);
  Cudd_RecursiveDeref(gbm,temp1);
  temp1 = temp;

  if (Cudd_IsConstant(temp1))
    res = (*(Cudd_V(temp1))).get_min();
  else {
    fprintf(stderr,"whoops - something wrong!\n");
    return -1;
  }
  return res;
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
// unfactored version
void POMDP::printBelief(DdNode *b, FILE *fd) {
  DdNode **tmp = beliefMarginal(b);
  printBelief(tmp,fd);
}
void POMDP::printBelief(DdNode **b, FILE *fd) {
  int j,i;
  DdNode *temp, *temp1;
  double tmpval;
  for (j=0; j<numorigvars; j++) {
    fprintf(fd,"%s ",orig_vars[j].name);
    for (i=0; i<orig_vars[j].nvals; i++) {
      temp1 = restrictVal(gbm,b[j],prime_vars,numvars,orig_vars,numorigvars,j,i);
      Cudd_Ref(temp1);
      //tmpval = 0.01*((int) (100*((*(Cudd_V(temp1))).get_max())));
      tmpval = (*(Cudd_V(temp1))).get_max();
      fprintf(fd," %s: %3.2f",orig_vars[j].valname[i],tmpval);
      Cudd_RecursiveDeref(gbm,temp1);
    }
    fprintf(fd,"\n");
  }
}
void POMDP::printPolicy(int *policy, int numbeliefs, DdNode ***beliefs, double *beliefValues)
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
    fprintf(stderr," value %f\n",beliefValues[b]);
    printBelief(beliefs[b]);
  }
}
#ifdef PRINTALLACTS
void POMDP::printPolicyHTML(char *fname, int *policy, int numbeliefs, DdNode ***beliefs, double *beliefValues, int hor, int ov)
{
  int b,a;
  int tmp;
  FILE *fd = fopen(fname,"w");
  fprintf(fd,"<HTML>\n<BODY>\n");

  for (b=0; b<numbeliefs; b++) {
    fprintf(fd,"<P><HR>\n<PRE>\n");
    fprintf(fd,"\nbest actions for belief %d are ",b);
    tmp = policy[b];
    a = 0;
    while (tmp >= 1) {
      if (tmp%2) 
	fprintf(fd," %s",actionlist[a].name);
      a++;
      tmp = tmp/2;
    }
    fprintf(fd," value %f\n",beliefValues[b]);
    printBelief(beliefs[b],fd);
    fprintf(fd,"</PRE>\nAlpha Values:\n");
    fprintf(fd,"<TABLE>\n<TR>\n");
    for (a=0; a<numactions; a++) 
      fprintf(fd,"<TD>%s</TD>",actionlist[a].name);
    fprintf(fd,"</TR>\n<TR>\n");
    for (a=0; a<numactions; a++) 
      fprintf(fd,"<TD><IMG SRC='obsreg%d_%d_%d.gif'></IMG></TD>\n",b,a,ov);
    fprintf(fd,"</TR>\n</TABLE>\n");
    fprintf(fd,"\nConditional Plans:\n");
    fprintf(fd,"<TABLE>\n<TR>\n");
    for (a=0; a<numactions; a++) 
      fprintf(fd,"<TD>%s</TD>",actionlist[a].name);
    fprintf(fd,"</TR>\n<TR>\n");
    for (a=0; a<numactions; a++) 
      fprintf(fd,"<TD><IMG SRC='condplan%d_%d_%d.gif'></IMG></TD>\n",b,a,ov);
    fprintf(fd,"</TR>\n</TABLE>\n");
  }
  fprintf(fd,"</BODY>\n</HTML>\n");
  fclose(fd);
}
#else
void POMDP::printPolicyHTML(char *fname, int *policy, int numbeliefs, DdNode ***beliefs, double *beliefValues, int hor, int ov)
{
  int b,a;
  int tmp;
  FILE *fd = fopen(fname,"w");
  fprintf(fd,"<HTML>\n<BODY>\n");

  for (b=0; b<numbeliefs; b++) {
    fprintf(fd,"<P><HR>\n<PRE>\n");
    fprintf(fd,"\nbest actions for belief %d are ",b);
    tmp = policy[b];
    a = 0;
    while (tmp >= 1) {
      if (tmp%2) 
	fprintf(fd," %s",actionlist[a].name);
      a++;
      tmp = tmp/2;
    }
    fprintf(fd," value %f\n",beliefValues[b]);
    printBelief(beliefs[b],fd);
    fprintf(fd,"</PRE><TABLE>\n<TR>\n");
    fprintf(fd,"<TD>Alpha Values</TD><TD>Conditional Plan</TD></TR>\n<TR>\n");
    fprintf(fd,"<TD><IMG SRC='obsreg%d_%d.gif'></IMG></TD>\n",b,ov);
    fprintf(fd,"<TD><IMG SRC='condplan%d_%d.gif'></IMG></TD>\n",b,ov);
    fprintf(fd,"</TR>\n</TABLE>\n");
  }
  fprintf(fd,"</BODY>\n</HTML>\n");
  fclose(fd);
}
#endif
void POMDP::printAlphas()
{ 
  int a;
  for (a=0; a<numalphas; a++) {
    fprintf(stderr,"alpha %d\n",a);
    Cudd_PrintDebug(gbm,alphas[a],1,2);
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
    MixGauss res;
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
  MixGauss nfp;
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
  if (primedvars)
    temp2 = sumOutAllPrimedOrigVars(temp1);
  else
    temp2 = sumOutAllOrigVars(temp1);
  Cudd_RecursiveDeref(gbm,temp1);
  temp1 = temp2;

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
  MixGauss *theMixGauss;
  MixGauss *fMixGauss = new MixGauss();
  int *thecube;
  int *thevals = new int[numvars];
  DdGen * theGen;
  DdNode *oCube, *temp1, *factorD;
  DdNode *theSpan = Zero;
  Cudd_Ref(theSpan);
  if (primedvars != 0) 
    primedvars = 1;
  theGen = Cudd_FirstCube(gbm,f,&thecube,&theMixGauss);
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
    fMixGauss->set(factor);
    factorD = Cudd_addConst(gbm,fMixGauss);
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
    gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
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
    MixGauss res;
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
      MixGauss res;
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

// iim input nx*ny ppm (RGB) image
// oim output nx*ny pgm (binary) edge image
void braindeadEdgeDetect(int nx, int ny, int *iim, char *oim)
{
  int i,j,index,cindex;
  for (i=0; i<nx*ny; i++) 
    oim[i] = 0;
  for (i=1; i<ny; i++) {
    for (j=0; j<nx-1; j++) {
      index = (i*nx+j);
      if (iim[index] != iim[((i-1)*nx+j)] || 
	  iim[index] != iim[(i*nx+j+1)]) {
	oim[index] = 255;
      }
    }
  }
}



// Gets the Voronoi diagram for the observation function o (which must be type 2)
void POMDP::getEdgeImage2D(int a, int o, char *edgeimage) 
{
  int i,j,index,l,maxmg;
  double maxval,val,sum(0), z[2];
  int *thecube;
  MixGauss *theMixGauss;
  MixGauss **theMGs = new MixGauss*[256];
  DdGen * theGen;
  char fname[256];
  int gen=1, nmgs=0;

  // first, collect all the different observation functions
  theGen = Cudd_FirstCube(gbm,totalObsFun[a][o],&thecube,&theMixGauss);
  while (gen) {
    theMGs[nmgs] = new MixGauss(*theMixGauss);
    nmgs++;
    gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
  }
  Cudd_GenFree(theGen);

  for (i=0; i<obsimage.ny; i++) {
    for (j=0; j<obsimage.nx; j++) {
      index = i*obsimage.nx+j;
      z[0] = j; z[1] = i;
      maxval = -1.0;
      for (l=0; l<nmgs; l++) {
	val = theMGs[l]->pdf(z);
	if (val > maxval) {
	  maxval = val;
	  maxmg = l;
	}
      }
      // this maxmg is nothing to do with alpha vectors - 
      // just needs to be different for each mixture model
      obsimage.pixels[index] = maxmg;
      obsimage.vals[index] = maxval;
    }
  }
  for (i=0; i<nmgs; i++) 
    delete theMGs[i];
  delete [] theMGs;
  // get the edges 
  braindeadEdgeDetect(obsimage.nx,obsimage.ny,obsimage.pixels,edgeimage);
}
void writeEdgeImagePGM(char *filename, char *eimage)
{
  FILE *stream = fopen(filename, "w" );
  if ( !stream )
    return;
  fprintf( stream, "P5\n%d %d\n255\n", obsimage.nx, obsimage.ny);
  fwrite(eimage, sizeof(char), obsimage.ny*obsimage.nx, stream);
  fclose(stream);
}
void writeEdgeImagePPM(char *filename, char *imname, char *eimage)
{
  FILE *stream;
  char *im = new char[obsimage.nx*obsimage.ny*3];
  // open the image for reading
  stream = fopen(imname,"r");
  if ( !stream )
    return;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(im, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
  fclose(stream);
  // add edges
  for (int i=0; i<obsimage.nx*obsimage.ny; i++) {
    if (eimage[i] != 0) {
      im[3*i] = im[3*i+1] = im[3*i+2] = 0;
      im[3*i] = 255;
    }
  }

  stream = fopen(filename, "w" );
  if ( !stream ) {
    fprintf(stderr,"could not open %s\n",filename);
    return;
  }

  fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
  fwrite(im, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
  fclose(stream);
}
double * readObsImageRaw(char *filename, int &nx, int &ny)
{
  // read im as raw doubles
  FILE *stream = fopen(filename, "r" );
  if ( !stream )
    return NULL;
  int nnx, nny;
  double *vv;
  fscanf(stream, "%*s\n%d %d\n%*d\n", &nnx, &nny);
  nx = nnx;
  ny = nny;
  vv = new double[nx*ny];
  fprintf(stderr," reading %d by %d image\n",nx,ny);
  fread(vv, ny*nx*sizeof(double), 1, stream);
  fclose(stream);
  return vv;
}
void writeObsImageRaw(char *filename)
{
  int i;
  // first, collect obsimage from obssamples
  double maxval, minval;
  for (i=0; i<obsimage.nx*obsimage.ny; i++) {
    //obsimage.vals[i] = obssamples[i]->val;
    if (i==0 || obsimage.vals[i] > maxval) 
      maxval = obsimage.vals[i];
    if (i==0 || obsimage.vals[i] < minval) 
      minval = obsimage.vals[i];
  }
  fprintf(stderr,"max/min vals for obsimage are %g  %g\n",maxval,minval);
  // save im as raw doubles
  int fd;
  if ( (fd = creat(filename,0644)) == -1 ) {
    fprintf(stderr,"could not open %s\n",filename);
    return;
  }
  write(fd,obsimage.vals,obsimage.nx*obsimage.ny*sizeof(double));
  close(fd);
  /*
  int nx, ny;
  double *v = new double[160*120];
  fd = open(filename,O_RDONLY);
  if (fd < 0) {
    fprintf(stderr,"could not open %s\n",filename);
    return;
  }
  read(fd,v,obsimage.ny*obsimage.nx*sizeof(double));

  // look for differences
  // see if its still the same
  for (i=0; i<totnumsamples; i++) {
    if (i==0 || v[i] > maxval) 
      maxval = v[i];
    if (i==0 || v[i] < minval) 
      minval = v[i];
  }
  fprintf(stderr,"max/min vals for obsimage %s are %g  %g\n",filename,maxval,minval);
  */
  /*
  for (i=0; i<totnumsamples; i++) {
    if (obsimage.vals[i] != v[i]) {
      fprintf(stderr,"%g ",obsimage.vals[i]-v[i]);//holy crap! pixel %d is %g should be %g\n",i,v[i],obsimage.vals[i]);
    }
  }
  */
  //delete [] v;
}

void writeObsImageMgs(char *filename)
{
  // first, collect obsimage from obssamples
  //for (int i=0; i<totnumsamples; i++) 
  //obsimage.pixels[i] = obssamples[i]->index;
  // save im as raw ints
  int fd;
  if ( (fd = creat(filename,0644)) == -1 ) {
    fprintf(stderr,"could not open %s\n",filename);
    return;
  }
  write(fd,obsimage.pixels,obsimage.nx*obsimage.ny*sizeof(int));
  close(fd);
  /*
  FILE *stream = fopen(filename, "w" );
  if ( !stream )
    return;
  fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
  fwrite(obsimage.pixels, sizeof(int), obsimage.ny*obsimage.nx, stream);
  fclose(stream);
  */
}

// computes the product of the observation probabilities for observation z           
DdNode * POMDP::computeObsProb(DdNode *lbaz, double *z)
{                                                                                    
  DdNode *sampo,*temp2,*temp1,*obsprob;                                              
  MixGauss mgsamp(0,2,z);                                                                   
  sampo = Cudd_addConst(gbm,&mgsamp);                                              
  Cudd_Ref(sampo);                                                                 
  obsprob = Cudd_addApply(gbm,computePDF,lbaz,sampo);                        
  Cudd_Ref(obsprob);          
  Cudd_RecursiveDeref(gbm,sampo);
  return obsprob;                                                                    
} 
void POMDP::printObsImage(char *fname, char *imname, int o, int a) {
  int i,j,k, index;
  int tmp, ap;
  unsigned char pixel;
  FILE *stream;
  int colmap[][3] =   {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,1,1},{0,1,1},{1,0,1}};
  unsigned char *cpreg = new unsigned char[obsimage.nx*obsimage.ny*3];
  unsigned char *im = new unsigned char[obsimage.nx*obsimage.ny*3];
  // open the image for reading
  stream = fopen(imname,"r");
  if ( !stream )
    return;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(im, sizeof(unsigned char), 3*obsimage.ny*obsimage.nx, stream);
  fclose(stream);

  double gray[3] = {0.3,0.59,0.11}; 

  for (i=0; i<obsimage.ny; i++) {
    for (j=0; j<obsimage.nx; j++) {
      index = i*obsimage.nx+j;
      // greyscale pixel
      pixel = (unsigned char) MAX(0,MIN(255,(im[index*3]*gray[0]+im[index*3+1]*gray[1]+im[index*3+2]*gray[2])));
      for (k=0; k<3; k++) {
	if (edgeimage[a][index2Dobsfuns[o]][index] != 0) 
	  cpreg[index*3+k] = 255;
	else {
	  cpreg[index*3+k] = pixel*colmap[obsimage.pixels[index]][k]; 
	}
      }
    }
  }
  stream = fopen(fname, "w" );
  if ( !stream ) {
    fprintf(stderr,"could not open %s\n",fname);
    return;
  }
  fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
  fwrite(cpreg, sizeof(unsigned char), 3*obsimage.ny*obsimage.nx, stream);
  fclose(stream);
  
  delete [] cpreg;
}
// efficient version
void POMDP::printConditionalPlanRegions2(char *fprefix, int *policy, int numbeliefs, DdNode ***beliefs, int o)
{
  int *pb = new int[numbeliefs];
  for (int i=0; i<numbeliefs; i++) {
    pb[i] = i;
  }
  printConditionalPlanRegions2(fprefix,policy,numbeliefs,pb,beliefs,o);
}
void POMDP::printConditionalPlanRegions2(char *fprefix, int *policy, int numbeliefs, int *pb, DdNode ***beliefs, int o)
{
  int b,a,ap,i,j,l,maxmg,k,tmp,index,bestalpha;
  double beliefValue, val, maxval;
  double *z = new double[2];
  DdNode *pbelief, *temp, *tmp1,*bb, *sampo;
  char *cpreg = new char[obsimage.nx*obsimage.ny*3];
  FILE *stream;
  char fname[256];
  int colmap[][3] =   {{255,0,0},{0,255,0},{0,0,255},{255,255,0},{255,0,255},{0,255,255},{255,255,255},{0,0,0}};

  // write the edge images
  for (a=0; a<numactions; a++) {
    // write edge images
    sprintf(fname,"%s/edges%d.pgm",fprefix,a);
    writeEdgeImagePGM(fname,edgeimage[a][index2Dobsfuns[o]]);
  }

  MixGauss **mgs;
  bestalpha = 0;
  int policyact;
  // first compute the observation probabilities
  DdNode ***obsprobs = new DdNode**[numactions];
  z[0] = 0; z[1] = 0;
  MixGauss mgsamp(0,2,z);  
  int stepsize = 1;
  // modulation so all chosen pixels are not in columns
  int modu;
  for (a=0; a<numactions; a++) {
    fprintf(stderr,"computing obs probs for action %d\n",a);
    obsprobs[a] = new DdNode *[obsimage.nx*obsimage.ny];
    modu = 0;
    for (i=0; i<obsimage.ny; i+=stepsize) {
      for (j=modu; j<obsimage.nx; j+=stepsize) {
	index = i*obsimage.nx+j;
	z[0] = j; z[1] = i;
	mgsamp.set_weights(z);
	sampo = Cudd_addConst(gbm,&mgsamp);                                              
	Cudd_Ref(sampo);                                                                 
	obsprobs[a][index] = Cudd_addApply(gbm,computePDF,totalObsFun[a][o],sampo);
	Cudd_Ref(obsprobs[a][index]);
	Cudd_RecursiveDeref(gbm,sampo);
      }
      modu = (modu+1)%stepsize;
    }
  }
  // now go over beliefs
  for (b=0; b<numbeliefs; b++) {
    tmp = policy[pb[b]];
    policyact = 0;
    while (tmp >= 1 && !(tmp%2)) {
      policyact++;
      tmp = tmp/2;
    }
#ifdef PRINTALLACTS
    for (a=0; a<numactions; a++) {
#else
      a = policyact;
#endif
      bb = bayesianUpdate(beliefs[pb[b]],a);
      // now go over all observations and figure out
      // which alpha would be selected (and hence which action)
      // assuming only this observation was observed
      modu = 0;
      // zero the image
      for (i=0; i<obsimage.ny; i++) {
	for (j=0; j<obsimage.nx; j++) {
	  index = i*obsimage.nx+j;
	  obsimage.pixels[index] = 0;
	  obsimage.vals[index] = 0;
	  for (k=0; k<3; k++) {
	    if (edgeimage[a][index2Dobsfuns[o]][index] != 0) 
	      cpreg[index*3+k] = 255;
	    else
	      cpreg[index*3+k] = 0; 
	  }
	}
      }
      for (i=0; i<obsimage.ny; i+=stepsize) {
	for (j=modu; j<obsimage.nx; j+=stepsize) {
	  index = i*obsimage.nx+j;
	  z[0] = j; z[1] = i;
	  // get baz at this value of z
	  temp = Cudd_addApply(gbm,Cudd_addTimes,bb,obsprobs[a][index]);
	  Cudd_Ref(temp);
	  pbelief = NULL;
	  normalizeFunction(temp,&pbelief,true);
	  Cudd_RecursiveDeref(gbm,temp);
	  // get best action at that belief
	  // gets action ap
	  bestalpha = getPolicy(pbelief,beliefValue);
	  Cudd_RecursiveDeref(gbm,pbelief);
	  
	  // tmp is actions to take
	  tmp = bestactions[bestalpha];
	  // choose first one
	  ap = 0;
	  while (tmp >= 1 && !(tmp%2)) {
	    ap++;
	    tmp = tmp/2;
	  }
	  // assign output image a fancy color for each action
	  for (k=0; k<3; k++) {
	    if (edgeimage[a][index2Dobsfuns[o]][index] == 0) 
	      cpreg[index*3+k] = colmap[ap][k]; 
	  }
	  obsimage.pixels[index] = bestalpha+1; 
	  obsimage.vals[index] = beliefValue; 
	}
	modu = (modu+1)%stepsize;
      }
#ifdef PRINTALLACTS
    // write out cpreg to a ppm image
      sprintf(fname,"%s/condplan%d_%d_%d.ppm",fprefix,pb[b],a,o);
      stream = fopen(fname, "w" );
      if ( !stream ) {
	fprintf(stderr,"could not open %s\n",fname);
	return;
      }
      fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
      fwrite(cpreg, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
      fclose(stream);
      
      // write obsimage.pixels & obsimage.vals
      sprintf(fname,"%s/obsreg%d_%d_%d.dat",fprefix,pb[b],a,o); 
      writeObsImageRaw(fname);
      sprintf(fname,"%s/obsreg%d_%d_%d.mgs",fprefix,pb[b],a,o); 
      writeObsImageMgs(fname);
    }
#else
    // write out cpreg to a ppm image
    sprintf(fname,"%s/condplan%s%d_%d.ppm",fprefix,pb[b],o);
    stream = fopen(fname, "w" );
    if ( !stream ) {
      fprintf(stderr,"could not open %s\n",fname);
      return;
    }
    fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
    fwrite(cpreg, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
    fclose(stream);
    
    // write obsimage.pixels & obsimage.vals
    sprintf(fname,"%s/obsreg%d_%d.dat",fprefix,pb[b],o); 
    writeObsImageRaw(fname);
    sprintf(fname,"%s/obsreg%d_%d.mgs",fprefix,pb[b],o); 
    writeObsImageMgs(fname);
#endif
  }
  for (a=0; a<numactions; a++) {
    modu = 0;
    for (i=0; i<obsimage.ny; i+=stepsize) {
      for (j=modu; j<obsimage.nx; j+=stepsize) {
	index = i*obsimage.nx+j;
	Cudd_RecursiveDeref(gbm,obsprobs[a][index]);
      }
      modu = (modu+1)%stepsize;
    }
    delete [] obsprobs[a];
  }
  delete [] obsprobs;
  
  delete [] cpreg;
  delete [] z;
}
void POMDP::printConditionalPlanRegions3(char *fprefix, int acttotake, DdNode *belief, int o) 
{
  int b,a,ap,i,j,l,maxmg,k,tmp,index,bestalpha;
  double beliefValue, val, maxval;
  double *z = new double[2];
  DdNode *pbelief, *temp, *tmp1,*bb, *sampo;
  unsigned char *cpreg = new unsigned char[obsimage.nx*obsimage.ny*3];
  unsigned char *cpreg2 = new unsigned char[obsimage.nx*obsimage.ny*3];
  unsigned char *mgs1 = new unsigned char[obsimage.nx*obsimage.ny*3];
  unsigned char *mgs2 = new unsigned char[obsimage.nx*obsimage.ny*3];
  FILE *stream;
  char fname[256];
  double colmap[][3] =   {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{0.5,0.75,2}};
  double **colmapmgs;
  colmapmgs = new double*[numalphas];
  //figure out how many values/color there are
  int  valspc = (int) ceil(pow(10.0,log10((double) numalphas)/3));
  double shade = 1.0/((double) valspc);
  fprintf(stderr,"%d shades/per color - %f - to represent %d colors  %lf\n",valspc,shade,numalphas,log10((double) numalphas));
  
  for (i=0; i<numalphas; i++) {
    colmapmgs[i] = new double[3];
    for (k=0; k<3; k++) 
      colmapmgs[i][k] = 1.0;
  }
  for (i=1; i<numalphas; i++) {
    for (k=0; k<3; k++) 
      colmapmgs[i][k] = colmapmgs[i-1][k];
    colmapmgs[i][0] -= shade;
    if (colmapmgs[i][0] < 0.0) {
      colmapmgs[i][0] = 1.0;
      colmapmgs[i][1] -= shade;
      if (colmapmgs[i][1] < 0.0) {
	colmapmgs[i][1] = 1.0;
	colmapmgs[i][2] -= shade;
	if (colmapmgs[i][2] < 0.0) 
	  colmapmgs[i][2] = 1.0;
      }
    }
  }
  
  
  MixGauss **mgs;
  bestalpha = 0;
  int policyact;
  // first compute the observation probabilities
  DdNode **obsprobs;
  z[0] = 0; z[1] = 0;
  MixGauss mgsamp(0,2,z);  
  int stepsize = 1;

  // read in the background image
  unsigned char *uim = new unsigned char[obsimage.nx*obsimage.ny*3];
  stream  = fopen("/h/23/jhoey/papers/pomdp/sinks.PPM","r");
  if ( !stream )
    return;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(uim, sizeof(unsigned char), 3*obsimage.ny*obsimage.nx, stream);
  fclose(stream);


  // modulation so all chosen pixels are not in columns
  int modu;
  a = acttotake;
  fprintf(stderr,"computing obs probs for action %d\n",a);
  obsprobs = new DdNode *[obsimage.nx*obsimage.ny];
  modu = 0;
  for (i=0; i<obsimage.ny; i+=stepsize) {
    for (j=modu; j<obsimage.nx; j+=stepsize) {
      index = i*obsimage.nx+j;
      z[0] = j; z[1] = i;
      mgsamp.set_weights(z);
      sampo = Cudd_addConst(gbm,&mgsamp);                                              
      Cudd_Ref(sampo);                                                                 
      obsprobs[index] = Cudd_addApply(gbm,computePDF,totalObsFun[a][o],sampo);
      Cudd_Ref(obsprobs[index]);
      Cudd_RecursiveDeref(gbm,sampo);
    }
    modu = (modu+1)%stepsize;
  }
  
  // now go over all observations and figure out
  // which alpha would be selected (and hence which action)
  // assuming only this observation was observed
  modu = 0;
  // zero the image
  for (i=0; i<obsimage.ny; i++) {
    for (j=0; j<obsimage.nx; j++) {
      index = i*obsimage.nx+j;
      obsimage.pixels[index] = 0;
      obsimage.vals[index] = 0;
      for (k=0; k<3; k++) 
	cpreg[index*3+k] = 0; 
    }
  }
  double gray[3] = {0.3,0.59,0.11}; 
  double pixel;
  // get updates belief on action acttotake
  bb = bayesianUpdate(belief,a);
  fprintf(stderr,"getting conditional plan\n");
  for (i=0; i<obsimage.ny; i+=stepsize) {
    for (j=modu; j<obsimage.nx; j+=stepsize) {
      index = i*obsimage.nx+j;
      z[0] = j; z[1] = i;
      // get baz at this value of z
      temp = Cudd_addApply(gbm,Cudd_addTimes,bb,obsprobs[index]);
      Cudd_Ref(temp);
      pbelief = NULL;
      normalizeFunction(temp,&pbelief,true);
      Cudd_RecursiveDeref(gbm,temp);
      // get best action at that belief
      // gets action ap
      bestalpha = getPolicy(pbelief,beliefValue);
      Cudd_RecursiveDeref(gbm,pbelief);
      
      // tmp is actions to take
      tmp = bestactions[bestalpha];
      // choose first one
      ap = 0;
      while (tmp >= 1 && !(tmp%2)) {
	ap++;
	tmp = tmp/2;
      }
      // assign output image a fancy color for each action
      for (k=0; k<3; k++) {
	// greyscale pixel
	pixel = uim[index*3]*gray[0]+uim[index*3+1]*gray[1]+uim[index*3+2]*gray[2];
	cpreg[index*3+k] = (unsigned char) MAX(0,MIN(255,pixel*colmap[ap%7][k]));
	if (edgeimage[a][index2Dobsfuns[o]][index] != 0) 
	  cpreg2[index*3+k] = 255;
	else
	  cpreg2[index*3+k] = (unsigned char) MAX(0,MIN(255,255*colmap[ap%7][k]));

	// mgs with image in background - bright colors only - not anymore
	mgs1[index*3+k] = (unsigned char) MAX(0,MIN(255,pixel*colmapmgs[bestalpha][k])); 
	// only mgs - many colors
	mgs2[index*3+k] = (unsigned char) MAX(0,MIN(255,255*colmapmgs[bestalpha][k])); 
      }
      obsimage.pixels[index] = bestalpha+1; 
      obsimage.vals[index] = beliefValue; 
    }
    modu = (modu+1)%stepsize;
  }
  // write out cpreg to a ppm image
  sprintf(fname,"%scondplanim.ppm",fprefix);
  writeImagePPM(fname,cpreg,obsimage.nx,obsimage.ny);

  sprintf(fname,"%scondplan.ppm",fprefix);
  writeImagePPM(fname,cpreg2,obsimage.nx,obsimage.ny);

  sprintf(fname,"%smgim.ppm",fprefix);
  writeImagePPM(fname,mgs1,obsimage.nx,obsimage.ny);

  sprintf(fname,"%smg.ppm",fprefix);
  writeImagePPM(fname,mgs2,obsimage.nx,obsimage.ny);

  
  // write obsimage.pixels & obsimage.vals
  sprintf(fname,"%sobsreg0_%d.dat",fprefix,o); 
  writeObsImageRaw(fname);
  sprintf(fname,"%sobsreg0_%d.mgs",fprefix,o); 
  writeObsImageMgs(fname);

  Cudd_RecursiveDeref(gbm,bb);
  modu = 0;
  for (i=0; i<obsimage.ny; i+=stepsize) {
    for (j=modu; j<obsimage.nx; j+=stepsize) {
      index = i*obsimage.nx+j;
      Cudd_RecursiveDeref(gbm,obsprobs[index]);
    }
    modu = (modu+1)%stepsize;
  }
  delete [] obsprobs;
  delete [] uim;
  delete [] mgs1;
  delete [] mgs2;
  delete [] cpreg;
  delete [] z;
}
void writeImagePPM(char *filename, unsigned char *im, int nx, int ny)
{
  FILE *stream = fopen(filename, "w" );
  if ( !stream )
    return;
  fprintf( stream, "P6\n%d %d\n255\n", nx, ny);
  fwrite(im, sizeof(unsigned char), 3*ny*nx, stream);
  fclose(stream);
}


void POMDP::printConditionalPlanRegions(char *fprefix, int *policy, int numbeliefs, DdNode ***beliefs, int o)
{
  int b,a,ap,i,j,l,maxmg,k,tmp,index,bestalpha;
  double beliefValue, val, maxval;
  double *z = new double[2];
  DdNode *pbelief, *temp, *tmp1;
  char *cpreg = new char[obsimage.nx*obsimage.ny*3];
  FILE *stream;
  char fname[256];
  int colmap[][3] =   {{255,0,0},{0,255,0},{0,0,255},{255,255,0},{255,0,255},{0,255,255},{255,255,255},{0,0,0}};

  // write the edge images
  for (a=0; a<numactions; a++) {
    // write edge images
    sprintf(fname,"%s/edges%d.pgm",fprefix,a);
    writeEdgeImagePGM(fname,edgeimage[a][index2Dobsfuns[o]]);
  }

  // first, compute the forward propagated belief states if it hasn't already been done
  DdNode ****lbaz;
  // this is used for the old sampling version
  //#ifndef SAMPLEOBS
  //  lbaz = baz;
  //#else
  lbaz = new DdNode***[numactions];
  for (i=0; i<numactions; i++) {
    lbaz[i] = new DdNode**[numorigobs];
    for (j=0; j<numorigobs; j++) 
      lbaz[i][j] = new DdNode *[numbeliefs];
  }
  
  
  for (a=0; a<numactions; a++) {
    for (b=0; b<numbeliefs; b++) {
      temp = bayesianUpdate(beliefs[b],a);
      for (j=0; j<numorigobs; j++) {
	// multiply by observation function
	lbaz[a][j][b] = Cudd_addApply(gbm,Cudd_addTimes,totalObsFun[a][j],temp);
	Cudd_Ref(lbaz[a][j][b]);
      }
      Cudd_RecursiveDeref(gbm,temp);
    }
  }
  //#endif
  MixGauss **mgs;
  //z[0] = 0; z[1] = 0;
  //pbelief = computeObsProb(lbaz[0][o][0],z);
  bestalpha = 0;
  int policyact;
  for (b=0; b<numbeliefs; b++) {
    tmp = policy[b];
    policyact = 0;
    while (tmp >= 1 && !(tmp%2)) {
      policyact++;
      tmp = tmp/2;
    }
#ifdef PRINTALLACTS
    for (a=0; a<numactions; a++) {
#else
      a = policyact;
#endif
      // now go over all observations and figure out
      // which alpha would be selected (and hence which action)
      // assuming only this observation was observed
      for (i=0; i<obsimage.ny; i++) {
	for (j=0; j<obsimage.nx; j++) {
	  index = i*obsimage.nx+j;
	  z[0] = j; z[1] = i;
	  // get baz at this value of z
	  temp = computeObsProb(lbaz[a][o][b],z);
	  pbelief = NULL;
	  normalizeFunction(temp,&pbelief,true);
	  Cudd_RecursiveDeref(gbm,temp);
	  // get best action at that belief
	  // gets action ap
	  bestalpha = getPolicy(pbelief,beliefValue);
	  Cudd_RecursiveDeref(gbm,pbelief);
	  
	  // tmp is actions to take
	  tmp = bestactions[bestalpha];
	  // choose first one
	  ap = 0;
	  while (tmp >= 1 && !(tmp%2)) {
	    ap++;
	    tmp = tmp/2;
	  }
	  // assign output image a fancy color for each action
	  for (k=0; k<3; k++) {
	    if (edgeimage[a][index2Dobsfuns[o]][index] != 0) 
	      cpreg[index*3+k] = 255;
	    else
	      cpreg[index*3+k] = colmap[ap][k]; 
	  }
	  obsimage.pixels[index] = bestalpha; 
	  obsimage.vals[index] = beliefValue; 
	}
      }
#ifdef PRINTALLACTS
    // write out cpreg to a ppm image
      sprintf(fname,"%s/condplan%d_%d_%d.ppm",fprefix,b,a,o);
      stream = fopen(fname, "w" );
      if ( !stream ) {
	fprintf(stderr,"could not open %s\n",fname);
	return;
      }
      fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
      fwrite(cpreg, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
      fclose(stream);
      
      // write obsimage.pixels & obsimage.vals
      sprintf(fname,"%s/obsreg%d_%d_%d.dat",fprefix,b,a,o); 
      writeObsImageRaw(fname);
      sprintf(fname,"%s/obsreg%d_%d_%d.mgs",fprefix,b,a,o); 
      writeObsImageMgs(fname);
    }
#else
    // write out cpreg to a ppm image
    sprintf(fname,"%s/condplan%d_%d.ppm",fprefix,b,o);
    stream = fopen(fname, "w" );
    if ( !stream ) {
      fprintf(stderr,"could not open %s\n",fname);
      return;
    }
    fprintf( stream, "P6\n%d %d\n255\n", obsimage.nx, obsimage.ny);
    fwrite(cpreg, sizeof(char), 3*obsimage.ny*obsimage.nx, stream);
    fclose(stream);
    
    // write obsimage.pixels & obsimage.vals
    sprintf(fname,"%s/obsreg%d_%d.dat",fprefix,b,o); 
    writeObsImageRaw(fname);
    sprintf(fname,"%s/obsreg%d_%d.mgs",fprefix,b,o); 
    writeObsImageMgs(fname);
#endif
  }
  delete [] cpreg;
  delete [] z;
}
void extractPrintSamples(FILE *fd, DdNode *dd)
{
  MixGauss *theMixGauss;
  int *thecube;
  DdGen * theGen;
  theGen = Cudd_FirstCube(gbm,dd,&thecube,&theMixGauss);
  int gen=1;
  double sum(0.0);
  double factor = 1;
  while (gen) {
    // see how many don't care's there are
    theMixGauss->printValMixWeights(fd);
    gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
  }
  gen = Cudd_GenFree(theGen);
}
void POMDP::discretizeObservationFunction(int a, int o, image & oimage) 
{
  int i,j,k,index,l,maxmg;
  double maxval,val,sum(0), z[2];
  double **ppos;
  int *thecube;
  MixGauss *theMixGauss;
  MixGauss **theMGs = new MixGauss*[256];
  DdGen * theGen;
  char fname[256];
  int gen=1, nmgs=0;

  // first, collect all the different observation functions
  theGen = Cudd_FirstCube(gbm,totalObsFun[a][o],&thecube,&theMixGauss);
  while (gen) {
    theMGs[nmgs] = new MixGauss(*theMixGauss);
    nmgs++;
    gen = Cudd_NextCube(theGen,&thecube,&theMixGauss);
  }
  Cudd_GenFree(theGen);

  ppos = new double*[nmgs];
  for (l=0; l<nmgs; l++) {
    ppos[l] = new double[nmgs];
    for (k=0; k<nmgs; k++) 
      ppos[l][k] = 0.0;
  }
  for (i=0; i<oimage.ny; i++) {
    for (j=0; j<oimage.nx; j++) {
      index = i*oimage.nx+j;
      z[0] = j; z[1] = i;
      k = oimage.pixels[index];
      for (l=0; l<nmgs; l++) {
	ppos[l][k] += theMGs[l]->pdf(z);
      }
    }
  }
  for (l=0; l<nmgs; l++) {
    sum = 0.0;
    for (k=0; k<nmgs; k++) 
      sum += ppos[l][k];
    for (k=0; k<nmgs; k++) 
      ppos[l][k] = ppos[l][k]/sum;
  }
  for (l=0; l<nmgs; l++) {
    //fprintf(stderr,"%s \n\n",theMGs[l]->toString());
    for (k=0; k<nmgs; k++) 
      fprintf(stderr,"%f ",ppos[l][k]);
    //fprintf(stderr,"\n\n\n");
    fprintf(stderr,"\n");
  }
  for (i=0; i<nmgs; i++) {
    delete theMGs[i];
    delete [] ppos[i];
  }
  delete [] ppos;
  delete [] theMGs;
}
