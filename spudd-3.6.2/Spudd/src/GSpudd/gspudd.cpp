#include <sys/time.h>
#include <sys/times.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>

#include "POMDP.h"
int main(int argc, char *argv[])
{
  int counter,newcounter,approx_flag,flag,spuddvalue,numwrite,numread;
  int i,j,k,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int approxMeth, optimalFrom,limit_type,reorderMeth,reorderApplyMeth,shuffleFlag,dotFlag,dompi;
  int templist[MAXVARS];
  char templog[MAXLINE];
  char basepath[MAXLINE];
  char *argflag, *argval, *outputfile, *dualfile;
  char *infile, *outpath, *inpath, *alphafile, *belieffile;
  double prune;
  double global_max_error,max_error,tmp_error;
  int bigadd;
  /* Allocate storage for timestamps */
  long curtime;
  unsigned seed;
  double maxApp,minApp;

  limit_type = LIMIT_SIZE;   
  long int maxusememory=0,usememory=0;
  
  DdNode *cube,*temperror,*tempo,*temp11,*offset,*temp21,*tempround;
  DdNode *ApproximationError;
  DdNode *temp, *temp1,*temp2,*one,*crap,*tempM,*maxADD;

  FILE *debug,*value,*action,*errf,*both,*stats,*log,*log2,*reward,*actions;
  FILE *crapFILE, *dataFile, *dataFile2;
  
  // intiialize the dd manager
  /* yyin is the lex global var for input file */
  /* default values first*/
  tolerance = 0.1;
  bigadd = BIGADD;
  int Max_Size = MAX_SIZE;
  double Max_Error = MAX_ERROR;
  limit_type = LIMIT_SIZE;    /* can be LIMIT_SIZE or LIMIT_ERROR */
  reorderApplyMeth = REORDERAPPLY_ALL;        /* can be REORDERAPPLY_(FIRSTFIVE | ALL)*/
  reorderMeth = REORDER_NONE;        /* can be REORDER_(SIFT | RANDOM |MINSPAN | NONE)*/
  approxMeth = APPROX_NONE;
  optimalFrom = OPTIMAL_NONE;
  int toleranceMode = TOLERANCE_FIXED;
  shuffleFlag = 0;
  dotFlag = 0;
  dompi = 0;
  double mtol = 0;

  int nsamples(1),mnb,numtrials(10);

  outputfile = NULL;
  bool havealphas = false;
  bool havebeliefs = false;
  bool dosampling = false;
  bool verbose = false;
  bool funky = false;
  bool startfromalpha = false;
  if (argc >= 5) {
    infile = *++argv;
    outpath = *++argv;
    argc = argc-2;
    while (--argc > 0 && (*++argv)[0] == '-') {
      argflag = ++argv[0];
      if (strcmp(argflag,"mnb") == 0) {
	mnb = atoi(*++argv);
      } else if (strcmp(argflag,"nsam") == 0) {
	nsamples = atoi(*++argv);
	dosampling = true;
      } else if (strcmp(argflag,"af") == 0) {
	alphafile = *++argv;
	havealphas = true;
      } else if (strcmp(argflag,"bf") == 0) {
	belieffile = *++argv;
	havebeliefs = true;
      } else if (strcmp(argflag,"nst") == 0) {
	numtrials = atoi(*++argv);
      } else if (strcmp(argflag,"v") == 0) {
	verbose = true;
	argc++;
      } else if (strcmp(argflag,"funky") == 0) {
	funky = true;
	argc++;
      } else if (strcmp(argflag,"saf") == 0) {
	startfromalpha = true;
	argc++;
      }
      argc -= 1;
    }
  } else {
    fprintf(stderr,"usage: testmdp <infile> <outpath>  [-nsam <nsamples>] [-mnb <mnb>] [-bf <belieffile>]  [-af <alphafile>] [-nst <nt>] [-funky] [-v]\n");
    fprintf(stderr,"\t<infile> input file name\n");
    fprintf(stderr,"\t<outpath> output path\n");
    fprintf(stderr,"\t<nsamples> number of samples\n");
    fprintf(stderr,"\t<mnb> min. number of beliefs to sample\n");
    fprintf(stderr,"\t<belieffile> file to get belief samples from (mnb is ignored)\n");
    fprintf(stderr,"\t<alphafile> file with alpha vectors (if this argument specified, then value iteration not performed)\n");
    fprintf(stderr,"\t<nt> number of simulation trials (Default 10)\n");
    fprintf(stderr,"\tfunky signals to do combination discrete/continuous simulation\n");
    fprintf(stderr,"\tv signals to be verbose in simulation\n");
    fprintf(stderr,"\tsaf signals to start from specified alpha vectors (not simulate)\n");
    exit(0);
  }
  // also reads in the MDP
  POMDP *mdp = new POMDP(infile,dosampling,bigadd);

  // randomize from the timer
  struct timeval sd;
  int rseed;
  gettimeofday(&sd,NULL);
  rseed = (int) sd.tv_sec;
  //rseed = 1104160409;
  //rseed = 1107231717;
  //rseed = 1110226380;
  srand(rseed);
  fprintf(stderr,"seeded random number generator with %d\n",rseed);

  // write out CPTs raw format
  ///mdp->writeCPTRaw("/tmp/rawcpt.txt");
  //mdp->writeMDPDot("/tmp/handw1");

  // for tiger problem to make everything positive
  //mdp->addToReward(100);

  //fprintf(stderr,"generating policy..");
  if (!havealphas || startfromalpha) {
    if (havealphas)
      mdp->readAlphasFile(alphafile);
    if (!havebeliefs) {
      mdp->generatePolicy(nsamples, mnb, outpath);
    } else {
      mdp->generatePolicy(nsamples, mnb, outpath,belieffile);
    }
    sprintf(basepath,"%s/stats.dat",outpath);
    stats = fopen(basepath,"w");
    Cudd_PrintInfo(gbm,stats);
    fclose(stats);
    sprintf(basepath,"%s/alphas.mdp",outpath);
    mdp->printAlphasFile(basepath);
  } else {
    mdp->readAlphasFile(alphafile);
    if (!havebeliefs) {
      fprintf(stderr,"must have beliefs file if alpha file specified!\n");
      exit(0);
    } else {
      mdp->readBeliefsFile(belieffile);
    }
    // select the beliefs you want to print here
    // they get generated in printConditionalPlansForBeliefs
    // so as long as the random seed is the same, the beliefs should be too
    int npb= mdp->numbeliefs;
    npb = 17;
    int *pbel = new int[npb];
    for (i=0; i<npb; i++) {
      pbel[i] = i;
    }
    pbel[0] = 0;    pbel[1] = 1;    pbel[2] = 2;
    pbel[3] = 10;    pbel[4] = 19;    pbel[5] = 20;
    pbel[6] = 33;    pbel[7] = 35;    pbel[8] = 46;
    pbel[9] = 48;    pbel[10] = 52;    pbel[11] = 60;
    pbel[12] = 70;    pbel[13] = 89;    pbel[14] = 131;
    pbel[15] = 216;    pbel[16] = 222; 
    //if (!funky) 
    //mdp->printConditionalPlansForBeliefs(pbel,npb,outpath);
    // if simulating a continuous POMDP, use 'false' or nothing as fourth (funky) argument,
    // otherwise, if simulating a discrete one using a continuous POMDP to generate samples, use 'true'
    double *rew = new double[numtrials];
    double totrew(0.0),rewstd,meanrew;
    if (!funky)
      mdp->allocate2Dstructures(0);
    for (i=0; i<numtrials; i++) {
      rew[i] = mdp->simulator(50,outpath,verbose, funky);
      totrew += rew[i];
      if (i>0)
	fprintf(stderr,"reward on %d trial: %f  total so far: %f  average per trial: %f\n",i,rew[i],totrew,totrew/((double) (i+1)));
    }
    meanrew = totrew/numtrials;
    rewstd = 0.0;
    for (i=0; i<numtrials; i++) 
      rewstd += (rew[i]-meanrew)*(rew[i]-meanrew);
    rewstd = sqrt(rewstd/(numtrials-1.0));
    sprintf(basepath,"%s/avgsimout.txt",outpath);
    crapFILE = fopen(basepath,"a");
    fprintf(crapFILE,"reward gathered total: %f  average: %f   std.dev: %f\n",totrew,totrew/numtrials,rewstd);
    //fprintf(crapFILE,"%f %f %f\n",totrew,totrew/numtrials,rewstd);
    fclose(crapFILE);
    
  }
  /*
  int nmix = 1;
  int fvdim = 1;
  double weight;
  double *mn = new double[fvdim];
  double **cv = new double*[fvdim];
  for (i=0; i<fvdim; i++) {
    cv[i] = new double[fvdim]; 
    for (k=0; k<fvdim; k++) 
      cv[i][k] = 0.0;
  }
  MixGauss *mgg2 = new MixGauss(nmix,fvdim);
  for (i=0; i<nmix; i++) {
    weight = ((double) rand())/((double) RAND_MAX+1.0);
    for (j=0; j<fvdim; j++) {
      mn[j] =  20*(((double) rand())/((double) RAND_MAX+1.0) - 0.5);
      cv[j][j] = 0.1;
    }
    mgg2->setGaussian(i,weight,mn,cv);
  }
  fprintf(stderr,"%s\n",mgg2->toString());

  weight = 1.0;
  cv[0][0] = 0.1;

  MixGauss **ns = new MixGauss*[2];
  MixGauss **nc = new MixGauss*[2];
  ns[0] = new MixGauss(nmix,fvdim);
  ns[1] = new MixGauss(nmix,fvdim);
  nc[0] = new MixGauss(nmix,fvdim);
  nc[1] = new MixGauss(nmix,fvdim);
  mn[0] =  -10.0;
  ns[0]->setGaussian(0,weight,mn,cv);
  mn[0] =  10.0;
  ns[1]->setGaussian(0,weight,mn,cv);
  mn[0] =  -0.1;
  nc[0]->setGaussian(0,weight,mn,cv);
  mn[0] =  0.1;
  nc[1]->setGaussian(0,weight,mn,cv);

  fprintf(stderr,"erfc(1) is %lf\n",erfc(1.0));
  double intval = ns[0]->integrate(1,100,0);
  fprintf(stderr,"integral is %f\n",intval);


  MixGauss ngg1((*ns[0])+(*ns[1]));
  MixGauss ngg2((*ns[0])+(*ns[1]));
  double weights[2];
  weights[0] = 0.81;
  weights[1] = 0.11;
  ngg1.set_weights(weights);
  weights[0] = 0.09;
  weights[1] = 0.19;
  ngg2.set_weights(weights);
  MixGauss nggs[2];
  nggs[0].copyFrom(ngg1);
  nggs[1].copyFrom(ngg2);
  intval = partitionIntegrateObsSpace(nggs,2);
  fprintf(stderr,"integral is %lf\n",intval);

  double *zc = new double[1024];
  fprintf(stderr,"computing intersections ....");
  int numints = ngg1.intersections(ngg2,zc);
  fprintf(stderr,"%d found: ",numints);
  for (i=0; i<numints; i++) 
    fprintf(stderr,"%f ",zc[i]);
  fprintf(stderr,"\n");

  DdNode **obsfun = new DdNode*[2];

  DdNode *ddmg = Cudd_addConst(gbm,ns[0]);
  Cudd_Ref(ddmg);
  temp2 = mdp->buildOneCubeOrig(0,0,true);
  temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,ddmg);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,ddmg);

  temp1 = Cudd_addApply(gbm,Cudd_addMinus,One,temp2);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,temp2);

  ddmg = Cudd_addConst(gbm,ns[1]);
  Cudd_Ref(ddmg);
  temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,ddmg);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);
  obsfun[0] = Cudd_addApply(gbm,Cudd_addPlus,temp2,temp);
  Cudd_Ref(obsfun[0]);
  Cudd_RecursiveDeref(gbm,temp2);
  Cudd_RecursiveDeref(gbm,temp);
  fprintf(stderr,"obsfun[0] is \n");
  Cudd_PrintDebug(gbm,obsfun[0],4,100);

  ddmg = Cudd_addConst(gbm,nc[0]);
  Cudd_Ref(ddmg);
  temp2 = mdp->buildOneCubeOrig(0,0,true);
  temp = Cudd_addApply(gbm,Cudd_addTimes,temp2,ddmg);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,ddmg);

  temp1 = Cudd_addApply(gbm,Cudd_addMinus,One,temp2);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,temp2);

  ddmg = Cudd_addConst(gbm,nc[1]);
  Cudd_Ref(ddmg);
  temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,ddmg);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp1);
  obsfun[1] = Cudd_addApply(gbm,Cudd_addPlus,temp2,temp);
  Cudd_Ref(obsfun[1]);
  Cudd_RecursiveDeref(gbm,temp2);
  Cudd_RecursiveDeref(gbm,temp);
  fprintf(stderr,"obsfun[1] is \n");
  Cudd_PrintDebug(gbm,obsfun[1],4,100);
  
  DdNode *nep = mdp->getNewPrime(0,0);
  Cudd_Ref(nep);
  tempM = Cudd_addApply(gbm,Cudd_addTimes,nep,obsfun[0]);
  Cudd_Ref(tempM);
  fprintf(stderr,"action 0 is \n");
  Cudd_PrintDebug(gbm,tempM,4,100);
  Cudd_RecursiveDeref(gbm,nep);

  // multiply by Reward (primed)
  nep = mdp->getRewardPrimed();
  Cudd_Ref(nep);
  temp1 = Cudd_addApply(gbm,Cudd_addTimes,nep,tempM);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,tempM);
  Cudd_RecursiveDeref(gbm,nep);
  tempM = temp1;
  fprintf(stderr,"after * reward its \n");
  Cudd_PrintDebug(gbm,tempM,4,100);
  
  // sum out primed vars
  temp1 = mdp->sumOutOnePrime(tempM,0);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,tempM);
  fprintf(stderr,"after sum its \n");
  Cudd_PrintDebug(gbm,temp1,4,100);

  // now, compute the function at a belief point
  // (init belief here for example)
  tempM = mdp->getInitBeliefPrimed();
  Cudd_Ref(tempM);
  
  temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,tempM);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,tempM);
  Cudd_RecursiveDeref(gbm,temp1);
  temp1 = temp2;
  fprintf(stderr,"at initial belief its \n");
  Cudd_PrintDebug(gbm,temp1,4,100);
  
  // now, we will normally have at least two of these for each belief point
  // so now we find the integration limits by finding the all intersection points
  // of each pair. So we find the region in which each such function is maximal,
  // then integrate each belief point over its region.
  // finally, we sum up all these integrations.
  
  
  // now, integrate over all observations
  MixGauss *limitmg = new MixGauss(2,1);
  limitmg->mixweights[0] = 0;  // normally this would be the low end
  limitmg->mixweights[1] = 0;  // the high end
  limitmg->val = 3;   // the type  (0 = normal, 1 = from -inf, 2 = to +inf, 3 = -inf:+inf)
  tempM = Cudd_addConst(gbm,limitmg);
  Cudd_Ref(tempM);
  temp2 = Cudd_addApply(gbm,integrate,temp1,tempM);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,tempM);
  fprintf(stderr,"after integration its \n");
  Cudd_PrintDebug(gbm,temp2,4,100);



  // add the reward
  nep = mdp->getReward();
  Cudd_Ref(nep);
  tempM = Cudd_addApply(gbm,Cudd_addPlus,temp2,nep);
  Cudd_Ref(tempM);
  Cudd_RecursiveDeref(gbm,temp2);
  fprintf(stderr,"after adding reward its \n");
  Cudd_PrintDebug(gbm,tempM,4,100);
  

  MixGauss *mggz = new MixGauss(nmix,fvdim);
  *mggz = (*mgg2)-(*mgg2);
  tempM = Cudd_addConst(gbm,mggz);
  Cudd_PrintDebug(gbm,tempM,4,100);
 
  */
  if (outputfile) {
    strcpy(basepath,outputfile);
  } else {
    sprintf(basepath,"./SPUDD");
  }  
}
