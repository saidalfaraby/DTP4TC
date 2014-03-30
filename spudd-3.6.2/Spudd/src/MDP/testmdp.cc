#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
//#define USE_POMDP
#ifndef USE_POMDP
#include "MDP.h"
#else
#include "POMDP.h"
#endif

//#define HANDWASHING
//#define OLDINTERFACE
int main(int argc, char *argv[])
{
  int counter,newcounter,approx_flag,flag,spuddvalue,numwrite,numread;
  int i,j,k,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int approxMeth, optimalFrom,limit_type,reorderMeth,reorderApplyMeth,shuffleFlag,dotFlag,dompi;
  int templist[MAXVARS];
  char templog[MAXLINE];
  char basepath[MAXLINE];
  char *argflag, *argval, *outputfile, *dualfile;
  char *infile, *outpath, *inpath;
  double prune;
  double global_max_error,max_error,tmp_error;
  int bigadd;
  /* Allocate storage for timestamps */
  long curtime;
  unsigned seed;
  Pair maxOPTP, maxPAIR;
  double maxApp,minApp;

  limit_type = LIMIT_SIZE;   
  long int maxusememory=0,usememory=0;
  
  DdNode *cube,*temperror,*tempo,*temp11,*offset,*temp21,*tempround;
  DdNode *ApproximationError;
  DdNode *temp1,*temp2,*one,*crap,*tempM,*maxADD;

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


  outputfile = NULL;
  
  if (argc >= 2) {
#ifdef OLDINTERFACE
    inpath = *++argv;
    infile = *++argv;
    outpath = *++argv;
    argc = argc-3;
#else
    infile = *++argv;
    argc = argc-1;
#endif
    while (--argc > 0 && (*++argv)[0] == '-') {
      argflag = ++argv[0];
      if (strcmp(argflag, "b") == 0)
	bigadd = atoi(*++argv);
      if (strcmp(argflag, "e") == 0) 
	prune = atof(*++argv);
      //Max_Error = atof(*++argv);
      else if (strcmp(argflag,"s") == 0)
	Max_Size = atoi(*++argv);
      else if (strcmp(argflag,"a") == 0) {
	argval = *++argv;
	if (strcmp(argval,"none") == 0)
	  approxMeth = APPROX_NONE;
	else if (strcmp(argval,"round") == 0)
	  approxMeth = APPROX_ROUNDOFF;
	else if (strcmp(argval,"allpair") == 0)
	  approxMeth = APPROX_ALLPAIRS;
      } else if (strcmp(argflag,"r") == 0) {
	argval = *++argv;
	if (strcmp(argval,"none") == 0)
	  reorderMeth = REORDER_NONE;
	else if (strcmp(argval,"sift") == 0)
	  reorderMeth = REORDER_SIFT;
	else if (strcmp(argval,"random") == 0)
	  reorderMeth = REORDER_RANDOM;
	else if (strcmp(argval,"minspan") == 0)
	  reorderMeth = REORDER_MINSPAN;
	else if (strcmp(argval,"exact") == 0)
	  reorderMeth = REORDER_EXACT;
	else if (strcmp(argval,"spsift") == 0)
	  reorderMeth = REORDER_SP_SIFT;
      } else if (strcmp(argflag,"ra") == 0) {
	argval = *++argv;
	if (strcmp(argval,"all") == 0)
	  reorderApplyMeth = REORDERAPPLY_ALL;
	else if (strcmp(argval,"ff") == 0)
	  reorderApplyMeth = REORDERAPPLY_FIRSTFIVE;
      } else if (strcmp(argflag,"m") == 0) {
	argval = *++argv;
	if (strcmp(argval,"error") == 0)
	  limit_type = LIMIT_ERROR;
	else if (strcmp(argval,"size") == 0) 
	  limit_type = LIMIT_SIZE;
      } else if(strcmp(argflag,"t") ==0){
	argval = *++argv;
	if(strcmp(argval,"sliding") == 0)
	  toleranceMode = TOLERANCE_SLIDING;
	else if(strcmp(argval,"fixed") == 0)
	  toleranceMode = TOLERANCE_FIXED;
      } else if (strcmp(argflag,"f") == 0) {
	argval = *++argv;
	if (strcmp(argval,"none") == 0)
	  optimalFrom = OPTIMAL_NONE;
	else if (strcmp(argval,"file") == 0)
	  optimalFrom = OPTIMAL_FILE;
	else if (strcmp(argval,"generate") == 0)
	  optimalFrom = OPTIMAL_GENERATE;
      } else if (strcmp(argflag,"sh") == 0) {
	shuffleFlag = 1;
      } else if (strcmp(argflag,"o") == 0) {
	outputfile = *++argv;
      } else if (strcmp(argflag,"dual") == 0) {
	dualfile = *++argv;
	dompi = 1;
      } else if (strcmp(argflag,"dd") == 0) {
	dotFlag = 1;
      } else if (strcmp(argflag,"tt") == 0) {
	mtol = atof(*++argv);
      }
      argc -= 1;
    }
  } else {
    fprintf(stderr,"usage: testmdp <infile> [-b  bigadd ] [-a approx-method] [-r reorder-method] [-ra reorder-apply-method] [-m approx-mode] [-s max_size] [-e max_error]  [-t tolerance_mode] [-f optimal-from]  [-sh] [-o <output file>] [-dual <dualfile>] [-dd]\n");
    fprintf(stderr,"\t<infile> input file name (will be used as a prefix for output files if none specified) \n");
    fprintf(stderr,"\tbigadd:  bigadd constant (default 10000) \n");
    fprintf(stderr,"\tapprox-method: can be none,round or allpair (default none)\n");
    fprintf(stderr,"\treorder-method: can be none, sift, random or minspan (default none)\n");
    fprintf(stderr,"\treorder-apply-method: can be all or ff (ff = reorder on first five iterations only  -- default all)\n");
    fprintf(stderr,"\tapprox-mode: can be error or size - only valid if -a flag specifies other than none (default size)\n");
    fprintf(stderr,"\tmax_size: value of max_size (default 10000)\n");
    fprintf(stderr,"\tmax_error: value of max_error (default 0.00001)\n");
    fprintf(stderr,"\ttolerance_mode:can be sliding or fixed (default fixed)\n");
    fprintf(stderr,"\toptimal_from: can be none, file or generate (default none) - if file, optimal value read from <outpath>\n");
    fprintf(stderr,"\t\tin the form <infile>-OptimalValue.ADD\n");
    fprintf(stderr,"\t-sh shuffle flag randomly shuffles variable order before value iteration\n");
    fprintf(stderr,"\t-o <output file> file to output dual ADD to\n");
    fprintf(stderr,"\t-dual <dual file> read the optimal value function from dualfile and start value iteration from there\n");
    fprintf(stderr,"\t-dd print output .dot files *instead of* dual ADD file\n");
       
    exit(0);
  }
#ifdef OLDINTERFACE
  dotFlag = 1;
  sprintf(basepath,"%s/%s",inpath,infile);
#else
  sprintf(basepath,"%s",infile);
#endif
  // also reads in the MDP
#ifndef USE_POMDP
  MDP *mdp = new MDP(basepath,bigadd);
#else
  POMDP *mdp = new POMDP(basepath,bigadd);
#endif
  /*
  double **T;
  int junk1, junk2;
  mdp->flattenT(&T,junk1,junk2);
  */
  fprintf(stderr,"generating policy..");
  // re-randomize after POMDP is generated
  //srand(32849239); 

  // or from the timer
  struct timeval sd;
  int rseed;
  gettimeofday(&sd,NULL);
  rseed = (int) sd.tv_sec;
  //rseed = 1094154109;
  //rseed = 1094397105;
  rseed = 1103995212;

  srand(rseed);
  fprintf(stderr,"seeded random number generator with %d\n",rseed);
#ifdef USE_POMDP
  //mdp->simulate(tolerance,100);
  mdp->generatePolicy(true);
#endif
  //mdp->pickPolicyAction();
  //mdp->modifiedPolicyIteration(false);
#ifndef USE_POMDP
  if (dompi) {
    mdp->readDualOptimal(dualfile);
    //mdp->modifiedPolicyIteration(true);
    mdp->generatePolicy(reorderMeth, reorderApplyMeth, tolerance, bigadd,  prune, Max_Size, approxMeth, limit_type, toleranceMode, shuffleFlag);
    //mdp->doshit(8);
  } else {
    mdp->generatePolicy(reorderMeth, reorderApplyMeth, tolerance, bigadd,  prune, Max_Size, approxMeth, limit_type, toleranceMode, shuffleFlag);
  }
#endif
#ifdef OLDINTERFACE
  strcpy(basepath,outpath);
  sprintf(basepath,"%s/%s",outpath,infile);
#else
  if (outputfile) {
    strcpy(basepath,outputfile);
  } else {
    sprintf(basepath,"./SPUDD");
  }  
#endif
#ifdef HANDWASHING
  mdp->simulateHandwashing("/project/coach/results/pob.dat");
#endif
#ifndef USE_POMDP
  //mdp->writeMDPDot(basepath);
  mdp->writePolicyToFile(basepath,dotFlag);
  if (approxMeth == APPROX_NONE) {
    strcat(basepath,"-stats.dat");
    mdp->printOptimalStats(basepath);
  } else {
    strcat(basepath,"-Approx-stats.dat");
    mdp->printApproximateStats(basepath);
  }
  //mdp->startServer();
  /*
  MDP *learnmdp = new MDP(mdp);
  learnmdp->reinforcementLearn(mdp,5000);

  learnmdp->writePolicyToFile(basepath,dotFlag);
  */
  /*

  learnmdp->writeMDPDot("learn");
  mdp->writeMDPDot("orig");
  fprintf(stderr,"\noriginal NewPrime\n");
  mdp->printNewPrime(stderr);
  fprintf(stderr,"\nlearned NewPrime\n");
  learnmdp->printNewPrime(stderr);
  int *mergeStates = new int[2];
  mdp->compareSubPol(0);
  */
  /*
  int nm = mdp->compareSubPol(7,mergeStates);
  fprintf(stderr,"number of merges : %d\n",nm);
  int *varvals = new int[2];
  varvals[0] = 3;
  varvals[1] = 0;
  int acttaken = mdp->consultOptimalPolicy(varvals);
  if (mtol > 0) {
    int *compmat;
    int nov;
    mdp->compareValueVals(0,mtol,&compmat,&nov);

    for (int i=0; i<nov; i++) {
      for (int j=0; j<nov; j++) 
	fprintf(stderr,"%d ",compmat[i*nov+j]);
      fprintf(stderr,"\n");
    }
    delete [] compmat;
    
  }
  */
#endif
  //delete mdp;
}
