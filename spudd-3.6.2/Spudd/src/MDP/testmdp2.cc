 #include "MDP.h"

int main(int argc, char *argv[])
{
  int counter,newcounter,approx_flag,flag,spuddvalue,numwrite,numread;
  int i,j,k,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int approxMeth, optimalFrom,limit_type,reorderMeth,reorderApplyMeth,shuffleFlag,dotFlag;
  int templist[MAXVARS];
  char templog[MAXLINE];
  char basepath[MAXLINE];
  char *argflag, *argval, *outputfile;
  char *infile, *outpath, *inpath;
  double prune,tolerance;
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
  
  
  outputfile = NULL;
  
  if (argc >= 2) {

    inpath = *++argv;
    infile = *++argv;
    outpath = *++argv;
    argc = argc-3;
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
      } else if (strcmp(argflag,"dd") == 0) {
	dotFlag = 1;
      }
      argc -= 1;
    }
  } else {
    fprintf(stderr,"usage: testmdp <infile> [-b  bigadd ] [-a approx-method] [-r reorder-method] [-ra reorder-apply-method] [-m approx-mode] [-s max_size] [-e max_error]  [-t tolerance_mode] [-f optimal-from]  [-sh] [-o <output file>] [-dd]\n");
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
    fprintf(stderr,"\t-dd print output .dot files *instead of* dual ADD file\n");
       
    exit(0);
  }

  dotFlag = 1;

  // also reads in the MDP
  MDP *mdp = new MDP(infile,bigadd);

  mdp->generatePolicy(reorderMeth, reorderApplyMeth, tolerance, bigadd,  prune, Max_Size, approxMeth, limit_type, toleranceMode, shuffleFlag);

  if (outputfile) {
    strcpy(basepath,outputfile);
  } else {
    sprintf(basepath,"./SPUDD");
  }    
  mdp->writePolicyToFile(basepath,dotFlag);

  if (approxMeth == APPROX_NONE) {
    strcat(basepath,"-stats.dat");
    mdp->printOptimalStats(basepath);
  } else {
    strcat(basepath,"-Approx-stats.dat");
    mdp->printApproximateStats(basepath);
  }
  //mdp->startServer();
  
}
