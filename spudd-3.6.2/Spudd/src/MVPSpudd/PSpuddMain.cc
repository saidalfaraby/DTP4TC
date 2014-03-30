#include "pspudd.hh"
//#define POLICYQUERY
//#define POLICYSERVER
int yyparse();
char firstchar;
double supnormOV_AV, supnormOV_APV;
DdNode **UPList; 

int main(int argc, char *argv[])
{
  int counter,newcounter,approx_flag,flag,spuddvalue,numwrite,numread;
  int i,j,k,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int approxMeth, optimalFrom,limit_type,reorderMeth,shuffleFlag,dotFlag;
  int templist[MAXVARS];
  char templog[MAXLINE];
  char basepath[MAXLINE];
  char *argflag, *argval, *outputfile;

  double global_max_error,max_error,tmp_error;
  double temps0,temps1,temps_garbage,temps_reorder;   
  int bigadd;
  /* Allocate storage for timestamps */
  struct timeval *Time0;  
  struct timeval *Time1;
  double temps;
  long curtime;
  unsigned seed;
  Pair maxOPTP, maxPAIR;
  double maxApp,minApp;

  limit_type = LIMIT_SIZE;      /* can be LIMIT_SIZE or LIMIT_ERROR */
  reorderMeth = REORDER_MINSPAN;        /* can be REORDER_(SIFT | RANDOM |MINSPAN | NONE)*/
  long int maxusememory=0,usememory=0;
  
  DdNode *cube,*temperror,*tempo,*temp11,*offset,*temp21,*tempround;
  DdNode *ApproximationError;
  DdNode *temp1,*temp2,*one,*crap,*tempM,*maxADD;

  FILE *debug,*value,*action,*errf,*both,*stats,*log,*log2,*reward,*actions;
  FILE *crapFILE, *dataFile, *dataFile2;
  
  // intiialize the dd manager
  gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,MAXMEM);
  /* yyin is the lex global var for input file */
  /* default values first*/
  bigadd = BIGADD;
  Max_Size = MAX_SIZE;
  Max_Error = MAX_ERROR;
  limit_type = LIMIT_SIZE;
  reorderMeth = REORDER_NONE;
  approxMeth = APPROX_NONE;
  optimalFrom = OPTIMAL_NONE;
  toleranceMode = TOLERANCE_FIXED;
  shuffleFlag = 0;
  dotFlag = 1;
  
  
  outputfile = NULL;
  
  if (argc >= 3) {
    
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
      } else if (strcmp(argflag,"dd") == 0) {
	dotFlag = 0;;
      }
      argc -= 1;
    }
  } else {
    fprintf(stderr,"usage: Spudd <inpath> <infile> <outpath> [-b  bigadd ] [-a approx-method] [-r reorder-method] [-m approx-mode] [-s max_size] [-e max_error]  [-t tolerance_mode] [-f optimal-from]  [-sh] [-o <output file>] [-dd]\n");
    fprintf(stderr,"\t<inpath> path of input files (terminated with /) \n");
    fprintf(stderr,"\t<infile> input file name (will be used as a prefix for output files) \n");
    fprintf(stderr,"\t<outpath> path for output files (terminated with /) \n");
    fprintf(stderr,"\tbigadd:  bigadd constant (default 10000) \n");
    fprintf(stderr,"\tapprox-method: can be none,round or allpair (default none)\n");
    fprintf(stderr,"\treorder-method: can be none, sift, random or minspan (default none)\n");
    fprintf(stderr,"\tapprox-mode: can be error or size - only valid if -a flag specifies other than none (default size)\n");
    fprintf(stderr,"\tmax_size: value of max_size (default 10000)\n");
    fprintf(stderr,"\tmax_error: value of max_error (default 0.00001)\n");
    fprintf(stderr,"\ttolerance_mode:can be sliding or fixed (default fixed)\n");
    fprintf(stderr,"\toptimal_from: can be none, file or generate (default none) - if file, optimal value read from <outpath>\n");
    fprintf(stderr,"\t\tin the form <outpath><infile>-OptimalValue.ADD\n");
    fprintf(stderr,"\t-sh shuffle flag randomly shuffles variable order before value iteration\n");
    fprintf(stderr,"\t-o <output file> file to output dual ADD to\n");
    fprintf(stderr,"\t-dd supresses dumping of output .dot files\n");
    

    
    exit(0);
  }
  
  setbuf(stdout,0);
  /* Structure to get the time of execution*/
  struct tms _t;


  pZero = new Pair(0.0);
  pOne = new Pair(1.0);
  
  One = Cudd_addConst(gbm,pOne);
  Cudd_Ref(One);
  Zero = Cudd_addConst(gbm,pZero);
  Cudd_Ref(Zero);
  
  /* We initialize the DdManager  */
  //gbm.makeVerbose();
  //gbm.makeTerse();
  //gbm.GarbageCollectionEnabled();
  //Cudd_DisableGarbageCollection(gbm);

  //Disables automatic dynamic reordering
  //gbm.AutodynDisable;
  Cudd_AutodynDisable(gbm);

  strcpy(basepath,inpath);
  strcat(basepath,infile);
  yyin = fopen(basepath,"r");

  setbuf(stdout,0); 
  
  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-error.dat");
  strcat(basepath,templog);
  errf = fopen(basepath,"w");

  // check for boolean input file
  mvinput=true;
  fscanf(yyin,"%c",&firstchar);
  if (firstchar == 'v') 
    mvinput = false; 

  // reset the input file
  rewind(yyin);


  /* Parse the input file */  
  numvars =0;
  numactions=0;
  numorigvars = 0;
  //yydebug = 1;
  if (mvinput) {
    yyparse();
  } else {
    if (!pparse(errf,&discount_factor,&horizon, 
	       &tolerance,&numvars,&numactions,vars,prime_vars,actionlist)) { 
      fprintf(errf,"parse error\n");
      Cudd_Quit(gbm);
      exit(0);
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


  /* convert the pruning percentage to error by using the span and discount */
  extR = get_extent(RewardD);
  //double discount_factor = ((*Cudd_V(discount)).get_min()) ;

  switch (toleranceMode) 
    {
    case TOLERANCE_FIXED:
      Max_Error = prune*discount_factor*extR/(1-discount_factor);
      break;
    case TOLERANCE_SLIDING:
      Max_Error = prune*extR;
      break;
    }

  
  //fprintf(stderr,"extent: %f discount: %f Max_Error: %f\n",extR,discount_factor,Max_Error);
  //Reward = ADD(&gbm,RewardD);
  
  /* CREATING THE CONSTANTS FOR THE DISCOUNTING AND THE PRIMING OF THE ACTION */
  toleranceP = new Pair(tolerance);
  //(*toleranceP).set(tolerance);
  
  
  Pair half;
  half.set(0.5);
  Half = Cudd_addConst(gbm,&half);
  Cudd_Ref(Half);

  Pair verysmall;
  verysmall.set(-500000000.0);
  VerySmall = Cudd_addConst(gbm,&verysmall);
  Cudd_Ref(VerySmall);
  
  Pair discountPair;
  discountPair.set(discount_factor);
  discount = Cudd_addConst(gbm,&discountPair);
  Cudd_Ref(discount);

  /*for (i=0; i< 2*numvars; i++)
    fprintf(stderr,"%d ",i);
  fprintf(stderr,"\n");
  for (i=0; i<2*numvars;i++)
    fprintf(stderr,"%d ",Cudd_ReadInvPerm(gbm,i));
  fprintf(stderr,"\n");

  */


  /* random seed */
  curtime = time(NULL);
  seed = (unsigned) curtime % INT_MAX;
  srand(seed);

  //srand();
  //fprintf(stderr,"seed: %d\n",seed);

  //Cudd_EnableReorderingReporting(gbm);


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
    
  }
  
  
  // Now, build up the Allprime trees - We want to build a minimum number of allprime
  // trees for each action, such that each contains less than bigadd nodes
  // We multiplys all the NewPrimed ADD until they reach a upper limit size...then we build another ADD.
  // We repeat that for each action.
  
  
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
	  fprintf(errf,"Memory required exceeds availabe %ld bytes. \n Try a smaller bigadd limit constant\n",MAXMEMHARD);
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
	  //  fprintf(stats,"The %d allprime diagram got too big - %d nodes\n",k,dagCount);
	  flag = 0;
	  /** NEW FEB 4 **/
	  Cudd_RecursiveDeref(gbm,temp2);
	}
      }
      
      //fprintf(stderr,"Assigned complete action diagram %d for action %d with %d nodes \n",k-1,i,dagCount);
      
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

  /* write out the input data as dot files */
  if (DUMPIN) {
    strcpy(basepath,outpath);
    strcat(basepath,"reward.dot");
    log = fopen(basepath,"w");
    Array[0]=RewardD;
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars+2],NULL,log);
    fclose(log);
    
    strcpy(basepath,outpath);
    strcat(basepath,"reward.html");
    log = fopen(basepath,"w");
    dumpHtmlRew(log,outpath);
    fclose(log);
    
    if (!mvinput) {
      for (i=0; i<numactions; i++) {
	strcpy(basepath,outpath);
	sprintf(templog,"action%d.html",i);
	strcat(basepath,templog);
	log2 = fopen(basepath,"w");
	fprintf(log2,"<html>\n\t<body><h1>Action: ");
	fprintf(log2,actionnames[i]);
	fprintf(log2,"</h1>\n<table rules=n>\n<tr>");
	for (j=0; j<numvars; j++) {
	  strcpy(basepath,outpath);
	  sprintf(templog,"action%d-%d.dot",i,j);
	  strcat(basepath,templog);
	  log = fopen(basepath,"w");
	  Array[0] = Actions[i][j];
	  My_DumpDot(gbm,1,Array,varnames,&lnames[j],NULL,log);
	  fclose(log);
	  fprintf(log2,"<td>\n<img src=http://www.research.att.com/~north/cgi-bin/webdot.cgi/http://www.cs.ubc.ca/spider/jhoey/spudd/");
	  fprintf(log2,outpath);
	  sprintf(templog,"action%d-%d.dot.gif>\n",i,j);
	  fprintf(log2,templog);
	  fprintf(log2,"</td>\n");
	}
	fprintf(log2,"</table>\n\n\t<hr>\n</body>\n</html>\n");
	fclose(log2);
      }
    }
      // re-order so all primed variables are last
      int *nlist = new int[2*numvars];
      int *orig_list = new int[2*numvars];
      for (i=0; i<numvars*2; i++) 
	// this is the original ordering
	orig_list[i] = Cudd_ReadInvPerm(gbm,i);
      for (i=0; i<numvars; i++) 
	nlist[i] = 2*i+1;
      for (i=numvars; i<2*numvars; i++) 
	nlist[i] = 2*(i-numvars);;
      
      
      int res = Cudd_ShuffleHeap(gbm,nlist);

      for (i=0; i<numactions; i++) {
	for (j=0; j<numorigvars; j++) {
	  strcpy(basepath,outpath);
	  sprintf(templog,"newprime%d-%d.dot",i,j);
	  strcat(basepath,templog);
	  log = fopen(basepath,"w");
	  Array[0] = NewPrime[i][j];
	  My_DumpDot(gbm,1,Array,varnames,NULL,NULL,log);
	  fclose(log);
	}
      }
      res = Cudd_ShuffleHeap(gbm,orig_list);
      delete [] nlist;
      delete [] orig_list;
   
  }
  int varid[numvars];
  for(i=0;i<numvars;i++)
    varid[i] = vars[i].number;
  
  char *Vnames[numvars];
  for(i=0;i<numvars;i++)
    Vnames[i] = varnames[i];



  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-dataO.dat"); 
  strcat(basepath,templog);
  dataFile = fopen(basepath,"a");

  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-dataC.dat"); 
  strcat(basepath,templog);
  dataFile2 = fopen(basepath,"a");
  

  fprintf(dataFile,"%d ",Max_Size);
  fprintf(dataFile,"%f ",Max_Error);
  fprintf(dataFile,"%d ",approxMeth);
  fprintf(dataFile,"%d ",reorderMeth);
  
  fprintf(dataFile2,"%d ",seed);
  fprintf(dataFile2,"%d ",Max_Size);
  fprintf(dataFile2,"%f ",Max_Error);
  fprintf(dataFile2,"%d ",approxMeth);
  fprintf(dataFile2,"%d ",reorderMeth);
  
  
    
  if (shuffleFlag) {
    fprintf(stderr,"Shuffling order...");
    shuffleRandom(); 
    fprintf(stderr,"done\n");
  }
  
    
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
      if (!reorderMinSpan(RewardD))
	fprintf(stderr,"******ERROR: Reorder minSpan failed \n");
      break;
    case REORDER_EXACT:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_EXACT, 0))
	fprintf(stderr,"*****ERROR: Reorder exact failed\n");
      break;
    case REORDER_NONE:
      break;
    }


  
  
  switch (approxMeth) 
    {
    case APPROX_NONE:      
      for (i=0; i<2*numvars;i++) 
	fprintf(dataFile,"%d ",Cudd_ReadInvPerm(gbm,i));
      fprintf(dataFile,"\n");
      OptimalValue  = spudd_Proc(reorderMeth,dataFile2);
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      //fprintf(stderr,"Generating Approximate policy...");
      // also sets ApproximateMidPointValue to be the mid-point aDD
      // and ApproximateValue to be the approximate value function (ranged)
      ApproximateValue = NewAspudd_Proc(limit_type,reorderMeth,approxMeth,dataFile, dataFile2);
      
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
      
      
      switch (optimalFrom) 
	{
	case OPTIMAL_FILE:
	  for(i=0;i<numvars;i++)
	    varid[i] = Cudd_ReadInvPerm(gbm,i);
	  OptimalValue = 
	    Dddmp_cuddAddLoad(gbm,DDDMP_VAR_MATCHIDS,NULL,varid,NULL,DDDMP_MODE_TEXT,basepath,crapFILE);
	  break;
	case OPTIMAL_GENERATE:
	  //fprintf(stderr,"Generating optimal policies...");
	  OptimalValue  = spudd_Proc(reorderMeth, dataFile2);
	  //fprintf(stderr,"done\n");
	  
	  break;
	}
      if (optimalFrom != OPTIMAL_NONE) {
	//fprintf(stderr,"Evaluating it...");
	ApproximatePolicyValue = evalPolicy(RewardD, ApproximatePolicy);
	//fprintf(stderr,"done\n");
       	// approximation error between unnormalized optimalvalue and approximatepolicy value
	ApproximationError = Cudd_addApply(gbm,Cudd_addMinus,OptimalValue,ApproximatePolicyValue);
	Cudd_Ref(ApproximationError);
	
	supnormOV_APV = (*Cudd_V(Cudd_addFindMax(gbm,ApproximationError ))).get_max();
	SUMSQUARE = 0.0;
	numSumSquare = 0;
	(void) My_addApply(SumSquareError,ApproximationError,ApproximationError);
	double policyValuess = SUMSQUARE;
	int policyValuenss = numSumSquare;

	// now, normalize the optimal value function
	// and get the error with the approximate mid point value function
	double maxOPT = (*Cudd_V(Cudd_addFindMax(gbm,OptimalValue))).get_max();
	double minOPT = (*Cudd_V(Cudd_addFindMin(gbm,OptimalValue))).get_min();
	
	maxOPTP.set(1/(maxOPT-minOPT));
	maxADD = Cudd_addConst(gbm,&maxOPTP);  
	Cudd_Ref(maxADD);
	  
	tempM = Cudd_addApply(gbm,Cudd_addTimes,OptimalValue,maxADD);
	Cudd_Ref(tempM);
	Cudd_RecursiveDeref(gbm,OptimalValue);
	Cudd_RecursiveDeref(gbm,maxADD);
	OptimalValue = tempM;
	Cudd_Ref(OptimalValue);
	Cudd_RecursiveDeref(gbm,tempM);
	
	ApproximationError = Cudd_addApply(gbm,Cudd_addMinus,OptimalValue,ApproximateMidPointValue);
	Cudd_Ref(ApproximationError);
	
	SUMSQUARE = 0.0;
	numSumSquare = 0;
	(void) My_addApply(SumSquareError,ApproximationError,ApproximationError);
	
	double aspan = get_error(ApproximateValue);
	fprintf(dataFile2,"Optimal extent: %f \n",(maxOPT-minOPT)); 
	fprintf(dataFile2,"Approximate mid-point extent: %f \n",(maxApp-minApp)); 
	fprintf(dataFile2,"approximate value span: %f \n",aspan); 
	fprintf(dataFile2,"supnorm OV-APV: %f\n",supnormOV_APV);
	fprintf(dataFile2,"policy value to optimal ssqe:  %f ",sqrt(policyValuess/((float) (policyValuenss-1))));
	fprintf(dataFile2,"mid-point value to optimal ssqe:  %f ",sqrt(SUMSQUARE/((float) (numSumSquare-1))));
	Cudd_RecursiveDeref(gbm,ApproximationError);
      }
      break;
    default:
      break;
     
    }
  //fprintf(stderr,"done!!!!!!!!!\n");
  fprintf(dataFile,"\n");
      fprintf(dataFile2,"\n");
  fflush(dataFile);
  fflush(dataFile2);
  fclose(dataFile);
  fclose(dataFile2);

  // write out dot files

  //UPList= (DdNode **)malloc(2*(sizeof(DdNode *)));
  
  //Cudd_PrintDebug(gbm,OptimalPolicy,4,100);
  

  switch (approxMeth) 
    {
    case APPROX_NONE: 
      if (dotFlag) {
	sprintf(basepath,"%s%s-OPT",outpath,infile);
	writeVPDot(OptimalValue,OptimalPolicy,basepath);
      }
      //save the value add and the policy to a file. 
      if (!outputfile) 
	sprintf(basepath,"%s%s-OPTDual.ADD",outpath,infile);
      else 
	sprintf(basepath,"%s",outputfile);
      numwrite = writeDualOptimal(basepath,OptimalPolicy,OptimalValue,vars,prime_vars,orig_vars);

      
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      if (dotFlag) {
	sprintf(basepath,"%s%s-APP",outpath,infile);
	writeVPDot(ApproximateValue,ApproximatePolicy,basepath);
      }

      //save the value add and the policy to a file.
      if (!outputfile) 
	sprintf(basepath,"%s%s-APPDual.ADD",outpath,infile);
      else 
	sprintf(basepath,"%s",outputfile);
      numwrite = writeDualOptimal(basepath,ApproximatePolicy,ApproximateValue,vars,prime_vars,orig_vars);


      break;
    }      
 
  
  // start policy querying app
#ifdef POLICYQUERY
  switch (approxMeth) 
    {
    case APPROX_NONE:      
      pQuery(gbm,OptimalValue,OptimalPolicy,orig_vars,vars,numvars,numorigvars);
      break;
    case APPROX_ROUNDOFF:
    case APPROX_ALLPAIRS:
      pQuery(gbm,ApproximateValue,ApproximatePolicy,orig_vars,vars,numvars,numorigvars);
      break;
    }    
#endif

#ifdef POLICYSERVER
  if (mvinput) {
    switch (approxMeth) 
      {
      case APPROX_NONE:      
	policyServe(gbm, OptimalPolicy, OptimalValue, numorigvars);
	break;
      case APPROX_ROUNDOFF:
      case APPROX_ALLPAIRS:
	policyServe(gbm, ApproximatePolicy, ApproximateValue, numorigvars);

	break;
      }    
  }


#endif




  /*********************************************************************************
  From here on down its just cleaning up and writing out some files
  *********************************************************************************/

 

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
	    Cudd_RecursiveDeref(gbm,ActionsPrimed[i][j]);  
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
      if (optimalFrom != OPTIMAL_NONE && optimalFrom != OPTIMAL_GENERATE) {
	Cudd_RecursiveDeref(gbm,ApproximationError);
	Cudd_RecursiveDeref(gbm,OptimalValue);
      }
    }
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
  
  if(ActionsPrimed){
    for(i=0;i<numactions;i++)
      free(ActionsPrimed[i]);
    free(ActionsPrimed);
  }
  
  if(NewPrime){
    for(i=0;i<numactions;i++)
      free(NewPrime[i]);
    free(NewPrime);
  }
  Cudd_Quit(gbm);
  return 1;
}
