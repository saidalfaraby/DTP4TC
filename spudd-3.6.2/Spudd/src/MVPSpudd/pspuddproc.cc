/***********************************************************************
  SPUDDPROC: Optimal value iteration with ADDS
  $Author: jhoey $
  $Date: 2003/01/16 22:44:00 $
  $Revision: 2.0 $
  $Log: pspuddproc.cc,v $
  Revision 2.0  2003/01/16 22:44:00  jhoey
  got everything working now with reading/writing, mutli-valued vars
  it also works with the original data style

  Revision 1.9  2003/01/13 19:47:29  jhoey
  getting policyQuery and policyServer to work

  Revision 1.8  2003/01/11 00:41:52  jhoey
  modified the policy generation thing so it doesn't keep all the Q
  functions every iteration of vi.
  also got adding a new policy to the policyQuery gui to work
  still have to make the makefile

  Revision 1.7  2003/01/10 22:16:36  staubin
  I created the first version of the load and store ADD functions....still some slight changes to be made.

  Revision 1.6  2003/01/10 03:10:16  jhoey
  working on policy query stuff

  Revision 1.5  2002/12/13 23:11:00  jhoey
  fixed round off method

  Revision 1.4  2002/12/12 21:55:33  jhoey
  *** empty log message ***

  Revision 1.3  2002/12/12 21:48:21  jhoey
  fixed dummy variable problem

  Revision 1.2  2002/12/10 17:54:17  jhoey
  fixed a few bugs

  Revision 1.1  2002/12/05 06:43:08  jhoey
  added mutli-valued variables to PSpudd

  Revision 1.4  2002/09/25 20:21:22  jhoey
  CLK_TCK fixed to comply with new standard

  Revision 1.3  2002/09/21 00:55:15  jhoey
  added action costs

  Revision 1.2  2001/04/12 21:10:08  jhoey
  corrected printing of constant nodes to use toString

  Revision 1.1  2001/04/06 23:30:08  staubin
  testing....

  Revision 1.3  2000/10/05 21:26:10  jhoey
  added cvs authorship stuff at header


   Code by: Jesse Hoey
            Robert St-Aubin
***********************************************************************/
#include "pspudd.hh"
#include <string.h>
#include <math.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>



DdNode *spudd_Proc(int reorderMeth, FILE *dataFile2)
{
  int test,counter;
  int i,j,k,tempo,tempn,tempfoo,index,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int templist[MAXVARS];
  int flag;
  long int maxusememory=0,usememory=0;
  DdNode *actionDD;
  DdNode *temp1,*temp2,*tempround;
  DdNode **Vcurrent,*Vpast,*Merged_Add;
  DdNode *Reward;
  char templog[MAXLINE];
  char basepath[MAXLINE];
  double temps0,temps1,temps;   
  struct timeval *Time0;       /* Timestamps */
  struct timeval *Time1;
  
  FILE *value,*action,*errf,*both,*stats,*log,*log2,*reward,*actions;
   
  /* Allocate storage for timestamps */
  //Time0 = (struct timeval*)malloc(sizeof(struct timeval));
  //Time1 = (struct timeval*)malloc(sizeof(struct timeval));
  
  /* Structure to get the time of execution*/
  struct tms _t;
  long clk_tck = sysconf(_SC_CLK_TCK);
  

    
  Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));

  for(i=0;i<numactions;i++) {
    Vcurrent[i] = RewardD;
    Cudd_Ref(Vcurrent[i]);
  }  
  test = 0;
  counter = 0;
  
  //We initialize our Add structures

   /* Get the starting up time */
  //gettimeofday(Time0, NULL );
  
  /* We get the starting up time */
  times(&_t);
  temps0 = 1.0*_t.tms_utime/clk_tck;
  
   
  Vpast = RewardD;
  Cudd_Ref(Vpast);
  
  Merged_Add = RewardD;
  Cudd_Ref(Merged_Add);
    
  Reward = RewardD;
  Cudd_Ref(Reward);


  //Opening the stats files. Amking sure that we don't exceed the maximum of memory available
  strcpy(basepath,outpath);
  strcat (basepath,infile);
  sprintf(templog,"-stats.dat"); 
  strcat(basepath,templog);
  stats = fopen(basepath,"w");
  fprintf(stats,"\n\n\n Name of input file: %s\n\n Discount factor is : %f \n\n\nTolerance is: %f \n horizon is: %f \n The BIGADD limit is set to: %d \n\n\n Target maximum memory use: %d bytes \n\n Hard Limit on memory use %d bytes\n\n\n",infile, discount_factor, tolerance, horizon, bigadd,MAXMEM, MAXMEMHARD);
  
  fprintf(stats,"\n\n");
  
  actionDD = Zero;	/* Stupid ref to maintain loop invariant */
  Cudd_Ref(actionDD);
  
  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  fprintf(stats,"Memory in use before start of iterations: %ld\n",usememory);
  fprintf(stats,"\n\n\n");
  fclose(stats);
  
  int lastiteration=0;

  /* VALUE ITERATION LOOP COMMENCES */
  while ((horizon == -1 && (test==0 || lastiteration == 1)) || (horizon >= 0 && counter < horizon)){ 
    

    /* KEEP TRACK OF THE NUMBER OF ITERATIONS */
    counter = counter + 1;
    //fprintf(stderr,"\niteration %d\n",counter);
    
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

    //dagCount = Cudd_DagSize(Merged_Add);
    //fprintf(stderr,"value function before multiply size is  %d\n",dagCount);

    /* Merged_Add maintains the current Value function*/
    Cudd_RecursiveDeref(gbm,Merged_Add);
    Merged_Add = VerySmall;
    Cudd_Ref(Merged_Add);
    
 
    /*Loop over the actions */
    for(i=0;i<numactions;i++) { 

      /* do the sum_t(Pr(s,a,t)*V(t)) */
#ifdef OPTIMIZED
      temp1 = primeReward(Reward,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
      temp1 = newPrimeReward(Reward,NewPrime[i]);
#endif

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

      // temp2 is the Q function for this action
      // old way to get policies was to save these Q functions
      // if there are lots of actions, this causes problems.
      /*
      Cudd_RecursiveDeref(gbm,Vcurrent[i]);
      Vcurrent[i] = temp2;
      Cudd_Ref(Vcurrent[i]);
      */

      /* maximize over actions */
      temp1=Cudd_addApply(gbm,Cudd_addMaximum,Merged_Add,temp2);
      Cudd_Ref(temp1);
      
      // new way - update the policy at each iteration
      // only need this on the final iteration
      if (lastiteration || (horizon >= 0 && counter == horizon)) {
	getActionAdds(&temp2,temp1,1,&actionDD,i+1);
      }

      Cudd_RecursiveDeref(gbm,temp2);
      Cudd_RecursiveDeref(gbm,Merged_Add);
      Merged_Add = temp1;
      
	  
      usememory = Cudd_ReadMemoryInUse(gbm);
      if (maxusememory < usememory)
	maxusememory = usememory;
      if (usememory > MAXMEMHARD) {
	fprintf(errf,"Memory required exceeds availabe %ld bytes. \nTry a smaller bigadd limit constant\n",MAXMEMHARD);
	Cudd_Quit(gbm);
	exit(0);
      }
    } // end of loop over actions
    /* COMPARE THE TWO MOST RECENT ADDs IN ORDER TO VERIFY THE CONVERGENCE  */
    test = Cudd_EqualSupNorm(gbm,Merged_Add,Vpast,toleranceP,0);

    if (test>0 || lastiteration)
      lastiteration++;
    
    /* attempt a reorder? */
    /* sometimes, we just want to reorder for a few of the first (5?)
       iterations - this actually gives the best performance, but
       is not implemented as  a flag yet 
    */
    //if (counter < 5) {
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
	if (!reorderMinSpan(Merged_Add))
	  fprintf(stderr,"******ERROR: Reorder minSpan failed \n");
	break;
      case REORDER_EXACT:
	if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_EXACT, 0))
	  fprintf(stderr,"*****ERROR: Reorder exact failed\n");
	break;
      case REORDER_NONE:
	break;
      default:
	break;
      }
    //}      // for the counter<5 loop - not implemented as a flag yet
  }   /* END OF VALUE ITERATION */

  /* We get the time of completion*/
  /*
    gettimeofday( Time1, NULL );
    temps = (double) ((Time1->tv_sec - Time0->tv_sec) * 1000000 + (Time1->tv_usec -Time0->tv_usec));
    temps = temps/1000000.0;
  */
  
  /* We get the time of completion*/
  times(&_t);
  temps1 = 1.0*_t.tms_utime/clk_tck;
  
  temps = temps1- temps0;
  /*
  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-OPToval.dot");
  strcat(basepath,templog);
  
  value = fopen(basepath,"w");
  
  Array[0] = Merged_Add;
  Cudd_Ref(Array[0]);
  My_DumpDot(gbm,1,Array,varnames,&lnames[numvars+1],NULL,value);  
  Cudd_RecursiveDeref(gbm,Array[0]);
  
  fclose(value);
  */


  // this has already been calculated now
  //actionDD = getActionAdd(Vcurrent, Merged_Add, One, Zero, numvars, numactions);
  // this code dumps out the binary action diagram
  /*
  
  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-OPTaction.dot");
  strcat(basepath,templog);
  
  action = fopen(basepath,"w");
  
  Array[0] = actionDD;
  Cudd_Ref(Array[0]);
  My_DumpDot(gbm,1,Array,varnames,&lnames[numvars+1],actionnames,action);  
  Cudd_RecursiveDeref(gbm,Array[0]);
  
  fclose(action);
  */
  
  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  
  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-stats.dat");
  strcat(basepath,templog);
  stats = fopen(basepath,"a");
  
  /* nodes in the Action ADD */
  dagCount = Cudd_DagSize(actionDD);
  dagLeafCount = Cudd_CountLeaves(actionDD);
  
  fprintf(dataFile2,"%d %d %d ",dagCount-dagLeafCount,dagLeafCount,dagCount);
  
  //Nodes in the equivalent value tree 
  treeCount = count_internal_nodes(Merged_Add);
  treeLeafCount = count_leaves(Merged_Add);

  fprintf(stats,"\n\n\nIterations to convergence %d\n",counter);
  fprintf(stats,"Final execution time: %8.4f  seconds\n\n\n\n\n",temps);
  fprintf(stats,"Memory usage presently: %ld bytes\n",usememory);
  fprintf(stats,"Maximum memory usage: %ld bytes\n\n\n\n",maxusememory);
  fprintf(stats,"Number of nodes in the action ADD: %d  internal nodes   %d leaves %d total nodes\n",
	  dagCount-dagLeafCount,dagLeafCount,dagCount);
  fprintf(stats,"Number of nodes in the equivalent tree: %d internal nodes  %d leaves  %d total nodes\n\n",
	  treeCount,treeLeafCount,treeCount+treeLeafCount); 
  
  // Nodes in the Value ADD
  dagCount = Cudd_DagSize(Merged_Add);
  dagLeafCount = Cudd_CountLeaves(Merged_Add);
  fprintf(stats,"Number of nodes in the value ADD: %d  internal nodes   %d leaves %d total nodes\n",
	   dagCount-dagLeafCount,dagLeafCount,dagCount);   
  fprintf(dataFile2,"%d %d %d ", dagCount-dagLeafCount,dagLeafCount,dagCount);   
  fprintf(dataFile2,"%d ",counter);
  fprintf(dataFile2,"%10.4f ",temps);
  //double max_error = get_error(MergedMin_Add);
  //double max_extent = get_extent(MergedMin_Add);
  fprintf(dataFile2,"0.0 0.0 ");

  //Nodes in the equivalent value tree 
  treeCount = count_internal_nodes(Merged_Add);
  treeLeafCount = count_leaves(Merged_Add);
  fprintf(stats,"Number of nodes in the value tree: %d internal nodes  %d leaves  %d total nodes\n\n",
	  treeCount,treeLeafCount,treeCount+treeLeafCount);
  
  Cudd_PrintInfo(gbm,stats);
  fclose(stats);
  
  // global OptimalPolicy
  OptimalPolicy = actionDD;
  return Merged_Add;
  //  return actionDD;
}
