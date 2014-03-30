#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>

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
  char *infile, *outpath, *inpath, *obsprobfile, *logfile;
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

  int numb;

  outputfile = NULL;
  
  if (argc >= 3) {
    infile = *++argv;
    obsprobfile = *++argv;
    numb = atoi(*++argv);
  }
  // log file to get the final state values only
  if (argc >= 5) {
    logfile = *++argv;
  } else {
    logfile = NULL;
  }
  sprintf(basepath,"%s",infile);
  // also reads in the POMDP
  POMDP *mdp = new POMDP(basepath);

  // or from the timer
  struct timeval sd;
  int rseed;
  gettimeofday(&sd,NULL);
  rseed = (int) sd.tv_sec;
  //rseed = 1094154109;
  rseed = 1094397105;
  srand(rseed);
  fprintf(stderr,"seeded random number generator with %d\n",rseed);
  //mdp->generatePolicy(true);

  //FILE *fid = fopen("/tmp/pfrac.dat","a");
  sprintf(basepath,"%s_bmon_nowf.dat",obsprobfile);
  
  double pfrac = mdp->simulateHandwashing(obsprobfile,basepath,numb,logfile);
  //  fprintf(fid,"%f\n",pfrac);
  //fclose(fid);
}




