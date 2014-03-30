#include "MDP.h"


// reads in a policy from filename
// queries the policy with state
// returns the action or -1 if fail
int mdpquery(char *filename, int *state)
{
  int optact(-1);
  DdNode *act, *val;

  DdNode **UPList= (DdNode **)malloc(2*(sizeof(DdNode *)));

  
  MDP *mdp = new MDP;

  //load from this file into dd: act and val
  int numread = mdp->readDualOptimal(gbm,filename,&UPList);
 
  if (numread == 0) {
    fprintf(stderr,"error reading dual policy from %s\n",filename);
    exit(1);
  } else {
    act = UPList[0];
    Cudd_Ref(act);

    val = UPList[1];
    Cudd_Ref(val);
    free(UPList);
  }

  optact = mdp->consultOptimalPolicy(act,val,state);
  
  return optact;
}


// test binary for mdpquery function
int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr,"usage: pquery <dualADD filename>\n");
    exit(1);
  }
  char * policyFile = *++argv;
  
  int nvars = 6;
  int *vvals = new int[nvars];
  for (int i=0; i<nvars; i++) 
    vvals[i] = 1;
  
  int gotact = mdpquery(policyFile,vvals);
  fprintf(stderr,"action was %d\n",gotact);
  return 1;
}
