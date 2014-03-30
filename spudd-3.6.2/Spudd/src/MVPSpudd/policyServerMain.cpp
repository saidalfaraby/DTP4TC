#include "pspudd.hh"
// these are in the homer project
#include "net2.h"
#include "networkSettings.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) {
  DdNode *act, *val;
  int i;

  if (argc < 2) {
    fprintf(stderr,"usage: pquery <dualADD filename>\n");
    exit(1);
  }
  char * policyFile= *++argv;


  // read in valueDD, actionDD, ovars, vars
  DdNode **UPList =  (DdNode **)malloc(2*(sizeof(DdNode *)));
  
  gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,MAXMEM);

  // load from this file into dd: act[numpolicies] and val[numpolicies]
  int numread = readDualOptimal(gbm,policyFile,&UPList);
  
  if (numread == 0) {
    fprintf(stderr,"error reading dual policy from %s\n",policyFile);
    exit(1);
  } else {
    act = UPList[0];
    Cudd_Ref(act);

    val = UPList[1];
    Cudd_Ref(val);
    free(UPList);
  }
  
  // start server
  policyServe(gbm,act,val, numorigvars);
  
  return 1;
}

