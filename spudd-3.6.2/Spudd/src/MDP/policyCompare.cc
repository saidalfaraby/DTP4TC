#include "MDP.h"
#include <stdlib.h>
#include <stdio.h>

void getMaxDiff(DdNode *dd1, DdNode *dd2, double & mmax, double & mmin) {
  DdNode *ttmp1, *ttmp2;
  ttmp1 = Cudd_addApply(gbm,Cudd_addMinus,dd1,dd2);
  Cudd_Ref(ttmp1);
  mmax = ((*Cudd_V(Cudd_addFindMax(gbm,ttmp1))).get_max());
  mmin = ((*Cudd_V(Cudd_addFindMin(gbm,ttmp1))).get_max());
}

int main(int argc, char **argv) {
  DdNode *act1, *val1, *act2, *val2;
  int i;
  char basepath[256];
  int bigadd(10000);
  
  if (argc < 2) {
    fprintf(stderr,"usage: pquery <MDP file> <dualADD filename> <dualADD filename>\n");
    exit(1);
  }
  // input MDP
  char *infile = *++argv;
  // two Dual ADD files
  char * policyFile1= *++argv;
  char * policyFile2= *++argv;


  // read in valueDD, actionDD, ovars, vars
  DdNode **UPList1 =  (DdNode **)malloc(2*(sizeof(DdNode *)));
  DdNode **UPList2 =  (DdNode **)malloc(2*(sizeof(DdNode *)));
  
  sprintf(basepath,"%s",infile);

  MDP *mdp = new MDP(basepath,bigadd);
  mdp->classifySupport();

  // load from this file into dd: act[numpolicies] and val[numpolicies]
  int numread = mdp->readDualOptimal(gbm,policyFile1,&UPList1);
  numread = mdp->readDualOptimal(gbm,policyFile2,&UPList2);
  
  if (numread == 0) {
    fprintf(stderr,"error reading dual policy from %s\n",policyFile1);
    exit(1);
  } else {
    act1 = UPList1[0];
    Cudd_Ref(act1);
    act2 = UPList2[0];
    Cudd_Ref(act2);

    val1 = UPList1[1];
    Cudd_Ref(val1);
    val2 = UPList2[1];
    Cudd_Ref(val2);

    free(UPList1);
    free(UPList2);
  }
  // now compare the two
  DdNode *actdiff = Cudd_addApply(gbm,Cudd_addMinus,act1,act2);
  Cudd_Ref(actdiff);
  Cudd_PrintDebug(gbm,actdiff,4,100);
  double mmax, mmin;
  getMaxDiff(val1,val2,mmax,mmin);
  fprintf(stderr,"max diff is %f\n min diff is %f\n",mmax,mmin);
  return 1;
}

