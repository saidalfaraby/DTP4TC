/***********************************************************************
  JAPRICODD: Approximate value iteration with ADDS
  $Author: jhoey $
  $Date: 2003/06/26 00:18:47 $
  $Revision: 2.4 $
  $Log: JAPRICODD.cc,v $
  Revision 2.4  2003/06/26 00:18:47  jhoey
  *** empty log message ***

  Revision 2.3  2003/06/25 21:47:29  jhoey
  *** empty log message ***

  Revision 2.2  2003/06/20 03:40:53  jhoey
  *** empty log message ***

  Revision 2.1  2003/01/31 00:38:58  jhoey
  removed some printing statements
  changed policy query a bit

  Revision 2.0  2003/01/16 22:43:58  jhoey
  got everything working now with reading/writing, mutli-valued vars
  it also works with the original data style

  Revision 1.7  2003/01/11 00:41:52  jhoey
  modified the policy generation thing so it doesn't keep all the Q
  functions every iteration of vi.
  also got adding a new policy to the policyQuery gui to work
  still have to make the makefile

  Revision 1.6  2002/12/14 02:04:13  jhoey
  *** empty log message ***

  Revision 1.5  2002/12/13 23:11:00  jhoey
  fixed round off method

  Revision 1.4  2002/12/12 21:55:33  jhoey
  *** empty log message ***

  Revision 1.2  2002/12/10 17:54:16  jhoey
  fixed a few bugs

  Revision 1.1  2002/12/05 06:43:06  jhoey
  added mutli-valued variables to PSpudd

  Revision 1.3  2002/09/25 20:21:21  jhoey
  CLK_TCK fixed to comply with new standard

  Revision 1.2  2002/09/21 00:55:14  jhoey
  added action costs

  Revision 1.1  2001/04/06 23:30:06  staubin
  testing....

  Revision 1.4  2001/03/21 04:36:14  jhoey
  JH: changed spanPair to span for terminal typedef stuff

  Revision 1.3  2000/10/05 21:26:08  jhoey
  added cvs authorship stuff at header


   Code by: Jesse Hoey
            Robert St-Aubin
***********************************************************************/
#include "pspudd.hh"
#include <limits.h>
DdNode *MergedMin_Add;
DdNode *VMinpast;
DdNode **VcurrentMin;
DdNode *RewardMax;
DdNode *nodeOne,*nodeTwo;
Pair tempOne,tempTwo;
double get_extent(DdNode *x) 
{
  double error;
  DdNode *mx,*mn;
  mx = Cudd_addFindMax(gbm,x);
  Cudd_Ref(mx);
  mn = Cudd_addFindMin(gbm,x);
  Cudd_Ref(mn);
  
  error = ((*(Cudd_V(mx))).get_max() - (*(Cudd_V(mn))).get_min());
  Cudd_RecursiveDeref(gbm,mx);
  Cudd_RecursiveDeref(gbm,mn);
  return error;
}

double get_error(DdNode *x){

  DdNode *errorAdd,*maxNode;
  double errorD;

  errorAdd = Cudd_addApply(gbm,getErrorAdd,x,One);
  Cudd_Ref(errorAdd);
  maxNode = Cudd_addFindMax(gbm,errorAdd);
  Cudd_Ref(maxNode);
  Cudd_RecursiveDeref(gbm,errorAdd);
  errorD = (*Cudd_V(maxNode)).get_max();
  //printf("Error:%f\n",errorD);
  Cudd_RecursiveDeref(gbm,maxNode);
  return errorD;
}

DdNode *getErrorAdd(DdManager *dd, DdNode **f, DdNode **x)
{
  Pair error;
  double errorD;
  DdNode *res;
  
  if (Cudd_IsConstant(*f)){
    errorD = (*Cudd_V(*f)).get_max() -(*Cudd_V(*f)).get_min();
    //printf("errorD:%f\n",errorD);
    error.set(errorD);
    res = Cudd_addConst(gbm,&error);
    return(res);
  } 
  return (NULL);     
} 

Pair get_span(DdNode *x) {
  Pair *res;
  res = (Pair *)malloc(sizeof(Pair));
  (*res).set_max((*Cudd_V(Cudd_addFindMax(gbm,x))).get_max());
  (*res).set_min((*Cudd_V(Cudd_addFindMin(gbm,x))).get_min());
  return *res;
}

double length_span(DdNode *x,DdNode *y){
  double length;
  double min,max;

  if((*Cudd_V(x)).get_min() <= (*Cudd_V(y)).get_min())
    min = (*Cudd_V(x)).get_min(); 
  else
    min = (*Cudd_V(y)).get_min(); 

  if((*Cudd_V(x)).get_max() <= (*Cudd_V(y)).get_max())
    max = (*Cudd_V(y)).get_max(); 
  else
    max = (*Cudd_V(x)).get_max(); 
      
  length = max-min;
  return length;
}

int close_enough(double val, double limit, double tol) 
{
  if(fabs(limit-val) < tol)
    return 1;
  else 
    return 0;
}

DdNode * My_cuddAddApplyRecur(DdNode * (*op)(DdManager *, DdNode **, DdNode **), 
			      DdNode * f, DdNode * g)
{
  DdManager *dd = gbm;
  DdNode *res,*fv, *fvn, *gv, *gvn, *T, *E;
  unsigned int ford, gord;
  unsigned int index;
  
  /* Check terminal cases. Op may swap f and g to increase the
   * cache hit ratio.
   */
  res = (*op)(dd,&f,&g);
  if (res != NULL) return(res);
  
  /* Recursive Step */
  ford = Cudd_ReadPerm(dd,f->index);
  gord = Cudd_ReadPerm(dd,g->index);
  if (ford <= gord) {
    index = f->index;
    fv = Cudd_T(f);
	fvn = Cudd_E(f);
    } else {
	index = g->index;
	fv = fvn = f;
    }
    if (gord <= ford) {
	gv = Cudd_T(g);
	gvn = Cudd_E(g);
    } else {
	gv = gvn = g;
    }

    T = My_cuddAddApplyRecur(op,fv,gv);
    if (T == NULL) 
      return(NULL);
    Cudd_Ref(T);

    E = My_cuddAddApplyRecur(op,fvn,gvn);
    if (E == NULL) {
      Cudd_RecursiveDeref(dd,T);
      return(NULL);
    }
    Cudd_Ref(E);

    //res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    res = Cudd_addIte(gbm,vars[(int) (index/2)].add_var,T,E);

    if (res == NULL) {
      Cudd_RecursiveDeref(dd, T);
      Cudd_RecursiveDeref(dd, E);
      return(NULL);
    }
    Cudd_Deref(T);
    Cudd_Deref(E);

    return(res);

} /* end of cuddAddApplyRecur */
DdNode *
My_addApply( DdNode * (*op)(DdManager *, DdNode **, DdNode **),DdNode * f,  DdNode * g)
{
    DdNode *res;
    DdManager *dd = gbm;
    res = My_cuddAddApplyRecur(op,f,g);
    return(res);

} /* end of Cudd_addApply */

/* 
converts the ADD f to a 0-1 ADD by setting every leaf in f 
which is equal to the constant c to 1, and everything else to 0 
For use with addApply
*/
DdNode *addAddExact(DdManager *dd, DdNode **f, DdNode **c)
{
  DdNode *res;
  
  if (Cudd_IsConstant(*f) && Cudd_IsConstant(*c)) {
    if (Cudd_V(*f)== Cudd_V(*c))
      res = Cudd_addConst(dd,pOne);
    else
      res = Cudd_addConst(dd,pZero);
    return (res);
  } 
  return (NULL);
}

DdNode *getTotalSpan(DdManager *dd, DdNode **f, DdNode **x)
{
  Pair pRes;
  //double pRes;
  DdNode *dMax, *dMin;
  
  
  if (Cudd_IsConstant(*f)){
    if ((*Cudd_V(*f)) == (*pOne)) {
      // *pRes = get_span(*x);
      dMax = Cudd_addFindMax(gbm, *x);
      Cudd_Ref(dMax);
      dMin = Cudd_addFindMin(gbm, *x);
      Cudd_Ref(dMin);
      //(*pRes).set_max((*Cudd_V(Cudd_addFindMax(gbm,*x))).get_max());
      //(*pRes).set_min((*Cudd_V(Cudd_addFindMin(gbm,*x))).get_min());
      pRes.set_max((*(Cudd_V(dMax))).get_max());
      pRes.set_min((*(Cudd_V(dMin))).get_min());

      Cudd_RecursiveDeref(gbm,dMax);
      Cudd_RecursiveDeref(gbm,dMin);
      avgError.span(pRes);
      return(Cudd_addConst(dd,&pRes));
      //delete(pRes);
      //return(One);
      //return(Cudd_addConst(dd,pZero));
    }
    //return(Cudd_addConst(dd,pZero));
    return(Zero);
  } 
  return (NULL);     
} 

DdNode *roundOffMMA(DdNode *res)
{
  DdGen *gen;
  DdNode *node, *tmpres, *tmp, *tmp2, *dRes, *thisLeaf;
  Pair tmpError;
  
  
  //tmpError = new Pair();

  dRes = Zero;
  Cudd_Ref(dRes);

  Cudd_ForeachNode(gbm,res,gen,node) {
    if (Cudd_IsConstant(node)) {
      avgError.set_max(-1.0*LARGEFLAG);
      avgError.set_min(1.0*LARGEFLAG);
      
      
      tmp = Cudd_addApply(gbm, addAddExact, res, node);
      Cudd_Ref(tmp);
     
      //This is only done for its effect (to get avgError)
      tmp2 = My_addApply(getTotalSpan, tmp, MergedMin_Add);
      Cudd_Ref(tmp2);
      
      tmpError = avgError;
      thisLeaf = Cudd_addConst(gbm,&tmpError);
      Cudd_Ref(thisLeaf);
      
      Cudd_RecursiveDeref(gbm,tmp2);
      tmp2 = Cudd_addApply(gbm, Cudd_addTimes, thisLeaf, tmp);
      Cudd_Ref(tmp2);
      
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,thisLeaf);
       
      tmp = Cudd_addApply(gbm, Cudd_addPlus, dRes, tmp2);
      Cudd_Ref(tmp);
      Cudd_RecursiveDeref(gbm,tmp2);
      Cudd_RecursiveDeref(gbm,dRes);
      dRes = tmp;
    }
  }
  return (dRes);
}

void size_approx(int limit_type) 
{
  DdNode *res, *temp, *tmp;
  DdNode *dRoe, *dLFlag, *dLFlag2, *dAvg, *z, *w;
  double oldroe = 0.0, olderoe;
  
  int count, oldcount;
  int iter1(0),iter2(0),maxiters(100);
  Pair pRoe;
    
  //just to make sure that it is bigger
  double roe = get_extent(MergedMin_Add)*2.0;

  if (limit_type == LIMIT_ERROR) {
    pRoe.set(Max_Error);
    fprintf(stderr,"Limit Error - Max_Error: %f \n",Max_Error);
  } else if (limit_type == LIMIT_SIZE) {
    pRoe.set(roe);
    // double precision digit limit
    // this should be a little bigger than unique->epsilon in the dd manager
    Max_Error = 2.0*(Cudd_ReadEpsilon(gbm)->get_min());
    //fprintf(stderr,"Limit Size - Max_Size: %d\n",Max_Size);
  }

  dRoe = Cudd_addConst(gbm,&pRoe);
  Cudd_Ref(dRoe);

  cntError = 0;

  res = Cudd_addApply(gbm, myRoundOff,MergedMin_Add, dRoe);
  Cudd_Ref(res);
  count = Cudd_DagSize(res);
  
  oldcount = -1;
  if (limit_type == LIMIT_ERROR) {
    /* replace MergedMin_Add with rounded off version */
    tmp = roundOffMMA(res);
    Cudd_RecursiveDeref(gbm,MergedMin_Add);
    MergedMin_Add = tmp;
    Cudd_RecursiveDeref(gbm,res);
  } else if (limit_type == LIMIT_SIZE) {

    /* while ADD is too small       - decrease the roe - increase the ADD size */
    // the maxiters thing is only to avoid infite loops which used to happen - 
    // probably should not be there
    oldroe = 0.0;
    while (iter1 < maxiters && count < Max_Size && (fabs(oldroe-roe) > Max_Error)) {
      oldcount = count;
      //fprintf(stderr,"count %d Max_Size %d  error %e max_error %e oldroe %e roe %e ",
      //      count,Max_Size,fabs(oldroe-roe),Max_Error,oldroe,roe);
      roe = roe - (roe - oldroe)/2.0;
      //fprintf(stderr,"new roe %e\n",roe);
      pRoe.set(roe);
      
      Cudd_RecursiveDeref(gbm,dRoe);
      dRoe = Cudd_addConst(gbm,&pRoe);
      Cudd_Ref(dRoe);
      
      cntError = 0;
      avgError.set(0.0);
      
      Cudd_RecursiveDeref(gbm,res);
      res = Cudd_addApply(gbm, myRoundOff,MergedMin_Add, dRoe);
      Cudd_Ref(res);
      
      count = Cudd_DagSize(res);
      if (count > Max_Size) 
	oldroe = roe*2.0 - oldroe;
      
      /* while the ADD is too big - increase the roe - decrease the ADD size */
      iter2=0;
      while (iter2 < maxiters && count > Max_Size & (fabs(oldroe-roe) > Max_Error)) {
	/* replace MergedMin_Add with rounded off version 
	   since we know it will be smaller than this */
	//fprintf(stderr,"--------count %d Max_Size %d  error %e max_error %e oldroe %e roe %e ",
	//	count,Max_Size,(oldroe-roe)/2.0,Max_Error,oldroe,roe);
	
	roe = roe + (oldroe - roe)/2.0;
	pRoe.set(roe);
	//fprintf(stderr,"new roe %f\n",roe);
	
	Cudd_RecursiveDeref(gbm,dRoe);
	dRoe = Cudd_addConst(gbm,&pRoe);
	Cudd_Ref(dRoe);


	Cudd_RecursiveDeref(gbm,res);
	res = Cudd_addApply(gbm, myRoundOff,MergedMin_Add, dRoe);
	Cudd_Ref(res);
 
	
	count = Cudd_DagSize(res);
	if (count <= Max_Size) {
	  olderoe = roe;
	  oldroe = roe*2.0 - oldroe;
	  //fprintf(stderr,"**** resetting count %d Max_Size %d  oldroe %e roe %e\n ",
	  //  count,Max_Size,oldroe,roe);	
	}
	else if (fabs(roe-oldroe) < Max_Error) {
	  roe = olderoe;
	  while (count > Max_Size) {
	    //fprintf(stderr,"?????????? count %d Max_Size %d  error %e max_error %e oldroe %f roe %f ",
	    //	    count,Max_Size,(oldroe-roe)/2.0,Max_Error,oldroe,roe);
	    pRoe.set(roe);
	    
	    Cudd_RecursiveDeref(gbm,dRoe);
	    dRoe = Cudd_addConst(gbm,&pRoe);
	    Cudd_Ref(dRoe);
	    
	    Cudd_RecursiveDeref(gbm,res);
	    res = Cudd_addApply(gbm, myRoundOff,MergedMin_Add, dRoe);
	    Cudd_Ref(res);
	    count = Cudd_DagSize(res);
	    roe += tolerance;
	  }
	  oldroe = roe;
	}
	iter2++;
	iter1=0;
      }
    }
    //fprintf(stderr,"roundoffMMA\n");
    tmp = roundOffMMA(res);
    Cudd_RecursiveDeref(gbm,MergedMin_Add);
    MergedMin_Add = tmp;
    Cudd_RecursiveDeref(gbm,res);
  }
}


/*expects two constants f,g, rounds off f to within trunc and returns the result */ 

DdNode *myRoundOff(DdManager *dd, DdNode **f, DdNode **roe)
{
  DdNode *tmp, *res, *F, *ROE;
  Pair value;

  F = *f; ROE = *roe;
  //divides so that it is between integers
  if (Cudd_IsConstant(F) && Cudd_IsConstant(ROE)) {
    tmp = Cudd_addDivide(dd, f, roe);
    Cudd_Ref(tmp);
    
    //since we are using integer round-off
    if ((floor((*Cudd_V(tmp)).get_min())) == (floor((*Cudd_V(tmp)).get_max()))) {
      avgError.span(*Cudd_V(F));
      cntError++;
      // why was this LARGEFLAG thing here????
      //value.set_min(2*LARGEFLAG + floor((*Cudd_V(tmp)).get_min()));
      //value.set_max(2*LARGEFLAG + ceil((*Cudd_V(tmp)).get_max()));
      value.set_min(floor((*Cudd_V(tmp)).get_min()));
      value.set_max(ceil((*Cudd_V(tmp)).get_max()));
      /* deal with the fence-sitters*/
      if(value.get_min() == value.get_max())
	value.set_max(value.get_max()+1);
      Cudd_RecursiveDeref(dd, tmp);
      res = Cudd_addConst(dd,&value);
      return (res);
    } else {
      Cudd_Deref(tmp);
      return (tmp);                // This added late thursday night in a rush - right?
    }
  }
  return (NULL);
}

/* 
   Convergence test to be used in conjuncture with Cudd_addApply
   Determines if the pair ADD has converged.  
   We are using a very conservative criterion: if all intervals overlap
   or are within tolerance eps of each other, we have converged.
*/

DdNode *Convergence_Test(DdManager * gbm, DdNode ** f, DdNode ** g)
{
 
  DdNode *F, *G;
  DdNode *oneD,*zeroD;
  F = *f; G = *g;
    
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    double min1 = (*Cudd_V(F)).get_min();
    double min2 = (*Cudd_V(G)).get_min();
    double max1 = (*Cudd_V(F)).get_max();
    double max2 = (*Cudd_V(G)).get_max();
    
    /*for the case where there is no approximation*/
    /*when there is approximation being done*/
    if( ((min1 <= min2) && (min2 <= max1)) || ((min2 <= min1) && (min1 <= max2)) || ((min2 >= max1) && ((min2 - max1) <= tolerance)) || ((max2 <= min1) && ((min1 - max2) <= tolerance))){
      oneD = Cudd_addConst(gbm,pOne);
      return (oneD);
    }else{
      zeroD = Cudd_addConst(gbm,pZero);
      return(zeroD);
    } 
  } 
  return(NULL);
} 
int reorderMinSpan(DdNode *addBase)
{
  double *infor;
  double s,st,se;
  int i,nt,ne,res;
  int *phase;
  int *reordList;
  DdNode *vRest, *oCube;
  DdNode *cVars;

  infor = (double *) malloc(numvars*sizeof(double));
  reordList = (int *) malloc(2*numvars*sizeof(int));
  phase = (int *) malloc(sizeof(int));

  for (i=0; i<2*numvars;i++)
    reordList[i] = i;
  
  s = get_extent(addBase);

  for (i=0;i<numvars;i++) {

    cVars = vars[i].add_var;
    Cudd_Ref(cVars);

    /* set this variable to 1 */
    phase[0] = 1;
    oCube = Cudd_addComputeCube(gbm,&cVars,phase,1);
    Cudd_Ref(oCube);


    vRest = Cudd_addRestrict(gbm,addBase,oCube);
    Cudd_Ref(vRest);
    Cudd_RecursiveDeref(gbm,oCube);

    st = get_extent(vRest);
    nt = Cudd_CountLeaves(vRest);
    Cudd_RecursiveDeref(gbm,vRest);

    /* set it to 0 */
    phase[0] = 0;
    oCube = Cudd_addComputeCube(gbm,&cVars,phase,1);
    Cudd_Ref(oCube);
    Cudd_RecursiveDeref(gbm,cVars);

    vRest = Cudd_addRestrict(gbm, addBase,oCube);
    Cudd_Ref(vRest);
    Cudd_RecursiveDeref(gbm,oCube);
    se = get_extent(vRest);
    ne = Cudd_CountLeaves(vRest);
    Cudd_RecursiveDeref(gbm,vRest);
    
    //infor[i] = (2*s-st-se)*nt/(s*(nt+ne));
    //infor[i] = -1.0*(log(1-st/s)+log(1-se/s))*nt/(nt+ne);
    infor[i] = 0;
    if (st > 0) 
      infor[i] += (nt/nt+ne)*log(st/s);
    else
      infor[i] -= 10000.0;
    if (se > 0)
      infor[i] += (ne/nt+ne)*log(se/s);
    else
      infor[i] -=  10000.0;

  }
  /* sort infor list and build reorder list (indices of variables which had
     smallest errors */

  infoSorter(infor,reordList);
  res = Cudd_ShuffleHeap(gbm,reordList);
  free(infor);
  free(reordList);
  free(phase);
  return(res);
  
}

void infoSorter(double *infoList, int *reorderList) 
{
  int i,j,minj;
  double maxI=-10000.0, minI;
  /* find maximum value for bookeeping*/
  for (i=0; i<numvars; i++) 
    if (infoList[i] > maxI)
      maxI = infoList[i];
 

  for (i=0; i<numvars; i++) {

    /* find minimum value */
    minI = maxI + 1.0;
    for (j=0; j<numvars; j++) 
      if (infoList[j] < minI) {
	minj = j;
	minI = infoList[j];
      }

    /* get this value out of the way*/
    infoList[minj] = maxI+1.0;

    /*update the unprimed variable order */

    reorderList[2*i+1] = 2*minj+1;
    reorderList[2*i] = 2*minj;
  }
}
//takes in a ADD and a constant ADD node and replaces the 2 nodes that I want by the ADD constant
DdNode *replacePairs(DdManager *dd, DdNode **f, DdNode **newC)
{
  DdNode *tmp, *res, *F, *ROE;
  double min_span,min_temp;
  
  F = *f; ROE = *newC;
  //Look for a constant that equals 
  //PairOne and PairTwo are the 2 constants that need to be replaced by the newC constant node
  if (Cudd_IsConstant(F)){
    if( ((*Cudd_V(F)) == tempOne) || ((*Cudd_V(F)) == tempTwo) ) {
      //if( F == nodeOne || F == nodeTwo ) {
      return ROE;
    }else{
      return F;
    } 
  }else{
    return(NULL);
  }
}

//removes an element from the constant list
void removeElement(list_cst *list,list_cst *element){
  list_cst *curr,*suivant;
  
  curr = list;
  suivant = list->next;
  while(suivant != NULL){
    if(suivant == element)
      {
	curr->next = suivant->next;
	Cudd_RecursiveDeref(gbm,suivant->add);
	free(suivant);
	break;
      }else{
	curr = suivant;
	suivant = suivant->next;
      }
  }
}
int getListSize(list_cst *list) {
  list_cst *curr;
  int cnt=0;
  curr = list;
  while (curr != NULL) {
    cnt++;
    curr = curr->next;
  }
  return cnt;
}
void printList(list_cst *list) {
  list_cst *curr;
  curr = list;
  while (curr != NULL) {
    fprintf(stderr,"%s\n",(*Cudd_V(curr->add)).toString());
    curr = curr->next;
  }
}
void printAllNodes(DdNode *mmadd) {
  DdGen *gen;
  DdNode *node;
  
  Cudd_ForeachNode(gbm,MergedMin_Add,gen,node) {
    if (Cudd_IsConstant(node)) { 
      fprintf(stderr,"%s\n",(*Cudd_V(node)).toString());
    }
  }
}

//This approximation method finds the 2 pairs with minimal span and approximates them
// by the global span of the two
void allPairsApprox(int limit_type){
  DdGen *gen;
  DdNode *node,*Approx,*temp;
  list_cst *const_list,*ins,*curr,*suivant;
  double min,max,min_span,min_temp;
  double current_min;
  list_cst *elem1, *elem2;
  int count;
  int typeFlag;

  //const_list  is a pointer to a constant node list 
  //At first the list is Null
  const_list = NULL;
  
  //We need to find all the constant nodes so that we include them in the constant node list
  Cudd_ForeachNode(gbm,MergedMin_Add,gen,node) {
    if (Cudd_IsConstant(node)) { // I want to put the pointer in a list by order of minimum
      
      //we allocate the memory for the new constant node in the list
      ins = (list_cst *)malloc(sizeof(list_cst));
      ins->next = NULL;
      ins->add = node;

      Cudd_Ref(ins->add);

      //fprintf(stderr,"%x %d %16.14lf %16.14lf\n",ins->add,ins->add->ref,(*Cudd_V(node)).get_min(),(*Cudd_V(node)).get_max());


      //we now insert in the list
      //we begin at the start of the list
      curr = const_list;
      suivant = const_list;
      
      //if the list is empty, we put this element at the beginning
      if (curr == NULL)
	{
	  const_list = ins;
	}
      //we put this element at his good position:(increasing order of the min)
      else{
	suivant = curr->next;
	//if the new constant has a smaller min than the current one, put it at the start of the list.
	if( (*Cudd_V(curr->add)).get_min() >= (*Cudd_V(ins->add)).get_min()){
	  ins->next = curr;
	  const_list = ins;
	}else{
	  //we put the new constant at the good position in the list
	  while(suivant !=NULL){
	    if((*Cudd_V(suivant->add)).get_min() < (*Cudd_V(ins->add)).get_min())
	      {
		curr = suivant;
		suivant = suivant->next;
	      }else{
		break;
	      }
	  }
	  ins->next = suivant;
	  curr->next = ins;
	}
      }
    }
  }
  //we are done ordering the constant nodes by increasing size of their min.
  //now we need to find which pairs to approximate.
  if(const_list == NULL)
    printf("There is a problem with my list.....\n");
  
    
  
  nodeOne = Zero;
  Cudd_Ref(nodeOne);
  nodeTwo = Zero;
  Cudd_Ref(nodeTwo);
  //Need to optimize it so that I don't run through the whole list n^2 time......
  //We run through all the constant node to find the minimum span.
  /*############################

    MIGHT WANT TO KEEP MORE THAN ONE MINIMUM.....TO LOOK INTO MAYBE

    #############################

  */

  min_span = 0.0;
  count = Cudd_DagSize(MergedMin_Add);
  if(limit_type == LIMIT_ERROR)
    typeFlag = 1;
  else if(limit_type == LIMIT_SIZE)
    typeFlag = 0;
  int shitty = 0;
  while ( ((count > Max_Size) && (typeFlag == 0)) || (((min_span/2.0) <= Max_Error) && typeFlag == 1)) {
    shitty++;
    //Set the span_min to a big number so that I find a new best span
    min_span = DBL_MAX;
    curr = const_list;
    elem1 = NULL;
    elem2 = NULL;
    while(curr != NULL){
      current_min = (*Cudd_V(curr->add)).get_min(); 
      suivant = curr->next;

      while(suivant != NULL){
	if( ((typeFlag == 1) && ((2*Max_Error) > (*Cudd_V(suivant->add)).get_min() - current_min)) || 
	    ((typeFlag == 0) && (min_span > (*Cudd_V(suivant->add)).get_min() - current_min)))  {
	  if(min_span > (min_temp = length_span(curr->add,suivant->add)))
	    {
	      min_span = min_temp;
	      elem1 = curr;
	      elem2 = suivant;
	    }
	  suivant = suivant->next;
	}else{
	  break;
	}
      }
      curr = curr->next;
    }
    if(((min_span/2.0) > Max_Error) && (typeFlag ==1))
      break;
    // may not have found anything so 
    /*
    if (elem1 == NULL || elem2 == NULL) {
      fprintf(stderr,"HEY HEY HEY HEY!");
      break;
    }
    */
    Cudd_RecursiveDeref(gbm,nodeOne);
    nodeOne = elem1->add;
    //fprintf(stderr,"nodeTwo's ref before dereffing is %d\n\n",nodeTwo->ref);
    Cudd_Ref(nodeOne);
    Cudd_RecursiveDeref(gbm,nodeTwo);
    //fprintf(stderr,"elem2's (%x) ref is %d elem1's ref is %d\n",elem2->add,elem2->add->ref,elem1->add->ref);
    nodeTwo = elem2->add;
    Cudd_Ref(nodeTwo);
    //fprintf(stderr,"nodeTwo's ref after reffing is %d\n",nodeTwo->ref);
    
    //We are done, now lets create the new constant that we will replace the two nodes with
    if((*Cudd_V(nodeOne)).get_min() <= (*Cudd_V(nodeTwo)).get_min())
      min = (*Cudd_V(nodeOne)).get_min(); 
    else
      min = (*Cudd_V(nodeTwo)).get_min(); 
    
    if((*Cudd_V(nodeOne)).get_max() <= (*Cudd_V(nodeTwo)).get_max())
      max = (*Cudd_V(nodeTwo)).get_max(); 
    else
      max = (*Cudd_V(nodeOne)).get_max(); 
    
    Pair newConstant;
    newConstant.set_min(min);
    newConstant.set_max(max);
    temp = Cudd_addConst(gbm,&newConstant);
    Cudd_Ref(temp);
  
    tempOne = *Cudd_V(nodeOne);
    tempTwo = *Cudd_V(nodeTwo);
    

    //Now we get the new Merged_Add 
    /*
    fprintf(stderr,"------- the add constants are\n");
    printAllNodes(MergedMin_Add);
    fprintf(stderr,"------- the list constants are\n");
    printList(const_list);
    fprintf(stderr,"------------------------------\n");
    fprintf(stderr,"replacing %s and %s with %s\n",tempOne.toString(),tempTwo.toString(),(*Cudd_V(temp)).toString());
    */
    Approx = Cudd_addApply(gbm,replacePairs,MergedMin_Add,temp);
    Cudd_Ref(Approx);
    Cudd_RecursiveDeref(gbm,MergedMin_Add);
    MergedMin_Add = Approx;


    //Cudd_Ref(MergedMin_Add);
    //Cudd_RecursiveDeref(gbm,Approx);

    count = Cudd_DagSize(MergedMin_Add);
    //int leafcount = Cudd_CountLeaves(MergedMin_Add);
    //fprintf(stderr,"leafcount is %d\n",leafcount);
    //Now update the list of constants.
    //We replace the element with the smallest min by the new span constant
    Cudd_RecursiveDeref(gbm,elem1->add);
    elem1->add = temp;
    
    //fprintf(stderr," ---------------- new elem1->add's (%x) new ref count is %d\n",elem1->add,elem1->add->ref);
    //Now we need to remove the other element from the list
    removeElement(const_list,elem2);
    /*
    fprintf(stderr,"------- the add constants are\n");
    printAllNodes(MergedMin_Add);
    fprintf(stderr,"------- the list constants are\n");
    printList(const_list);
    fprintf(stderr,"------------------------------\n");
    */
    //leafcount = getListSize(const_list);
    //fprintf(stderr,"nodes in list %d\n\n",leafcount);
  }
 }
DdNode *My_addMaximum(DdManager * dd, DdNode ** f, DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;

    if (F == DD_MINUS_INFINITY(dd)) return(G);
    if (G == DD_MINUS_INFINITY(dd)) return(F);

    if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
        if ((*Cudd_V(F)) >= (*Cudd_V(G))) {
          return(F);
        } else if ((*Cudd_V(G)) > (*Cudd_V(F))) {
          return(G);
        } else {
          Pair newPair;
          double mn,mx;
          double fmx = (*Cudd_V(F)).get_max();
          double fmn = (*Cudd_V(F)).get_min();
          double gmx = (*Cudd_V(G)).get_max();
          double gmn = (*Cudd_V(G)).get_min();
          if (fmx >= gmx)
            mx = fmx;
          else 
            mx = gmx;
          if (fmn >= gmn)
            mn = fmn;
          else
            mn = gmn;

          newPair.set_both(mn,mx);
          DdNode *res = Cudd_addConst(dd,&newPair);
          return (res);
          
        }
    }
    return(NULL);

} /* end of My_addMaximum */

DdNode *NewAspudd_Proc(int limit_type, int reorderMeth, int approxMeth, FILE *dataFile, FILE *dataFile2)
{
  int newcounter,test;
  int counter,approx_flag,flag;
  int i,j,k,nodeCount,dagCount,treeCount,treeLeafCount,dagLeafCount;
  int templist[MAXVARS];
  long int maxusememory=0,usememory=0;
  double global_max_error,max_error,tmp_error,tol;   
  char templog[MAXLINE];
  char basepath[MAXLINE];
  DdNode *MidPoint_Add,*temp1,*temp2,*tmep3,*actionAdd2;
  DdNode *oldMerged;
  DdNode **V0vector,*testD,*testing;
  double temps0,temps1,temps;   
  struct timeval *Time0;       /* Timestamps */
  struct timeval *Time1;

  /* Structure to get the time of execution*/
  struct tms _t;
  long clk_tck = sysconf(_SC_CLK_TCK);

  FILE *value,*action,*errf,*both,*stats,*log,*log2,*reward,*actions;

  int *varorder;
  int oldCount;
  varorder = (int*) malloc(2*numvars*sizeof(int));
  
  VcurrentMin =(DdNode **)malloc(numactions*sizeof(DdNode *));
  for(i=0;i<numactions;i++)
    {
      VcurrentMin[i] = RewardD;
      Cudd_Ref(VcurrentMin[i]);
    }
    
  // Vpast keep tracks of the previous value ADD, initialized at V0
  VMinpast = RewardD; 
  Cudd_Ref(VMinpast);   

  /* Allocate storage for timestamps */
  Time0 = (struct timeval*)malloc(sizeof(struct timeval));
  Time1 = (struct timeval*)malloc(sizeof(struct timeval));

  
  //Merged keep tracks of the current value ADD, initialized to RewardD
  MergedMin_Add = RewardD; 
  Cudd_Ref(MergedMin_Add); 

  //ADD RewardMax = RewardD;
  RewardMax = RewardD;
  Cudd_Ref(RewardMax);

 
  // Convergence test
  test = 0;
  
  //Number of iteration to convergence
  counter = 0;
  
  //counter for the number of approximation performed
  newcounter = 0;
  
  //Global approximation error
  global_max_error = 0.0;
  
  // Used in calculating the number of approximation performed
  approx_flag = 0;

  actionAdd2 = RewardD;
  Cudd_Ref(actionAdd2);

  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-Approx-stats.dat");
  strcat(basepath,templog);
  stats = fopen(basepath,"w");
  fprintf(stats,"\n\n\n Name of input file: %s\n\n Discount factor is : %f \n\n\nTolerance is: %f \n horizon is: %f \n The BIGADD limit is set to: %d \n  The MAX_SIZE is set to: %d\n\n Target maximum memory use: %d bytes \n\n Hard Limit on memory use %d bytes\n\n\n",infile, discount_factor, tolerance, horizon, bigadd,MAX_SIZE,MAXMEM, MAXMEMHARD);
  fprintf(stats,"\n\n");
  

  usememory = Cudd_ReadMemoryInUse(gbm);
  if (maxusememory < usememory)
    maxusememory = usememory;
  fprintf(stats,"Memory in use before start of iterations: %ld\n",usememory);
  fprintf(stats,"\n\n\n");
  fclose(stats);
    

   /* Get the starting up time */
  // gettimeofday(Time0, NULL );
   /* Get the starting up time */
  times(&_t);
  temps0 = 1.0*_t.tms_utime/clk_tck;
 

  
  // while the value ADD has not converged to the optimal value
  while ((horizon == -1 && test==0) || (horizon >= 0 && counter < horizon)){
    /* WE KEEP TRACK OF THE NUMBER OF ITERATIONS */
    counter = counter + 1;

    //fprintf(stderr,"iteration: %d\n",counter);
    //Cudd_PrintDebug(gbm,MergedMin_Add,3,3);

    Cudd_RecursiveDeref(gbm,VMinpast); 
    VMinpast = MergedMin_Add; 
    Cudd_Ref(VMinpast); 
     
    Cudd_RecursiveDeref(gbm,RewardMax); 
    RewardMax = Cudd_addSwapVariables(gbm,MergedMin_Add,Array1,Array2,numvars); 
    Cudd_Ref(RewardMax); 
    
    Cudd_RecursiveDeref(gbm,MergedMin_Add);
    MergedMin_Add = VerySmall;
    Cudd_Ref(MergedMin_Add);

    

    /* WE  COMPUTE THE ADDs FOR EACH VARIABLES WITH THEIR PRIMED ADDs AND THEN 
	 WE ADD EVERYTHING */

    for(i=0;i<numactions;i++) { 
#ifdef OPTIMIZED
      temp1 = primeReward(RewardMax,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
      temp1 = newPrimeReward(RewardMax,NewPrime[i]);
#endif
      temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
      Cudd_Ref(temp2);
      
      Cudd_RecursiveDeref(gbm,temp1);
      temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
      Cudd_Ref(temp1);
      
      Cudd_RecursiveDeref(gbm,temp2);
      
      // add the action costs. 
      temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp1,actionCost[i]);
      Cudd_Ref(temp2);

      Cudd_RecursiveDeref(gbm,temp1);
      
      // totally useless
      /*
      Cudd_RecursiveDeref(gbm,VcurrentMin[i]);
      VcurrentMin[i] = temp2;
      Cudd_Ref(VcurrentMin[i]);       
      */


      temp1=Cudd_addApply(gbm,My_addMaximum,MergedMin_Add,temp2);
      Cudd_Ref(temp1);
      
      Cudd_RecursiveDeref(gbm,temp2);
      Cudd_RecursiveDeref(gbm,MergedMin_Add);
      MergedMin_Add = temp1;      
    }

    /* Nodes in the Value ADD*/
    dagCount = Cudd_DagSize(MergedMin_Add);
    dagLeafCount = Cudd_CountLeaves(MergedMin_Add);
    tmp_error = get_error(MergedMin_Add);
    
    //fprintf(stderr,"nodes in value ADD after summation %d\n",dagCount);

    //Cudd_ReduceHeap(gbm,CUDD_REORDER_RANDOM,20);
    tol = -1.0;
    if (limit_type == LIMIT_ERROR && (counter > -1)) //We always prune 
      tol = ERROR_TOL;  //tmp_error > Max_Error) //&&	!close_enough(tmp_error,Max_Error,ERROR_TOL)) 
    else if (limit_type == LIMIT_SIZE && 
	     dagCount > Max_Size) // &&  !close_enough((double) dagCount, Max_Size, SIZE_TOL))
      tol = SIZE_TOL;
    //fprintf(stdout,"Dead nodes at end of loop - before approx: %d\n",Cudd_ReadDead(gbm));


    /************ NEW CODE */
    /* save old variable order - should be a list where the 
       ith entry is the index of the variable at the ith level */

    // approximate
    if (tol > 0) {

      approx_flag = 1;
 
      switch (reorderMeth)
	{
	case REORDER_SIFT:
	  if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_SIFT, 0))
	    fprintf(stderr,"*****ERROR: Reorder sift failed \n");
	  break;
	case REORDER_RANDOM:
	  if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_RANDOM, 0))
	    fprintf(stderr,"*****ERROR: Reorder random failed \n");
	  break;
	case REORDER_MINSPAN:
	  if (!reorderMinSpan(MergedMin_Add))
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

      /* Nodes in the Value ADD*/
      dagCount = Cudd_DagSize(MergedMin_Add);
      //      fprintf(stderr,"nodes in value ADD after reorder %d\n",dagCount);

      // We need to change tolerance if it is in the sliding mode
      switch (toleranceMode) 
	{
	case TOLERANCE_SLIDING:
	  Max_Error = prune*extR*((pow(discount_factor,((double)counter +1))-1)/(discount_factor-1));
	  //printf("we change the tolerance\n");
	  break;
	default:
	  break;
	}

      /* Nodes in the Value ADD*/
      dagCount = Cudd_DagSize(MergedMin_Add);
      max_error = get_error(MergedMin_Add);
      
      if(limit_type != LIMIT_SIZE || dagCount > Max_Size) {
	// Our two approximation techniques
	switch (approxMeth) 
	  {
	  case APPROX_ROUNDOFF:
	    size_approx(limit_type);
	    break;
	  case APPROX_ALLPAIRS:
	    allPairsApprox(limit_type);
	    break;
	  default:
	    break;
	  }
      }
      dagCount = Cudd_DagSize(MergedMin_Add);
      //      fprintf(stderr,"nodes in value ADD after approx %d\n",dagCount);
      
    }
    tmp_error= 0.0;

    usememory = Cudd_ReadMemoryInUse(gbm);

    if (maxusememory < usememory)
      maxusememory = usememory;
    if (usememory > MAXMEMHARD) {
      fprintf(errf,"Memory required exceeds availabe %ld bytes. \nTry a smaller bigadd limit constant\n",MAXMEMHARD);
      Cudd_Quit(gbm);
      exit(0);
    }
    
    /* WE COMPARE THE TWO MOST RECENT ADDs IN ORDER TO VERIFY THE CONVERGENCE  */
    
    testD = Cudd_addApply(gbm,Convergence_Test,VMinpast,MergedMin_Add);
    Cudd_Ref(testD);
    if((*Cudd_V(testD)).get_min() == 1.0)
      {
	Cudd_RecursiveDeref(gbm,testD);
	test = 1;
      }else{
	Cudd_RecursiveDeref(gbm,testD);
      } 
  }
  /* We get the time of completion*/
  times(&_t);
  temps1 = 1.0*_t.tms_utime/clk_tck;
  temps = temps1- temps0;
  for (i=0; i<2*numvars;i++)
    fprintf(dataFile,"%d ",Cudd_ReadInvPerm(gbm,i));

  /* Now back up once more after finding the midpoint ADD */

  MidPoint_Add = Cudd_addApply(gbm, addMean, MergedMin_Add, One);
  Cudd_Ref(MidPoint_Add);
  
  Cudd_RecursiveDeref(gbm,RewardMax); 
  RewardMax = Cudd_addSwapVariables(gbm,MergedMin_Add,Array1,Array2,numvars); 
  Cudd_Ref(RewardMax); 

  Cudd_RecursiveDeref(gbm,MidPoint_Add);
  MidPoint_Add = VerySmall;
  Cudd_Ref(MidPoint_Add);
  
  /* WE  COMPUTE THE ADDs FOR EACH VARIABLES WITH THEIR PRIMED ADDs AND THEN 
     WE ADD EVERYTHING */
  dagCount = Cudd_DagSize(RewardMax);
  for(i=0;i<numactions;i++) { 
#ifdef OPTIMIZED
    temp1 = primeReward(RewardMax,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
    temp1 = newPrimeReward(RewardMax,NewPrime[i]);
#endif

    temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
    Cudd_Ref(temp2);
    
    Cudd_RecursiveDeref(gbm,temp1);

    temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
    Cudd_Ref(temp1);   
    Cudd_RecursiveDeref(gbm,temp2);
    
    // add the action costs. 
    temp2 = Cudd_addApply(gbm,Cudd_addPlus,temp1,actionCost[i]);
    Cudd_Ref(temp2);

    Cudd_RecursiveDeref(gbm,temp1);


    // temp2 is the midpoint Q add for this action
    temp1 = Cudd_addApply(gbm, addMean, temp2, One);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp2);

    Cudd_RecursiveDeref(gbm,VcurrentMin[i]);
    VcurrentMin[i] = temp1;
    Cudd_Ref(VcurrentMin[i]);
       

    
    temp2=Cudd_addApply(gbm,Cudd_addMaximum,MidPoint_Add,temp1);
    Cudd_Ref(temp2);    
    Cudd_RecursiveDeref(gbm,temp1);

    Cudd_RecursiveDeref(gbm,MidPoint_Add);
    MidPoint_Add = temp2;
    
  }
  
  /* get the approximate policy and output it */
  actionAdd2 = getActionAdd(VcurrentMin, MidPoint_Add, numactions);

  /* Output the stats */

  strcpy(basepath,outpath);
  strcat(basepath,infile);
  sprintf(templog,"-Approx-stats.dat");
  strcat(basepath,templog);
  stats = fopen(basepath,"a");

  double maxVal = (*Cudd_V(Cudd_addFindMax(gbm,MergedMin_Add))).get_max();
  
  dagCount = Cudd_DagSize(actionAdd2);
  dagLeafCount = Cudd_CountLeaves(actionAdd2); 
 
  fprintf(stats,"\n\n\nIterations to convergence %d\n",counter);
  fprintf(stats,"Final execution time: %8.4f  seconds\n\n\n\n\n",temps);
  fprintf(stats,"Memory usage presently: %ld bytes\n",usememory);
  fprintf(stats,"Maximum memory usage: %ld bytes\n\n\n\n",maxusememory);
  fprintf(stats,"Number of nodes in the action DD: %d  internal nodes   %d leaves %d total nodes\n",
	  dagCount-dagLeafCount,dagLeafCount,dagCount);
  fprintf(stats,"Number of nodes in the equivalent tree: %d internal nodes  %d leaves  %d total nodes\n\n",
	  treeCount,treeLeafCount,treeCount+treeLeafCount); 
  fprintf(stats,"Maximum in the Approximate Value Function:%f\n",maxVal);
  //fprintf(stderr,"Maximum in the Approximate Value Function:%f\n",maxVal);
 
  fprintf(dataFile2,"%d %d %d ",dagCount-dagLeafCount,dagLeafCount,dagCount);
 
  // Nodes in the Value ADD
  dagCount = Cudd_DagSize(MergedMin_Add);
  dagLeafCount = Cudd_CountLeaves(MergedMin_Add);
  fprintf(stats,"Number of nodes in the value ADD: %d  internal nodes   %d leaves %d total nodes\n",
	   dagCount-dagLeafCount,dagLeafCount,dagCount);   
  fprintf(dataFile2,"%d %d %d ",
	   dagCount-dagLeafCount,dagLeafCount,dagCount);   
  fprintf(dataFile2,"%d ",counter);
  fprintf(dataFile2,"%10.4f ",temps);
  max_error = get_error(MergedMin_Add);
  double max_extent = get_extent(MergedMin_Add);
  fprintf(dataFile2,"%10.4f %10.4f ", max_error, max_extent);

  Cudd_PrintInfo(gbm,stats);
  fclose(stats);
  ApproximateMidPointValue = MidPoint_Add;
  ApproximatePolicy = actionAdd2;

  return MergedMin_Add;
}
