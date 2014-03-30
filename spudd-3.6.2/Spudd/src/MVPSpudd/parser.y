// bison parser for multi-valued spudd mdps
// written by Jesse Hoey
// November 2002
%{
#include "pspudd.hh"
  //#define YYSTYPE double
  #define YYDEBUG 1
  int yylex( void );
  int yyerror( const char *s);
  void error(const char *s, const char *v);
  int curr_primed_ovar, curr_ovar, curr_oval;
  DdNode *subADDlist[MAXVARS][MAXVALS];
  DdNode **goodState;
  double cpt[MAXVALS];
  double pvsum;
  int valindex, level, pvcount, conj, pvindex;
  bool doing_reward;
  int levelIndex[MAXVARS];
  int branchCount[MAXVARS];
%}
%union {
  double val;
  int ival;
  DdNode * dnode;
  DdNode ** p_dnode;
  double *darray;
  char *nme;
}
%type <dnode> constadd
%type <dnode> primeadd
%type <dnode> currentadd
%type <dnode> theadd
%type <dnode> subadd
%type <dnode> con_dis_add
%token <nme> NAME 
%token <val> REAL
%token <ival> INTGR
%token OPP CLP NAME VARIABLES DISCOUNT TOLERANCE REWARD VALUE ACTION ENDACTION
%token DISJ   
%token CONJ   
%token OSB    
%token CSB    
%token HORIZON
%token VAL    

%%
input: 
mdp;

mdp:  
varilist actslist reward disc tolhor;

varilist: 
OPP VARIABLES varlist CLP {
  // now we have the variables, so allocate for the NewPrime Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i;
  NewPrime = (DdNode ***)malloc(MAXACT*(sizeof(DdNode **)));
  for(i=0;i<MAXACT;i++)
    NewPrime[i] = (DdNode **)malloc(numorigvars*(sizeof(DdNode*)));
  // build the good-state ADDs
  buildGoodStateADDs();
  doing_reward = false;
  pvindex = 0;
  pvsum = 0.0;
};

varlist: 
vardec | varlist vardec;

vardec: 
OPP varname vallist CLP {
  newADDVar(orig_vars,numorigvars);
  numorigvars++;
};

varname: 
NAME { 
  orig_vars[numorigvars].name = strdup($1);
};

vallist: 
vallist NAME { 
  addVar(orig_vars,numorigvars,$2);
} | NAME NAME { 
  addVar(orig_vars,numorigvars,$1);
  addVar(orig_vars,numorigvars,$2);
} ;

actslist:
actslist action | action;

action:
simpleaction {
  pvcount = 0;
  numactions++;
};

simpleaction:
actionname actioncost acttreelist ENDACTION {  
  if (pvcount != numorigvars)
    error("missing primed variable in action",actionlist[numactions].name);
};

actionname:
ACTION NAME {
  /* name of action is $2, increment action counter */
  //  fprintf(stderr,"action name %s\n",$2); 
  actionlist[numactions].name = strdup($2);
};

actioncost: /*empty*/ {
  // the action cost is zero
  setActionCost(numactions,0.0);
} | REAL {
  setActionCost(numactions,$1);
};

acttreelist: 
acttreelist acttree | acttree;

acttree: 
NAME {
  /* primed variable is $1 theadd is its action tree*/
  //fprintf(stderr,"primed var %s\n",$1);
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar($1)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",$1);
    exit(0);
  }
  level = 0;
  // the corresponding binary variables are vars[i] and prime_vars[i]
  // where i:orig_vars[curr_primed_ovar].var1index...orig_vars[curr_primed_ovar].var1index+orig_vars[curr_primed_ovar].nbvars-1
} theadd {
  // the new prime diagram is that returned in theadd
  // its $3 because the above action is also counted
  NewPrime[numactions][curr_primed_ovar] = $3;
  Cudd_Ref(NewPrime[numactions][curr_primed_ovar]);
  Cudd_RecursiveDeref(gbm,$3);
  removeDummy(&(NewPrime[numactions][curr_primed_ovar]),curr_primed_ovar);
  //fprintf(stderr,"NewPrime[%d][%d] is :\n",numactions,curr_primed_ovar);
  //Cudd_PrintDebug(gbm,NewPrime[numactions][curr_primed_ovar],4,100);
  pvcount++;
};


theadd: 
OSB CONJ {
  conj = 1;
} con_dis_add CSB {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);
} | 
OSB DISJ {
  conj = 0;
} con_dis_add CSB {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);
} |
OPP NAME {
  /* current root node is $2 */
  //fprintf(stderr,"root node %s\n",$2);
  if ((curr_ovar = findOVar($2)) < 0) 
    error("could not find variable",$2);
  levelIndex[level] = curr_ovar;
  branchCount[level] = 0;
}
currentadd CLP {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);
  //fprintf(stderr,"add built for root %s\n",$2);
  //Cudd_PrintDebug(gbm,$$,4,100);
  if ((curr_ovar = findOVar($2)) < 0) 
    error("could not find variable",$2);
  if (branchCount[level] != orig_vars[curr_ovar].nvals) {
    fprintf(stderr,"%d branches missing from variable %s primed variable %s in action %s\n",
	    orig_vars[curr_ovar].nvals-branchCount[level],
	    orig_vars[curr_ovar].name,orig_vars[curr_primed_ovar].name,actionlist[numactions].name);
    exit(-1);
  }
} | 
OPP primeadd CLP {
  $$ = $2;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$2);

  // check to make sure it matches the primed variables number
  if (!doing_reward) {
    if (pvindex != orig_vars[curr_primed_ovar].nvals) 
      error("Some cpt values missing from primed variable",orig_vars[curr_primed_ovar].name);
    if (fabs(pvsum - 1.0) > 1e-5) 
      error("Probabilities don't sum to one for primed variable",orig_vars[curr_primed_ovar].name);
  } else if (doing_reward && pvindex != 1) {
    error("only one reward value per state","please!");
  }
  // reset for next time
  pvindex = 0;
  pvsum = 0.0;
};

con_dis_add:
con_dis_add theadd {
  int i;
  DdNode *temp;
  if (conj == 1)
    temp = Cudd_addApply(gbm,Cudd_addTimes,$1,$2);
  else
    temp = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,$1);
  Cudd_RecursiveDeref(gbm,$2);
  $$ = temp;
}| 
theadd {
  $$ = $1;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
};

currentadd:
currentadd subadd {
  $$ = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
  Cudd_RecursiveDeref(gbm,$2);
} | 
subadd {
  $$ = $1;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
};

// subadd is the add rooted at the curr_oval branch of the
// current original root variable. 
// semantic value of subadd is that (recursively parsed) add
// rooted at the cube of the corresponding binary variables to the
// value of that branch
subadd: 
OPP NAME {
  //fprintf(stderr,"branch %s level increasing to %d\n",$2,level+1);
  level++;
} theadd CLP {
  level--;
  //fprintf(stderr,"parsed branch %s\n",$2); 
  if ((curr_oval = findOVal(levelIndex[level],$2)) < 0) {
    fprintf(stderr,"could not find value %s in variable %s action %s\n",
	    $2,orig_vars[levelIndex[level]].name,actionlist[numactions].name);
    exit(-1);
  }
  $$ = buildCubeCPT(vars,curr_oval,orig_vars[levelIndex[level]],$4);
  Cudd_RecursiveDeref(gbm,$4);
  branchCount[level]++;
};

primeadd:
primeadd constadd {
  $$ = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
  Cudd_RecursiveDeref(gbm,$2);
}  |
constadd {
  $$ = $1;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
};

constadd:
REAL {
  Pair *value = new Pair();
  (*value).set($1);
  DdNode *temp = Cudd_addConst(gbm,value);
  Cudd_Ref(temp);
  if (!doing_reward) {
    //fprintf(stderr,"found a constant %f at pvindex %d variable %s\n",$1,pvindex,orig_vars[curr_primed_ovar].name);
    // buildCubeCPT does its own reffing
    $$ = buildCubeCPT(prime_vars,pvindex,orig_vars[curr_primed_ovar],temp);
    Cudd_RecursiveDeref(gbm,temp);
  } else {
    $$ = temp;
  }    
  pvindex++;
  pvsum += $1;
};
  

reward:
REWARD { /*fprintf(stderr, "doing reward\n"); */pvindex = 0; level = 0; doing_reward = true; } theadd {
  RewardD = $3;
  Cudd_Ref(RewardD);
  Cudd_RecursiveDeref(gbm,$3);
  removeAllDummys(&RewardD);
};


disc: 
DISCOUNT  REAL { 
  discount_factor = $2; 
};
tolhor:
tol | hor;

tol: 
TOLERANCE REAL { 
  tolerance = $2; 
  horizon = -1.0;
};

hor:
HORIZON REAL {
  tolerance = 0.0;
  horizon = int($2);
};

%%

#include "lex.yy.c"
// function to add in zeros for all the dummy variables
// do this by mutliplying each NewPrime diagram by 
// a 'goodstate' ADD for each variable. The goodstate ADD
// for an original (mv) variable is the sum of all cubes over
// binary variables corresponding to valid original values 
// and 0 to all others (which correspond to dummy values)
void removeDummy(DdNode **np, int cpovar) {
   DdNode *tmp,*tmp2;
  tmp = Cudd_addApply(gbm,Cudd_addTimes,goodState[cpovar],np[0]);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,np[0]);
  np[0] = tmp;
}
void removeAllDummys(DdNode **np) {
  int v;
  for (v=0; v<numorigvars; v++) 
    removeDummy(np,v);
}
// function to build the goodState ADDs once at the start
void buildGoodStateADDs() {
  int a,v,va,pv,i,j;
  goodState = new DdNode*[numorigvars];
  DdNode *tmp,*tmp2;
  // construct the goodState ADDs for each varable
  for (v=0; v < numorigvars; v++) {
    goodState[v] = Zero;
    Cudd_Ref(goodState[v]);
    for (va=0; va<orig_vars[v].nvals; va++) {
      tmp = buildCubeCPT(vars,va,orig_vars[v],One);
      tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,goodState[v]);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,goodState[v]);
      goodState[v] = tmp2;
    }
    //fprintf(stderr,"goodState for %d is \n",v);
    //Cudd_PrintDebug(gbm,goodState[v],4,100);
  }
}

void setActionCost(int nac, double val) {
  //  fprintf(stderr,"Action cost is %f for action %d\n",val,nac);
  Pair *value = new Pair();
  (*value).set((double)(-1.0*fabs(val)));
  actionCost[nac] = Cudd_addConst(gbm,value);
  Cudd_Ref(actionCost[nac]);
}
// builds an add  rooted at a cube over the 
// variables in v corresponding to ovar's coval branch
// with add at that leaf
DdNode * buildCubeCPT(rnum *v, int coval, onum ovar, DdNode *add) {
  // the returned  is a DdNode* referring to the add at this branch
  //build the cube for this branch
  int nbv = ovar.nbvars;
  int nbvals = int(pow(2.0,nbv));
  int *phase = new int[nbv];
  int i;
  for (i=0; i<nbv; i++)
    phase[i] = 0;
  //phase should be nbvals-coval in binary
  int tmp = nbvals-coval-1;
  i=nbv-1;

  while (tmp > 0 && i >=0) {
    phase[i--] = tmp%2;
    tmp = tmp/2;
  }
  /*
  fprintf(stderr,"variable %s coval %d nbv: %d phase: ",ovar.name,tmp,nbv);
  for (i=0; i<nbv; i++)
    fprintf(stderr,"%d ",phase[i]);
  fprintf(stderr,"\n");
  */
  DdNode **arrayofvars = new DdNode*[nbv];
  for (i=0; i<nbv; i++) {
    arrayofvars[i] = v[ovar.var1index+i].add_var;
    Cudd_Ref(arrayofvars[i]);
  }
  DdNode *branch = Cudd_addComputeCube(gbm,arrayofvars,phase,nbv);
  Cudd_Ref(branch);
  for (i=0; i<nbv; i++) 
    Cudd_RecursiveDeref(gbm,arrayofvars[i]);

  //fprintf(stderr,"built branch:\n");
  //Cudd_PrintDebug(gbm,branch,4,100);

  DdNode *cubecpt = Cudd_addApply(gbm,Cudd_addTimes,branch,add);
  Cudd_Ref(cubecpt);
  Cudd_RecursiveDeref(gbm,branch);

  //fprintf(stderr,"built cube:\n");
  //Cudd_PrintDebug(gbm,cubecpt,4,100);

  delete [] phase;
  delete [] arrayofvars;
  return cubecpt;
}


int findOVal(int ovar, const char *findval) {
  int oval = 0;
  while (oval < orig_vars[ovar].nvals && strcmp(orig_vars[ovar].valname[oval],findval) != 0)
    oval++;
  if (oval < orig_vars[ovar].nvals)
    return oval;
  else 
    return -1;
}
// returns the index of orignal variable findname, or -1 if not found
int findOVar(const char *findname) {
  int ovar = 0;
  while (ovar < numorigvars && strcmp(orig_vars[ovar].name,findname) != 0)
    ovar++;
  if (ovar < numorigvars)
    return ovar;
  else 
    return -1;
}
// adds a new multi-valued variable to a list of such variables
void addVar(onum *ovar, int vnum, char * varn) {
  ovar[vnum].valname[ovar[vnum].nvals] = strdup(varn);
  ovar[vnum].nvals++;
}

// constructs a new ADD variable for the novar original multi-valued variable
void newADDVar(onum *ovar, int novar) {
  char tmp[MAXLEN];
  
  // figure out how many new variables to add
  // I realise this could be done by taking logs and ceil, but this is safer
  int i,bvars(1),bnvals(2);
  while (bnvals < ovar[novar].nvals) {
    bvars++;
    bnvals = bnvals*2;
  }
  ovar[novar].nbvars = bvars;
  ovar[novar].var1index = numvars;
  for (i=0; i<bvars; i++) {
    vars[numvars].orig_num = novar;
    vars[numvars].number = numvars;
    prime_vars[numvars].number = numvars;
    sprintf(tmp,"%s%d",ovar[novar].name,i);
    vars[numvars].name = strdup(tmp);
    sprintf(tmp,"%sr",vars[numvars].name);
    prime_vars[numvars].name = strdup(tmp);
    prime_vars[numvars].add_var = Cudd_addNewVar(gbm);
    Cudd_Ref(prime_vars[numvars].add_var);
    vars[numvars].add_var = Cudd_addNewVar(gbm);
    Cudd_Ref(vars[numvars].add_var);
    numvars++;
  }
}

int yyerror (const char *s)  /* Called by yyparse on error */
{
  fprintf (stderr,"%s\n", s);
  exit(-1);
}
void error(const char *s, const char *v) {
  fprintf(stderr, "%s %s\n",s,v);
  exit(-1);
}

