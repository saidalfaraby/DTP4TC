// bison parser for multi-valued spudd mdps
// written by Jesse Hoey
// November 2002
%{

  //#define USE_POMDP
#ifndef USE_POMDP
  #include "MDP.h"
#else
  #include "POMDP.h"
#endif



  //#define YYSTYPE double
  #define MAXDDS 100000
  //#define REDUCE_PRECISION 1000
  //#define PARSERDEBUG 1
  //#define PARSERDEBUG2 1
  int yylex( void );
  int yyerror( const char *s);
  void error(const char *s, const char *v);
  int curr_primed_ovar, curr_ovar, curr_oval, renormvar;
  DdNode *dds[MAXDDS];
  char *ddnames[MAXDDS];
  DdNode **goodState;
  DdNode **goodStatep;
  double cpt[MAXVALS];
  double pvsum;
  int valindex, level, pvcount, conj, pvindex, numdds;
  int oldconj[64], conjlevel;
  bool doing_reward, doing_dd, doing_obs, doing_action, unnormalized;
  int levelIndex[MAXVARS];
  // tells if a level is primed variable
  bool levelPrime[MAXVARS];
  int branchCount[MAXVARS];
#ifndef USE_POMDP
  extern MDP *__theMDP;
#else
  extern POMDP *__theMDP;
#endif
  map<const char*, int, ltstr> ddindices;
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
%type <dnode> constnode
%type <dnode> primeadd
%type <dnode> currentadd
%type <dnode> theadd
%type <dnode> subadd
%type <dnode> con_dis_add
%type <nme> rootnodename
%token <nme> NAME 
%token <val> REAL
%token <ival> INTGR
%token OPP CLP NAME VARIABLES DISCOUNT TOLERANCE REWARD VALUE ACTION ENDACTION
%token OBSERVATIONS OBSERVE ENDOBSERVE
%token BELIEF ENDBELIEF
%token DISJ   
%token CONJ   
%token RENORM
%token OSB    
%token CSB    
%token HORIZON
%token VAL COST PRIME
%token UNNORM
%token STARTDD ENDDD
%%
input: 
mdp;

mdp:  
varilist obslist ibelief unn actslist reward disc tolhor {
  ddindices.clear();
};

varilist: 
OPP VARIABLES varlist CLP {
  // now we have the variables, so allocate for the NewPrime Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i,j;
  __theMDP->NewPrime = (DdNode ***)malloc(MAXACT*(sizeof(DdNode **)));
  for(i=0;i<MAXACT;i++)
    __theMDP->NewPrime[i] = (DdNode **)malloc(__theMDP->numorigvars*(sizeof(DdNode*)));
  // build the good-state ADDs
  buildGoodStateADDs(__theMDP->orig_vars,__theMDP->vars,__theMDP->prime_vars,__theMDP->numorigvars);
  doing_reward = false;
  doing_action = false;
  doing_obs = false;
  doing_dd = false;
  pvindex = 0;
  pvsum = 0.0;

  numdds = 0;
  // we also want to create a bunch of 'generic' dds like ones for variables that
  // stay the same, which will be named SAME<varname>
  for (i=0; i<__theMDP->numorigvars; i++) {
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** about to parse the dd for %d *****\n",i);
#endif
    ddnames[numdds] = new char[128];
    sprintf(ddnames[numdds],"SAME%s",__theMDP->orig_vars[i].name);
    dds[numdds] = __theMDP->buildSameDD(i);
    Cudd_Ref(dds[numdds]);
    // store in hash table
    ddindices[ddnames[numdds]] = numdds;
    
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** parsed the dd %s at %d*****\n",ddnames[numdds],ddindices[ddnames[numdds]]);
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
#endif
    numdds++;

    // and a bunch of dds - one for each value of each variable
    for (j=0; j<__theMDP->orig_vars[i].nvals; j++) {
      ddnames[numdds] = new char[128];
      sprintf(ddnames[numdds],"%s%s",__theMDP->orig_vars[i].name,__theMDP->orig_vars[i].valname[j]);
      dds[numdds] = __theMDP->buildOneCubeOrig(i,j,true);
      Cudd_Ref(dds[numdds]);
      ddindices[ddnames[numdds]] = numdds;
      numdds++;
    }
  }
};

varlist: 
vardec | varlist vardec;

vardec: 
OPP varname vallist CLP {
  __theMDP->newADDVar();
  __theMDP->numorigvars++;
};

varname: 
NAME { 
  __theMDP->orig_vars[__theMDP->numorigvars].name = strdup($1);
#ifdef PARSERDEBUG
  fprintf(stderr,"name of %d orig_var %s\n",__theMDP->numorigvars,__theMDP->orig_vars[__theMDP->numorigvars].name);
#endif
};

vallist: 
vallist NAME { 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,$2);
} | NAME NAME { 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,$1);
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,$2);
} ;

obslist: /*empty*/ | 
OPP OBSERVATIONS oblist CLP {
#ifdef USE_POMDP
  // now we have the observations, so allocate for the ObsFun Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i;
  __theMDP->ObsFun = (DdNode ***)malloc(MAXACT*(sizeof(DdNode **)));
#ifdef PARSERDEBUG
  fprintf(stderr,"allocating for %d original observations in MDP\n", __theMDP->numorigobs);
#endif
  for(i=0;i<MAXACT;i++)
    __theMDP->ObsFun[i] = (DdNode **)malloc(__theMDP->numorigobs*(sizeof(DdNode*)));
  __theMDP->computeTotalNumObs();
  // also allocate for the beliefSTate
  __theMDP->beliefState = (DdNode **) malloc(__theMDP->numorigvars*sizeof(DdNode *));
  // build the good-state ADDs
  // do we have to do this for observations?
  //buildGoodStateADDs(__theMDP->orig_vars,__theMDP->vars,__theMDP->prime_vars,__theMDP->numorigvars);
#endif  
};
oblist:
obdec | oblist obdec;

obdec:
OPP obsname obsvallist CLP {
#ifdef USE_POMDP
  __theMDP->newADDObs();
  __theMDP->numorigobs++;
#endif
};
obsvallist:
obsvallist NAME {
#ifdef USE_POMDP
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,$2);
#endif
} | NAME NAME { 
#ifdef USE_POMDP
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,$1);
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,$2);
#endif
};
obsname:
NAME {
#ifdef USE_POMDP
  __theMDP->orig_obs[__theMDP->numorigobs].name = strdup($1);
#ifdef PARSERDEBUG
  fprintf(stderr,"name of %d orig_obs %s\n",__theMDP->numorigobs,__theMDP->orig_obs[__theMDP->numorigobs].name);
#endif
#endif
};
ibelief:
/* empty */ |
BELIEF {
  // initialize initBelief
#ifdef USE_POMDP
  __theMDP->initBelief = One;
  Cudd_Ref(__theMDP->initBelief);
#ifdef PARSERDEBUG
  fprintf(stderr,"initialized initBelief to One\n");
#endif
#endif
} belieflist ENDBELIEF {
  pvcount = 0;
};

belieflist:
belieflist belief | belief;


belief:
NAME {
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,$1)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",$1);
    exit(0);
  }
#ifdef PARSERDEBUG
  fprintf(stderr,"found belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
#endif
  level = 0;
} theadd {
#ifdef USE_POMDP
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,$3,4,100);
#endif
#endif
  __theMDP->beliefState[curr_primed_ovar]  = $3;
  Cudd_Ref(__theMDP->beliefState[curr_primed_ovar]);
  DdNode *temp1 =  Cudd_addApply(gbm,Cudd_addTimes,__theMDP->initBelief,$3);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,$3);
  Cudd_RecursiveDeref(gbm,__theMDP->initBelief);
  __theMDP->initBelief = temp1;
#ifdef PARSERDEBUG
  fprintf(stderr,"Initial Belief is currently\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,__theMDP->initBelief,4,100);
#endif
#endif
#endif
};

unn: /* empty */ {
  unnormalized = false;
} | UNNORM {
  unnormalized = true;
};

actslist:
actslist actionordd | actionordd;

actionordd:
action | thedd;

thedd:
ddname {
  conjlevel = 0;
} theadd ENDDD {
  dds[numdds] = $3;
  Cudd_Ref(dds[numdds]);
  Cudd_RecursiveDeref(gbm,$3);
  // store in hash table
  ddindices[ddnames[numdds]] = numdds;
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** parsed the dd %s at %d*****\n",ddnames[numdds],ddindices[ddnames[numdds]]);
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
#endif
  numdds++;
  doing_dd = false;
};

ddname:
STARTDD NAME {
  ddnames[numdds] = strdup($2);
  doing_dd = true;
};

action:
actionname oldactioncost acttreelist observation actioncost ENDACTION {  
  if (pvcount != __theMDP->numorigvars)
    error("missing primed variable in action",__theMDP->actionlist[__theMDP->numactions].name);
  pvcount = 0;
  __theMDP->numactions++;
  doing_action = false;
};

actionname:
ACTION NAME {
  /* name of action is $2, increment action counter */
#ifdef PARSERDEBUG
  fprintf(stderr,"------------- action name %s\n",$2); 
#endif
  __theMDP->actionlist[__theMDP->numactions].name = strdup($2);
  doing_action = true;
};

oldactioncost: /*empty*/ {
  // the action cost is zero
  setActionCost(__theMDP->actionCost,__theMDP->numactions,0.0);
} | REAL {
  setActionCost(__theMDP->actionCost,__theMDP->numactions,$1);
};


actioncost: /* empty */ {
  // the action cost is zero
  //fprintf(stderr,"action cost 0\n");
  setActionCost(__theMDP->actionCost, __theMDP->numactions,0.0);
  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);
} | COST { pvindex = 0; level = 0; doing_reward = true; } theadd {
  // assign the action cost
  __theMDP->actionCost[__theMDP->numactions] = $3;
  Cudd_Ref(__theMDP->actionCost[__theMDP->numactions]);
  Cudd_RecursiveDeref(gbm,$3);

  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);

  // TEMPORARILY REMOVE THIS FOR PRINTING OUT
#ifdef COSTSPRIMED
  removeAllDummysp(__theMDP->actionCost+__theMDP->numactions,__theMDP->numorigvars);
#else
  removeAllDummys(__theMDP->actionCost+__theMDP->numactions,__theMDP->numorigvars);
#endif
  // multiply by -1 (since its a cost)
  Pair *value = new Pair(-1.0);
  DdNode *temp = Cudd_addConst(gbm,value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,__theMDP->actionCost[__theMDP->numactions],temp);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,__theMDP->actionCost[__theMDP->numactions]);
  __theMDP->actionCost[__theMDP->numactions] = temp2;
  doing_reward = false;
};

observation: /*empty*/ | 
OBSERVE {
  doing_obs = true;
} obsfunlist {
  doing_obs = false;
} ENDOBSERVE;

obsfunlist:
obsfunlist obsfun | obsfun;

obsfun:
NAME {
#ifdef USE_POMDP
  // find this variable in orig_obs
  if ((curr_primed_ovar = findOVar(__theMDP->orig_obs,__theMDP->numorigobs,$1)) < 0) {
    fprintf(stderr,"could not find original observation %s\n",$1);
    exit(0);
  }
#endif
} theadd {
#ifdef USE_POMDP
  // observation function in theadd ($3)
  __theMDP->ObsFun[__theMDP->numactions][curr_primed_ovar] = $3;
  Cudd_Ref(__theMDP->ObsFun[__theMDP->numactions][curr_primed_ovar]);
#endif
};
acttreelist: 
acttreelist acttree | acttree;

acttree: 
NAME {
  /* primed variable is $1 theadd is its action tree*/
#ifdef PARSERDEBUG
   fprintf(stderr,"primed var %s\n",$1);
#endif
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,$1)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",$1);
    exit(0);
  }
  level = 0;
  conjlevel = 0;
  // the corresponding binary variables are vars[i] and prime_vars[i]
  // where i:orig_vars[curr_primed_ovar].var1index...orig_vars[curr_primed_ovar].var1index+orig_vars[curr_primed_ovar].nbvars-1
} theadd {
  // the new prime diagram is that returned in theadd
  // its $3 because the above action is also counted
  removeDummyp(&($3),curr_primed_ovar);
  DdNode *tmp2;
  DdNode *temp = sumOutPrime($3,__theMDP->orig_vars+curr_primed_ovar,__theMDP->prime_vars);
  Cudd_Ref(temp);


  if (unnormalized) {
    // renormalize - this can add in Dummy states!
#ifdef PARSERDEBUG2
    fprintf(stderr,"dividing this:\n");
    pdd($3);
    fprintf(stderr,"by this (previous summed over %d):\n",curr_primed_ovar);
    pdd(temp);
#endif    
    tmp2 = Cudd_addApply(gbm,Cudd_addDivide,$3,temp);

    Cudd_Ref(tmp2);

    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,$3);
    temp = tmp2;
    // remove dummy states
    //removeAllDummys(&temp,__theMDP->numorigvars);
    //removeAllDummysp(&temp,__theMDP->numorigvars);


#ifdef REDUCE_PRECISION
    // reduce the precision of the final diagram
    Pair *prec = new Pair();
    (*prec).set(REDUCE_PRECISION);
    DdNode *predd = Cudd_addConst(gbm,prec);
    Cudd_Ref(predd);
    tmp2 = Cudd_addApply(gbm,reducePrecision,temp,predd);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,predd);
    temp = tmp2;
    // and normalize again
    tmp2 = sumOutPrime(temp,__theMDP->orig_vars+curr_primed_ovar,__theMDP->prime_vars);
    Cudd_Ref(tmp2);
    tmp1 = Cudd_addApply(gbm,Cudd_addDivide,temp,tmp2);
    Cudd_Ref(tmp1);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,tmp2);
    temp = tmp1;
#endif



    __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = temp;
    Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
    
    

  } else {
    // check if normalized
    // if normalized - temp should be exactly 1 only for non-dummy values
    // so, remove dummys, take 1-temp and then remove the dummys again - 
    // the result should be Zero
    removeAllDummys(&temp,__theMDP->numorigvars);
    removeAllDummysp(&temp,__theMDP->numorigvars);
    tmp2 = Cudd_addApply(gbm,Cudd_addMinus,One,temp);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,temp);
    temp = tmp2;
    removeAllDummys(&temp,__theMDP->numorigvars);
    removeAllDummysp(&temp,__theMDP->numorigvars);
    if (temp != Zero) {
      fprintf(stderr,"CPT for primed variable %s is not normalized:\n",$1);
      Cudd_PrintMinterm(gbm,$3);
      fprintf(stderr,"curr_primed_ovar: %d\n summed diagram:\n",curr_primed_ovar);
      Cudd_PrintDebug(gbm,temp,4,100);
      exit(-1);
    } else {
      __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = $3;
      Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
      Cudd_RecursiveDeref(gbm,$3);
    }
  } 
  Cudd_RecursiveDeref(gbm,temp);
  // removeDummy(&(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]),curr_primed_ovar);


#ifdef PARSERDEBUG
  fprintf(stderr,"the add for %d action %d ovar is \n",__theMDP->numactions, curr_primed_ovar);
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar],4,100);
#endif
#endif
  pvcount++;
};



theadd: 
OSB CONJ {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 1;
  conjlevel++;
} con_dis_add CSB {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);

#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
} | 
OSB DISJ {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 0;
  conjlevel++;
} con_dis_add CSB {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);
#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
} |
OSB RENORM NAME PRIME {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 2;
  conjlevel++;
} theadd CSB {
  //renormalize $6 over variable $3
  if ((renormvar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,$3)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",$3);
    exit(0);
  }
  DdNode *temp = NULL;
  __theMDP->normalizeFunction($6,&temp,renormvar);
  $$ = temp;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,$6);
#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
} |
OPP rootnodename {
  /* current root node is $2 */
#ifdef PARSERDEBUG
  fprintf(stderr,"********************************  root node %s\n",$2);
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,$2)) < 0) 
    error("could not find variable",$2);
  levelIndex[level] = curr_ovar;
  branchCount[level] = 0;
} 
currentadd CLP {
  $$ = $4;
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$4);
#ifdef PARSERDEBUG
  fprintf(stderr,"***************************** the add is\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,$$,4,100);
#endif
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,$2)) < 0) 
    error("could not find variable",$2);
  if (branchCount[level] != __theMDP->orig_vars[curr_ovar].nvals) {
    fprintf(stderr,"%d branches missing from variable %s primed variable %s in action %s\n",
	    __theMDP->orig_vars[curr_ovar].nvals-branchCount[level],
	    __theMDP->orig_vars[curr_ovar].name,__theMDP->orig_vars[curr_primed_ovar].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
} | 
constnode |
OPP primeadd CLP {

  // check to make sure it matches the primed variables number
  if (!doing_reward) {
    if (doing_obs) {
#ifdef USE_POMDP
      if (pvindex != __theMDP->orig_obs[curr_primed_ovar].nvals) 
	error("Some cpt values missing from observation",__theMDP->orig_obs[curr_primed_ovar].name);
#endif
    } else {
      if (pvindex != __theMDP->orig_vars[curr_primed_ovar].nvals) 
	error("Some cpt values missing from primed variable",__theMDP->orig_vars[curr_primed_ovar].name);
    }
    $$ = $2;
  } else if (doing_reward && pvindex != 1) {
    error("only one reward value per state","please!");
  } else {
    $$ = $2;
  }
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$2);

  // reset for next time
  pvindex = 0;
  pvsum = 0.0;
};

rootnodename:
NAME {
#ifdef PARSERDEBUG
  fprintf(stderr,"parsing root node name %s\n",$1);
#endif
  $$ = $1;
  levelPrime[level] = false;
}
| NAME PRIME {
#ifdef PARSERDEBUG
  fprintf(stderr,"saw a prime on variable %s\n",$1);
#endif
  $$ = $1;
  levelPrime[level] = true;
};

con_dis_add:
con_dis_add theadd {
  int i;
  DdNode *temp;
  if (conj == 1) {
#ifdef PARSERDEBUG
    fprintf(stderr,"multiplying these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addTimes,$1,$2);
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"adding these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  }
#ifdef PARSERDEBUG
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,$1,1,2);
  Cudd_PrintDebug(gbm,$2,1,2);
#endif
#endif
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
#ifdef PARSERDEBUG
  fprintf(stderr,"branch %s level increasing to %d\n",$2,level+1);
#endif
  level++;
} theadd CLP {
  level--;
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed branch %s\n",$2); 
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,$4,4,100);
#endif
#endif
  if ((curr_oval = findOVal(__theMDP->orig_vars, levelIndex[level],$2)) < 0) {
    fprintf(stderr,"could not find value %s in variable %s action %s\n",
	    $2,__theMDP->orig_vars[levelIndex[level]].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
  // need to check if curr_ovar is a primed variable here
  if (!levelPrime[level] && !doing_obs)
    $$ = buildCubeCPT(__theMDP->vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],$4);
  else
    $$ = buildCubeCPT(__theMDP->prime_vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],$4);
  Cudd_RecursiveDeref(gbm,$4);
#ifdef PARSERDEBUG
  fprintf(stderr,"branch CPT is\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,$$,4,100);
#endif
#endif
  branchCount[level]++;
};

// old style primeadd - just a list of probabilities
// must be at least two!
primeadd:
primeadd constadd {
  $$ = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
  Cudd_RecursiveDeref(gbm,$2);
}  |
constadd constadd {
  $$ = Cudd_addApply(gbm,Cudd_addPlus,$1,$2);
  Cudd_Ref($$);
  Cudd_RecursiveDeref(gbm,$1);
  Cudd_RecursiveDeref(gbm,$2);
};

constnode:
OPP REAL CLP {
  Pair *value = new Pair();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  (*value).set(((double) ((int) ($2*REDUCE_PRECISION)))/(1.0*REDUCE_PRECISION));
#else
  (*value).set($2);
#endif
  //fprintf(stderr,"setting %f to %f\n",$2,((double) ((int) ($2*100)))/100.0);
  $$ = Cudd_addConst(gbm,value);
  Cudd_Ref($$);
} | 
OPP NAME CLP {
  // search for NAME in list of defined DDs
  // by looking in hash table
#ifdef PARSERDEBUG
  fprintf(stderr,"********************** looking for dd %s\n",$2);
#endif
  map<const char*, int, ltstr>::iterator ddit = ddindices.find($2);
  if (ddit==ddindices.end()) {
    fprintf(stderr,"could not find defined DD %s in list\n",$2);
    exit(-1);
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** Found defined DD %s in list at %d\n",$2,(*ddit).second);
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[(*ddit).second],1,2);
#endif
#endif
    $$ = dds[(*ddit).second];
    Cudd_Ref($$);
  }
};

constadd:
REAL {

  Pair *value = new Pair();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  (*value).set(((double) ((int) ($1*REDUCE_PRECISION)))/(1.0*REDUCE_PRECISION));
#else
  (*value).set($1);
#endif
  //fprintf(stderr,"setting %f to %f\n",$1,((double) ((int) ($1*100)))/100.0);


  DdNode *temp = Cudd_addConst(gbm,value);
  delete value;
  Cudd_Ref(temp);
  if (!doing_reward) {
#ifdef PARSERDEBUG
    fprintf(stderr,"found a constant %f at variable %s\n",$1,__theMDP->orig_vars[curr_primed_ovar].name);
#endif
    // buildCubeCPT does its own reffing
    if (doing_obs) {
#ifdef USE_POMDP
      $$ = buildCubeCPT(__theMDP->obs,pvindex,__theMDP->orig_obs[curr_primed_ovar],temp);
#endif
    } else {
      $$ = buildCubeCPT(__theMDP->prime_vars,pvindex,__theMDP->orig_vars[curr_primed_ovar],temp);
    }
    Cudd_RecursiveDeref(gbm,temp);
  } else {
    //fprintf(stderr," found a constant %f pvindex %d\n",$1,pvindex);
    $$ = temp;
  }    
  pvindex++;
  pvsum += $1;
};
  

reward:
REWARD { 
#ifdef PARSERDEBUG
  fprintf(stderr, "doing reward\n");
#endif
  pvindex = 0; level = 0; doing_reward = true; 
} theadd {
  __theMDP->RewardD = $3;
  Cudd_Ref(__theMDP->RewardD);
  Cudd_RecursiveDeref(gbm,$3);

  __theMDP->RewardDNoDummy =   __theMDP->RewardD; 
  Cudd_Ref(__theMDP->RewardDNoDummy);

  removeAllDummys(&(__theMDP->RewardD),__theMDP->numorigvars);
};

disc: 
DISCOUNT  REAL { 
  Pair discountPair;
  discountPair.set($2);
  __theMDP->discount = Cudd_addConst(gbm,&discountPair);
  Cudd_Ref(__theMDP->discount);
} | DISCOUNT {
  doing_reward = true;
} theadd {
  __theMDP->discount = $3;
  Cudd_Ref(__theMDP->discount);
};

tolhor:
tol | hor;

tol: 
TOLERANCE REAL { 
  tolerance = $2; 
  __theMDP->horizon = -1.0;
};

hor:
HORIZON REAL {
  tolerance = 0.0;
  __theMDP->horizon = int($2);
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
// same for primed variables
void removeDummyp(DdNode **np, int cpovar) {
   DdNode *tmp,*tmp2;
  tmp = Cudd_addApply(gbm,Cudd_addTimes,goodStatep[cpovar],np[0]);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,np[0]);
  np[0] = tmp;
}
void removeAllDummys(DdNode **np, int nv) {
  int v;
  for (v=0; v<nv; v++) 
    removeDummy(np,v);
}
void removeAllDummysp(DdNode **np, int nv) {
  int v;
  for (v=0; v<nv; v++) 
    removeDummyp(np,v);
}
// add a 1 to np for each dummy state
void addDummyStates(DdNode **np, int nv) {
  int v;
  DdNode *tmp,*tmp2;
  for (v=0; v<nv; v++) {
    tmp = Cudd_addApply(gbm,Cudd_addMinus,One,goodState[v]);
    Cudd_Ref(tmp);
    tmp2 = Cudd_addApply(gbm,Cudd_addPlus,np[0],tmp);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,tmp);
    Cudd_RecursiveDeref(gbm,np[0]);
    np[0] = tmp2;
  }
    
}
// function to build the goodState ADDs once at the start
void buildGoodStateADDs(onum *ov, rnum *vars, rnum *pvars, int nov) {
  int a,v,va,pv,i,j;
  goodState = new DdNode*[nov];
  goodStatep = new DdNode*[nov];
  DdNode *tmp,*tmp2;
  // construct the goodState ADDs for each varable
  for (v=0; v < nov; v++) {
    goodState[v] = Zero;
    Cudd_Ref(goodState[v]);
    for (va=0; va< ov[v].nvals; va++) {
      tmp = buildCubeCPT(vars,va,ov[v],One);
      tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,goodState[v]);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,goodState[v]);
      goodState[v] = tmp2;
    }
    // build prime good state ADDs too
    goodStatep[v] = Zero;
    Cudd_Ref(goodStatep[v]);
    for (va=0; va< ov[v].nvals; va++) {
      tmp = buildCubeCPT(pvars,va,ov[v],One);
      tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,goodStatep[v]);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,goodStatep[v]);
      goodStatep[v] = tmp2;
    }
    //fprintf(stderr,"goodState for %d is \n",v);
#ifdef PARSERDEBUG2
    //Cudd_PrintDebug(gbm,goodState[v],4,100);
#endif

  }
}

void setActionCost(DdNode **ac, int nac, double val) {
  //  fprintf(stderr,"Action cost is %f for action %d\n",val,nac);
  Pair *value = new Pair();
  (*value).set((double)(-1.0*fabs(val)));
  ac[nac] = Cudd_addConst(gbm,value);
  Cudd_Ref(ac[nac]);
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
#ifdef PARSERDEBUG2
  //Cudd_PrintDebug(gbm,branch,4,100);
#endif

  DdNode *cubecpt = Cudd_addApply(gbm,Cudd_addTimes,branch,add);
  Cudd_Ref(cubecpt);
  Cudd_RecursiveDeref(gbm,branch);

  //fprintf(stderr,"built cube:\n");
#ifdef PARSERDEBUG2
  //Cudd_PrintDebug(gbm,cubecpt,4,100);
#endif

  delete [] phase;
  delete [] arrayofvars;
  return cubecpt;
}


int findOVal(onum *ov, int ovar, const char *findval) {
  int oval = 0;
  while (oval < ov[ovar].nvals && strcmp(ov[ovar].valname[oval],findval) != 0)
    oval++;
  if (oval < ov[ovar].nvals)
    return oval;
  else 
    return -1;
}


// returns the index of orignal variable findname, or -1 if not found
int findOVar(onum *ov, int nov, const char *findname) {
  int ovar = 0;
  while (ovar < nov && strcmp(ov[ovar].name,findname) != 0)
    ovar++;
  if (ovar < nov)
    return ovar;
  else 
    return -1;
}
// adds a new multi-valued variable to a list of such variables
void addVar(onum *ovar, int vnum, char * varn) {
  ovar[vnum].valname[ovar[vnum].nvals] = strdup(varn);
#ifdef PARSERDEBUG
  fprintf(stderr,"adding %d value (%s) to var %d\n",ovar[vnum].nvals,ovar[vnum].valname[ovar[vnum].nvals],vnum);
#endif
  ovar[vnum].nvals++;
}
//renormalizes a cube of primed variables
// result is not reffed
// does not deref cube
DdNode * renormalizeCube(DdNode *cube, double pvsum) {
  // renormalize the sum - divide by pvsum
  Pair *value = new Pair(1.0/pvsum);
  DdNode *temp = Cudd_addConst(gbm,value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,cube,temp);
  Cudd_RecursiveDeref(gbm,temp);
  return temp2;
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

