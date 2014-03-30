/***********************************************************************
   parser for boolean style data
   converted for use with mdp class
   this just requires having a global pointer to the mdp object being
   parsed. 
   
   NOTE: this parser grew directly from rev. 2.0 of parser.c in Spudd/src/MVPSpudd
  original parser.c code by: Jesse Hoey & Robert St-Aubin
  these revisions by : Jesse Hoey April 2003

***********************************************************************/


#include<string.h>
#include<stdlib.h>
#include<stdio.h>

#include "MDP.h"
#include "mdpparser.tab.h"
extern MDP *__theMDP;

int parseVars(FILE *errf, int *numvars, rnum *vars, rnum *prime_vars)
{
  int token,index=0;
  char temptext[MAXLEN];

  if (yylex() == OPP) {
    while ((token = yylex()) == NAME) {
      strcpy(temptext,yytext);
      vars[index].number = index;
      prime_vars[index].number = index;
      vars[index].name = strdup(temptext);
      prime_vars[index].name = strdup(strcat(temptext,"r"));
      index++;
      (*numvars)++;
    } 
    if (token == CLP) 
      return 1;
    else 
      fprintf(errf,"expected ) at end of variable list\n");
  } else {
    fprintf(errf,"expected variable list of the form\n");
    fprintf(errf,"variables (a b c)\n");
  }
  return 0;
}

int findVar(int *index, int numvars, char *vname, rnum *vars)
{
  int i;
  for (i=0; i<numvars; i++) {
    if (strcmp(vars[i].name,vname) == 0) {
      *index = i;
      return 1;
    }
  }
  return 0;
}

DdNode *parseADD(FILE *errf,int numvars, rnum *vars, int firstGone)
{
  DdNode *theADD,*leftADD,*rightADD;
  char temptext[MAXLEN];
  int vindex,token;
  Pair *value = new Pair();

  if (firstGone || yylex() == OPP) {
    if ((token = yylex()) == REAL) {
      (*value).set((double) 1.0*atof(yytext));
      theADD = Cudd_addConst(gbm,value);
      Cudd_Ref(theADD);
    } else if (token == NAME) {
      strcpy(temptext,yytext);
      if (findVar(&vindex, numvars, temptext, vars)) {
	if ((leftADD = parseADD(errf, numvars, vars,0)) &&
	    (rightADD = parseADD(errf, numvars, vars,0))) {
	  theADD = Cudd_addIte(gbm,vars[vindex].add_var,leftADD,rightADD);
	  // This was the problem discovered by Carlos Guestrin's MDP - 
	  // these should be Derefs, not recursive
	  Cudd_Ref(theADD);
	  Cudd_RecursiveDeref(gbm,leftADD);
	  Cudd_RecursiveDeref(gbm,rightADD);
	} else {
	  return NULL;
	}
      } else {
	fprintf(errf,"variable %s not found in list of variables\n",temptext);
	return NULL;
      }
    }
    if (yylex() == CLP) {
      return theADD;
    } else {
      return NULL;
    }
  } else {
    fprintf(errf,"Error in specified tree structure\n");
    return NULL;
  }
}

DdNode *parseDiagram(FILE *errf, int numvars, rnum *vars, int firstTok) {
  
  int token3,conj = 1;
  DdNode *cumAct, *temp, *tempAct;
  char temptext[MAXLEN];
  
  if (firstTok == OSB || (firstTok != OPP && ((token3 = yylex()) == OSB)))
    if ((token3 = yylex()) == CONJ)
      conj = 1;
    else if (token3 == DISJ) 
      conj = 0;
    else {
      fprintf(errf,"Only conjunction (*) or disjusnction (+) allowed\n");
      return NULL;
    }
  else if (firstTok == OPP || token3 == OPP) {
    if ((tempAct = parseADD(errf, numvars,vars,1))) 
      return tempAct;
    else 
      return NULL;
  }
  /* now must parse a set of diagrams */
  if (conj == 1) 
    cumAct = Cudd_ReadOne(gbm);
  else
    cumAct = Cudd_ReadZero(gbm);
  Cudd_Ref(cumAct);

  while ((token3 = yylex())  != CSB) {
    tempAct = parseDiagram(errf,numvars,vars,token3);
    if (conj == 1)
      temp = Cudd_addApply(gbm,Cudd_addTimes,cumAct,tempAct);
    else
      temp = Cudd_addApply(gbm,Cudd_addPlus,cumAct,tempAct);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,cumAct);
    Cudd_RecursiveDeref(gbm,tempAct);
    cumAct = temp;
  }
  if (token3 == CSB)
    return cumAct;
  else
    return NULL;
}
int parseActions(FILE *errf, int numvars, int *numactions, rnum *vars, actt *actionlist)
{
  int tokenac,token1,token2,token3,aindex=0,vindex;
  DdNode *tempAct;
  char temptext[MAXLEN];
  Pair *value = new Pair();

  while ((token1 = yylex()) == ACTION) {
    if ((token2 = yylex()) == NAME) {
      strcpy(temptext,yytext);
      actionlist[aindex].name = strdup(temptext);
      
      // action costs
      if ((token3 = yylex()) == REAL || token3 == INTGR) {
	(*value).set((double)(-1.0*fabs(atof(yytext))));
	__theMDP->actionCost[aindex] = Cudd_addConst(gbm,value);
	token3 = yylex();
      }  else {
	//(*value).set(0.0);
	//	actionCost[aindex] = Cudd_addConst(gbm,value);
	__theMDP->actionCost[aindex] = Cudd_ReadZero(gbm);
      }
      Cudd_Ref(__theMDP->actionCost[aindex]);
      while (token3 == NAME) {
	strcpy(temptext,yytext);
	if (findVar(&vindex, numvars, temptext, vars)) {
	  if (tempAct = parseDiagram(errf,numvars,vars,0)) {
	    __theMDP->Actions[aindex][vindex] = tempAct;
	    Cudd_Ref(__theMDP->Actions[aindex][vindex]);
	    Cudd_RecursiveDeref(gbm, tempAct);
	  } else {
	    fprintf(errf,"failed parse of a probability diagram\n");
	    return 0;
	  }
	} else {
	  fprintf(errf,"action %d, variable %s was not in the variable list\n",aindex,temptext);
	  return 0;
	}
	token3 = yylex();
      }
      if (token3 == ENDACTION) {
	aindex++;
      } else {
	fprintf(errf,"expected endaction keyword at end of action %d\n",aindex);
	return 0;
      }
    } else {
      fprintf(errf,"expected action name after action %d\n",aindex);
      return 0;
    }
  }
  if (token1 == REWARD) {
    *numactions = aindex;
    return 1;
  } else {
    fprintf(errf,"expected reward keyword\n");
    return 0;
  }
}

int pparse(FILE *errf, double *discount, double *horizon, double *tolerance,  int *numvars, int *numactions, rnum *vars, rnum *primevars, actt *actionlist)
{
  int i,k,token;
  DdNode *temprew;

  if(yylex() == VARIABLES) {
    if (parseVars(errf,numvars, vars, primevars)) { 

      /* Create unprimed variables */
      /*
      for (i=0; i< (*numvars); i++) {
	vars[i].add_var = Cudd_addNewVar(gbm);
	Cudd_Ref(vars[i].add_var);
      }

      for (i=0; i<(*numvars); i++) {
	primevars[(*numvars)-i-1].add_var = Cudd_addNewVarAtLevel(gbm,0);
	Cudd_Ref(primevars[(*numvars)-i-1].add_var);
      }
      */
      
      for (i=0; i< (*numvars); i++) {
	primevars[i].add_var = Cudd_addNewVar(gbm);
	Cudd_Ref(primevars[i].add_var);
	vars[i].add_var = Cudd_addNewVar(gbm);
	Cudd_Ref(vars[i].add_var);
      }
      

      /*Allocate memory for actions table*/
      __theMDP->Actions  = (DdNode ***) malloc(MAXACT*sizeof(DdNode**)); 
      
      for (k=0; k<MAXACT; k++)
	__theMDP->Actions[k] = (DdNode **) malloc((*numvars)*sizeof(DdNode*));
  

      if (parseActions(errf, *numvars, numactions, vars, actionlist)) {
	if ((temprew = parseDiagram(errf, *numvars, vars,0))) {
	  __theMDP->RewardD = temprew;

	  if (yylex() == DISCOUNT) {
	    if (yylex() == REAL) 
	      *discount = atof(yytext);
	    else {
	      fprintf(errf,"expected discount factor after discount keyword\n");
	      return 0;
	    }
	  }  else {
	    fprintf(errf,"expected discount or horizon keyword\n");
	    return 0;
	  }
	  if ((token = yylex()) == TOLERANCE) {
	    if (yylex() == REAL) {
	      *horizon = -1.0;
	      *tolerance = atof(yytext);
	    } else {
	      fprintf(errf,"expected tolerance factor after threshold keyword\n");
	      return 0;
	    }
	  } else if (token == HORIZON) {
	    if (yylex() == REAL) {
	      *horizon = atof(yytext);
	      *tolerance = 0.0;
	    } else {
	      fprintf(errf,"expected horizon after horizon keyword\n");
	      return 0;
	    }
	  } else {
	    fprintf(errf,"expected tolerance or horizon keyword\n");
	    return 0;
	  }
	} else {
	  fprintf(errf,"error parsing reward diagram\n");
	  return 0;
	}
      } else {
	fprintf(errf,"error parsing action diagrams\n");
	return 0;
      }
    } else {
      fprintf(errf,"error parsing variable list\n");
    }
  } else {
    fprintf(errf,"expected variable list\n");
  }
  return 1;
}



