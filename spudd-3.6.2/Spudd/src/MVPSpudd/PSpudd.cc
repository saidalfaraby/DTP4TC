/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ASPUDD: Approximate version of Spudd.....
  $Author: jhoey $
  $Date: 2003/03/12 05:59:25 $
  $Revision: 2.5 $
  $Log: PSpudd.cc,v $
  Revision 2.5  2003/03/12 05:59:25  jhoey
  *** empty log message ***

  Revision 2.4  2003/02/16 23:50:48  jhoey
  fixed a small error

  Revision 2.3  2003/02/11 23:53:04  jhoey
  *** empty log message ***

  Revision 2.2  2003/02/11 20:59:58  staubin
  Changes to debug the problem with the uploading of the files in pquerY

  Revision 2.1  2003/02/07 00:07:43  staubin
  Changes to update the manager after reading in a file...If the support does not contain all
  the variables that were in the manager at the time of the dump, I add them.

  Revision 2.0  2003/01/16 22:43:58  jhoey
  got everything working now with reading/writing, mutli-valued vars
  it also works with the original data style

  Revision 1.12  2003/01/15 19:34:22  staubin
  SLight changes made to the functions to handle the difference between numvars and numorigvars.

  Revision 1.11  2003/01/15 19:15:16  staubin
  *** empty log message ***

  Revision 1.10  2003/01/15 18:21:54  staubin
  Changes made so that my upload and store functions work like we wanted to
  It now saves all the structures necessary and it also fills the add_var field of
  the vars rnum structures....

  So the setup, once a DualOptimal.ADD is uploaded is the same as in Spudd

  Revision 1.9  2003/01/13 19:47:29  jhoey
  getting policyQuery and policyServer to work

  Revision 1.8  2003/01/11 00:41:52  jhoey
  modified the policy generation thing so it doesn't keep all the Q
  functions every iteration of vi.
  also got adding a new policy to the policyQuery gui to work
  still have to make the makefile

  Revision 1.7  2003/01/10 22:16:36  staubin
  I created the first version of the load and store ADD functions....still some slight changes to be made.

  Revision 1.6  2002/12/14 02:04:13  jhoey
  *** empty log message ***

  Revision 1.5  2002/12/13 23:29:02  jhoey
  *** empty log message ***

  Revision 1.4  2002/12/13 23:11:00  jhoey
  fixed round off method

  Revision 1.3  2002/12/12 21:48:21  jhoey
  fixed dummy variable problem

  Revision 1.1  2002/12/05 06:43:06  jhoey
  added mutli-valued variables to PSpudd

  Revision 1.4  2002/09/25 20:35:20  jhoey
  small changes for gcc 3.1
  updated dumpdot to make red lines for else arcs

  Revision 1.3  2002/09/25 20:21:21  jhoey
  CLK_TCK fixed to comply with new standard

  Revision 1.2  2002/09/21 00:55:14  jhoey
  added action costs

  Revision 1.1  2001/04/06 23:30:06  staubin
  testing....

  Revision 1.6  2001/03/21 04:36:14  jhoey
  JH: changed spanPair to span for terminal typedef stuff

  Revision 1.5  2000/11/09 17:58:54  jhoey
  Changed input argument -t slide to -t sliding

  Revision 1.4  2000/10/05 21:17:43  jhoey
  added funky log stuff at top of file



  Code by: Jesse Hoey
                 Robert St-Aubin
  
  HOW TO PRINT THE DOT FILES
  
  To produce a postscript of a dot file in solaris:
  /isd/users/dcurrie/dot/sun4/bin/dot -Tps test.dot -o test.ps   
 

  ***** VARIABLE ORDERING INFORMATION  *****************
     
      The variables are ordered alternating primed and unprimed variables. The
      LEVEL is also referred to as the order - and is changed by the
      reordering routines. The NAMES are our names for the variables - and do
      not change. The INDEX is the index, e.g. given as an argument to 
      Cudd_ReadPerm (which gives the level) - which also does not change.
      *NEW *
      LEVEL: 0              1          2          3    --------------> 2*numvars-2           2*numvars-1
      
      NAMES: primevars[0] vars[0] primevars[1] vars[1] -------------> primevars[numvars-1] vars[numvars-1]

      INDEX: 0              1          2          3    -------------> 2*numvasr-2            2*numvars-1


      So, to go  FROM         TO             USE
                INDEX        LEVEL         Cudd_ReadPerm
                LEVEL        INDEX         Cudd_ReadInvPerm
                INDEX        NAMES(nodes)  Cudd_ReadVars
                NAMES        INDEX         Cudd_NodeReadIndex

      This is what the variable reordering looked like originally (in the old
      version of spudd)

      LEVEL: 0 --------------------> numvars-1 numvars ----------> 2*numvars-1
      
      NAMES: 0 --------------------> numvars-1 0 ----------------> numvars-1
               primed_vars[i].add_var             vars[i].add_var
      
      INDEX: 2*numvars-1 <---------- numvars   0 ----------------> numvars-1

      **********************************************************************

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
#include "pspudd.hh"
#include <sys/times.h>

/* Global Variables  */

DdNode *VerySmall,*Zero,*One,*discount,*Half;
DdNode **Array1,**Array2,**Array,**ArraySum;
DdNode ***ActionsPrimed, ***NewPrime,***Allprime,***Actions;
DdNode *RewardD,*MergedMin;
DdNode *actionCost[MAXACT];
DdNode *ApproximatePolicy,*ApproximateValue, *ApproximatePolicyValue, *ApproximateMidPointValue;
DdNode *OptimalValue, *OptimalPolicy;

int allPlist[MAXACT][MAXVARS];
rnum prime_vars[MAXVARS], vars[MAXVARS];
onum orig_vars[MAXVARS];

rnum p_vars[MAXVARS],p_prime_vars[MAXVARS];
onum p_orig_vars[MAXVARS];

actt actionlist[MAXACT];

char **varnames, **actionnames, **lnames;
char *outpath,*inpath,*infile;
char linebuf[MAXLINE],name[MAXLINE],primedname[MAXLINE];

double tolerance,off,discount_factor,horizon,prune,extR;

int Leaf_Counter;  //counter used in the approximation algorithm.....
int numvars,numactions,numorigvars; //global variables of the domain
int toleranceMode;
bool mvinput;
int bigadd;
int Max_Size;
double Max_Error,SUMSQUARE;
double MAX_DIFF;
int numSumSquare;
Pair avgError;
Pair *pZero;
Pair *pOne;
Pair *toleranceP;


int cntError;

//Global Manager
/* Initialize the DdManager  */
DdManager *gbm;

int includedIn(DdManager *dd,int nsuppvars,char **Support){

  int i,j;

  for(i=0;i<numvars;i++){
    for(j=0;j<nsuppvars;j++){
      if(!strcmp(Support[j],vars[i].name))
	break;
      else if(j==nsuppvars-1)
	{
	Cudd_addIthVar(dd,vars[i].add_var->index);
	//fprintf(stderr,"We are adding a var:%d\n",vars[i].add_var->index);
	}
    }
  }
  return 1;
}



int fillRnumDdNode(DdManager *dd)
{
  int i;
  
  for(i=0;i<numvars;i++){
    //variables[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,2*i+1));
    //prime_variables[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,2*i));
    vars[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,vars[i].add_var->index));
    prime_vars[i].add_var = Cudd_ReadVars(dd,Cudd_ReadPerm(dd,prime_vars[i].add_var->index));
  }
  return 1;
}


//returns 0 if there is a discrepancy between the structures..... 1 if everything is fine.
int compareRnumStruct(rnum varsP[],rnum varsT[])
{
  int i;

  for(i=0;i<numvars;i++){

    if( strncmp(varsP[i].name,varsT[i].name,128) || (varsP[i].number != varsT[i].number) || (varsP[i].orig_num != varsT[i].orig_num) || (varsP[i].add_var->index != varsT[i].add_var->index)) 
      return 0;
  }
  return 1;
}

//Compare two structures of type onum....returns 1 if is it the same, 0 otherwise
//NOTE it assumes numvars is a global variable.....we might have to add it as a parameter if it is not the case
int compareOnumStruct(onum orig_varsP[],onum orig_varsT[])
{
  int i,j;

  for(i=0;i<numorigvars;i++){
    if( strncmp(orig_varsP[i].name,orig_varsT[i].name,128) || (orig_varsP[i].nvals != orig_varsT[i].nvals) || (orig_varsP[i].nbvars != orig_varsT[i].nbvars) || (orig_varsP[i].var1index != orig_varsT[i].var1index)){ 
      return 0;
    }
    for(j=0;j<(orig_varsP[i].nvals);j++)
      if(strcmp(orig_varsP[i].valname[j],orig_varsT[i].valname[j])){
	return 0;
      }
  }
  return 1;
}



//We copy the read structures (structIn), into the one that we will use throught the rest of the program (structOut)
int structCopy(rnum *varsIn,rnum *primeVarsIn,onum *origVarsIn)
{

  int i,j,nvals;
  
  for(i=0;i<numvars;i++){
 
    vars[i].name = strdup(varsIn[i].name);
    vars[i].number = varsIn[i].number;
    vars[i].orig_num = varsIn[i].orig_num;
    vars[i].add_var = (DdNode *)malloc(sizeof(DdNode));
    vars[i].add_var->index = varsIn[i].add_var->index; 
    
    prime_vars[i].name = strdup(primeVarsIn[i].name);
    prime_vars[i].number = primeVarsIn[i].number;
    prime_vars[i].orig_num = primeVarsIn[i].orig_num;
    prime_vars[i].add_var = (DdNode *)malloc(sizeof(DdNode));
    prime_vars[i].add_var->index = primeVarsIn[i].add_var->index;
  }

  for(i=0;i<numorigvars;i++){
    orig_vars[i].name = strdup(origVarsIn[i].name);
    orig_vars[i].nvals = origVarsIn[i].nvals;
    
    nvals = origVarsIn[i].nvals;
    for(j=0;j<nvals;j++){
      orig_vars[i].valname[j] = strdup(origVarsIn[i].valname[j]);
    }
    
    orig_vars[i].nbvars = origVarsIn[i].nbvars;
    orig_vars[i].var1index = origVarsIn[i].var1index;
  }
  
  return 1;
}

//Saves a binary files of vars, prime_vars and orig_vars structures and appends the two Optimal ADDS (value and policy)
//NOTE: the binary version of dddmp_cuddAddArrayStore is not written yet....so we save in text format.
//NEED

int writeDualOptimal(char *binfilename,DdNode *action,DdNode *value,rnum *varsT,rnum *prime_varsT,onum *orig_varsT)
{
  FILE *fp;
  int i,j,res;
  DdNode **List;
  char *basepath = "./";
  size_t length;

  
  //now let's try to save them in one file as a vector of ADD
  List  = (DdNode **)malloc(2*(sizeof(DdNode *)));
  
  List[0] = action;
  Cudd_Ref(action);
  List[1] = value;
  Cudd_Ref(value);

  //Names of the ADDs in the file saved on disk
  char **addnames = (char **) malloc(2*(sizeof(char *)));
  for(i=0;i<2; i++)
    addnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
  strcpy(addnames[0],"OptimalPolicy");
  strcpy(addnames[1],"OptimalValue");
  

  int *varid = (int *) malloc(2*numvars*(sizeof(int)));
  for(i=0;i<numvars; i++){
    varid[2*i] = varsT[i].number;
    varid[2*i+1] = prime_varsT[i].number;
  }

  /* prepares binary file for writing */ 
  if ( (fp = fopen(binfilename, "wb")) != NULL ) 
    { 
      
    
      //I put all the structures in the bin file that I am opening
      
      //We write the number of variables in the domain......
      res = fwrite(&numvars, sizeof(int), 1, fp); 
      res += fwrite(&numorigvars,sizeof(int),1,fp);
      
      
      //We write the number of actions in the domain
      res += fwrite(&numactions , sizeof(int), 1, fp); 
      //write the action names
      for(i=0;i<numactions;i++){
	length = strlen(actionnames[i]) + 1 ;
	res += fwrite(&length,sizeof(size_t),1,fp);
	res += fwrite(actionnames[i],sizeof(char),length,fp);
      }
      //I write every rnum_ structure...need to add the DdNode * later
      for(i=0;i<numvars;i++){
	length = strlen(vars[i].name) + 1 ;
	res += fwrite(&length, sizeof(size_t), 1, fp); 
	res += fwrite(varsT[i].name, sizeof(char), length, fp); 
	res += fwrite(&varsT[i].number, sizeof(int), 1, fp); 
	res += fwrite(&varsT[i].orig_num, sizeof(int), 1, fp); 
	//I am printing the index field of the DdNode * in the file now.....
	res += fwrite(&varsT[i].add_var->index, sizeof(DdHalfWord), 1, fp);
	

      }

      //I write every rnum_ structure...need to add the DdNode * later
      for(i=0;i<numvars;i++){
	length = strlen(prime_varsT[i].name) + 1;
	res += fwrite(&length, sizeof(size_t), 1, fp); 
	res += fwrite(prime_varsT[i].name, sizeof(char), length, fp); 
	res += fwrite(&prime_varsT[i].number, sizeof(int), 1, fp); 
	res += fwrite(&prime_varsT[i].orig_num, sizeof(int), 1, fp); 
	//I am printing the index field of the DdNode * in the file now.....
	res += fwrite(&prime_varsT[i].add_var->index, sizeof(DdHalfWord), 1, fp);
      }

      //I write every onum_ structure...need to add the DdNode * later
      for(i=0;i<numorigvars;i++){
	length = strlen(orig_varsT[i].name) + 1;
	res += fwrite(&length, sizeof(size_t), 1, fp); 
	res += fwrite(orig_varsT[i].name, sizeof(char), length, fp); 
	res += fwrite(&orig_varsT[i].nvals, sizeof(int), 1, fp); 
	for(j=0;j<orig_varsT[i].nvals;j++){
	  length = strlen(orig_varsT[i].valname[j]) + 1;
	  res += fwrite(&length, sizeof(size_t), 1, fp); 
	  res += fwrite(orig_varsT[i].valname[j], sizeof(char), length, fp); 
	}
	res += fwrite(&orig_varsT[i].nbvars, sizeof(int), 1, fp); 
	res += fwrite(&orig_varsT[i].var1index, sizeof(int), 1, fp); 
      }
            

      //We write the ADD in TEXT format
      Dddmp_cuddAddArrayStore(gbm,"DualAdds",2,List,addnames,varnames,varid,DDDMP_MODE_TEXT,DDDMP_VARNAMES,basepath,fp);
   
      fclose(fp); 

    }
  
  
  if(addnames){
    for(i=0;i<2;i++){
      free(addnames[i]);
    }
    free(addnames);
  }
  
  Cudd_RecursiveDeref(gbm,List[0]);
  Cudd_RecursiveDeref(gbm,List[1]);

  free(varid);
  free(List);

  return res;
}

//This function reads from a binary file with the rnum structures and the Dual ADDs
//It needs to do error checking for the case where it compares the variables with what is in the current gbm

int readDualOptimal(DdManager *dd,char *binfilename,DdNode ***DualOptimal)
{
  FILE *fp;
  int i,j,res,compres,tmp;
  
  int numa,numv,numov; //variables used until bounding the values to numvars, numactions, numorigvars
  char *temp;
  int *varid;
  char **actionREAD;
  char *basepath = "./";
  rnum varsREAD[MAXVARS];
  rnum prime_varsREAD[MAXVARS];
  onum orig_varsREAD[MAXVARS];
  size_t length;
  
  //variables to use in the HeaderLoad thing.....
  int nvars;
  int nsuppvars;
  int nRoots;
  char **orderedVarNames;
  char ** suppVarNames;
  int *varIDs;
  Dddmp_DecompType ddtype;


  char **addnames = (char **) malloc(2*(sizeof(char *)));
  for(i=0;i<2; i++)
    addnames[i] = (char *)malloc(MAXLINE*(sizeof(char)));
  
  strcpy(addnames[0],"OptimalPolicy");
  strcpy(addnames[1],"OptimalValue");
 

  /* prepares binary file for reading */ 
  if ( (fp = fopen(binfilename, "rb")) != NULL ) { 
    //read the number of variables in the domain
    res = fread(&numv,sizeof(int),1,fp);
    //read the number of original variables in the domain
    res += fread(&numov,sizeof(int),1,fp);
    //read the number of actions
    res += fread(&numa,sizeof(int),1,fp);
    
    //I am in the case where I don't even need to read the file
    if( ((numv != numvars) && (numvars != 0)) || ((numa != numactions) && (numactions != 0)) || ((numov != numorigvars) && (numorigvars != 0 ))) {
      //      fprintf(stderr,"******ERROR: PROBLEM WITH THE DOMAIN: not the same as the one uploaded already\n");
      return 0;
    }

    //read the actionnames and fill the structure with it
    actionREAD = (char**)malloc(numa *(sizeof(char *)));
    

    for(i=0;i<numa;i++){
      res += fread(&length, sizeof(size_t), 1, fp);
      actionREAD[i] = (char*)malloc(length);
      res += fread(actionREAD[i],sizeof(char),length,fp);
    }
      
    for(i=0;i<numv;i++){
      res += fread(&length, sizeof(size_t), 1, fp);
      varsREAD[i].name = (char*)malloc(length);
      res += fread(varsREAD[i].name, sizeof(char), length, fp); 
      res += fread(&varsREAD[i].number, sizeof(int),1, fp); 
      res += fread(&varsREAD[i].orig_num, sizeof(int), 1, fp); 
      //I am reading the index field of the DdNode * in the file now.....
      varsREAD[i].add_var = (DdNode *)malloc(sizeof(DdNode));
      res += fread(&varsREAD[i].add_var->index, sizeof(DdHalfWord), 1, fp);
    }

    for(i=0;i<numv;i++){
      res += fread(&length, sizeof(size_t), 1, fp); 
      prime_varsREAD[i].name = (char*)malloc(length);
      res += fread(prime_varsREAD[i].name, sizeof(char), length, fp); 
      res += fread(&prime_varsREAD[i].number, sizeof(int),1, fp); 
      res += fread(&prime_varsREAD[i].orig_num, sizeof(int), 1, fp); 
      //I am printing the index field of the DdNode * in the file now.....
      prime_varsREAD[i].add_var = (DdNode *)malloc(sizeof(DdNode));
      res += fread(&prime_varsREAD[i].add_var->index, sizeof(DdHalfWord), 1, fp);
    }
    
    for(i=0;i<numov;i++){
      res += fread(&length, sizeof(size_t), 1, fp); 
      orig_varsREAD[i].name = (char*)malloc(length);
      res += fread(orig_varsREAD[i].name, sizeof(char), length, fp);
      res += fread(&orig_varsREAD[i].nvals, sizeof(int),1, fp); 
      for(j=0;j<orig_varsREAD[i].nvals;j++){
	res += fread(&length, sizeof(size_t), 1, fp); 
	orig_varsREAD[i].valname[j] = (char*)malloc(length);
	res += fread(orig_varsREAD[i].valname[j], sizeof(char), length, fp); 
      } 
      res += fread(&orig_varsREAD[i].nbvars, sizeof(int), 1, fp); 
      res += fread(&orig_varsREAD[i].var1index, sizeof(int), 1, fp); 
    }

    //If the varsCOMP is NULL, then I need to write what I read from the file in it.
    //We do this before reading the ADDs since we might not need to read in the ADDs
        
    //This is the first domain that I enter
    if(numvars == 0){
	
      //      fprintf(stderr,"We are uploading the initial domain\n");
      //We just copy the read structures into the struct that we use
      numvars = numv;
      numorigvars = numov;
      numactions = numa;
      
      //Put the actionnames in the structure
      actionnames = (char**)malloc(numa *(sizeof(char *)));
      for(i=0;i<numa;i++)
	actionnames[i] = strdup(actionREAD[i]);
      
      //copy the structures that I just read into the permanent structures
      structCopy(varsREAD,prime_varsREAD,orig_varsREAD);
      
      //prepare to rewind the file
      //obtains the current value of the file position indicator for FILE *fp
      long filePosition = ftell(fp);
      
      //Here I read the header of the dump file to make sure that the manager had all the info needed.
      Dddmp_cuddHeaderLoad(&ddtype,&nvars,&nsuppvars,&nRoots,&orderedVarNames,&suppVarNames,&varIDs,binfilename,fp);
      
      //I include the variables that are not in the support.
      if(nsuppvars != numvars)
	includedIn(dd,nsuppvars,suppVarNames);
      
      //I go back to where I was before calling cuddHeaderLoad
      fseek(fp,filePosition,0);

    }else{   
      //There already is a domain that was uploaded, now I need to compare this new one
      //I am checking if the structures are the same, otherwise I need to return 0

      //      fprintf(stderr,"We are uploading a second structure\n");
      //compare actionnames
      for(i=0;i<numa;i++){
	if(strcmp(actionREAD[i],actionnames[i])){
	  fprintf(stderr,"THE TWO STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
	}
      }

      //compare varsREAD-vars
      if(!compareRnumStruct(vars,varsREAD))  {
	//fprintf(stderr,"THE TWO VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
      }

      //compare prime_varsREAD-prime_vars
      if(!compareRnumStruct(prime_vars,prime_varsREAD)) {
	//fprintf(stderr,"THE TWO PRIME_VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	  return 0;
      }


      //compare orig_varsREAD-orig_vars
      if(!compareOnumStruct(orig_vars,orig_varsREAD)) {
	//fprintf(stderr,"THE TWO ORIG_VARS STRUCTURES THAT YOU ARE UPLOADIN ARE NOT THE SAME \n");
	return 0;
      }
      
    }
    varid = (int *) malloc(2*(numvars)*(sizeof(int)));
    for(i=0;i<numvars; i++){
      varid[2*i] = varsREAD[i].number;
      varid[2*i+1] = prime_varsREAD[i].number;
    }
    //Loading the array of ADDs that form the graph structure to look.
    Dddmp_cuddAddArrayLoad(dd,DDDMP_ROOT_MATCHNAMES,addnames,DDDMP_VAR_MATCHIDS,orderedVarNames,varid,varid,DDDMP_MODE_TEXT,basepath,fp,DualOptimal);
    fclose(fp); 
      
  }
  
  //We fill the Rnum structures with the (DdNode *) that we are looking for....seems to work
  fillRnumDdNode(dd);

  if(addnames){
    for(i=0;i<2;i++){
      free(addnames[i]);
    }
    free(addnames);
  }
  
  if(actionREAD){
    for(i=0;i<numa;i++)
      free(actionREAD[i]);
    free(actionREAD);
  }
  
  free(varid);
  
  return res;
}




//return the maximum of two integers
int max(int i, int j)
{
  if( i< j) 
    return j;
  else
    return i;
}

// we want this function to mutliply 
DdNode * newPrimeTimes(DdManager *dd, DdNode **f, DdNode **g) {
  DdNode *res;
  DdNode *F, *G;
  double gval,fval;
  double value;
  Pair pVal;
  
  
  F = *f; G = *g;
  return F;
}
DdNode * MyPlus(DdManager * dd, DdNode ** f, DdNode ** g)
{
  DdNode *res;
  DdNode *F, *G;
  double gval,fval;
  double value;
  Pair pVal;
  
  
  F = *f; G = *g;
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    gval = (*Cudd_V(G)).get_min();
    fval = (*Cudd_V(F)).get_min();
    
    if (fval == 0.0) return(G);
    
    return(F);
  }
  return(NULL);

}

DdNode * addMean(DdManager * dd, DdNode ** f, DdNode ** g)
{
    DdNode *F, *res;
    double fmn,fmx;
    double value;
    Pair pVal;


    F = *f; 
    if (Cudd_IsConstant(F)) {
      fmn = (*Cudd_V(F)).get_min();
      fmx = (*Cudd_V(F)).get_max();
      
      pVal.set((fmn+fmx)/2.0);
      res = Cudd_addConst(gbm,&pVal);
      return(res);
    }
    return(NULL);

}

/* PickOneAction not used anymore since we don't use
the binary representation - JH OCT 3 00 */
/* this function will work for terminal nodes - g is the action
to pick (a number starting at 0) and the f DdNode should be
the terminal nodes of the action diagram: 101110101010 type things
The result is a 1 if the first 1 in the action diagram terminal node
is in the position specified by that action, 0 otherwise 
returns NULL if f and g are not constant*/
DdNode * PickOneAction(DdManager * gbm, DdNode ** f, DdNode ** g)
{
    DdNode *res;
    DdNode *F, *G;
    double gval, fval;
    int act=0;
    
    F = *f; G = *g;


    if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
      gval = (*Cudd_V(G)).get_min();
      fval = (*Cudd_V(F)).get_min();
      if (fval == 0.0) {
	fprintf(stderr,"No action detected - abort\n");
	exit(0);
      }
      while (fmod(fval,10.0) != 1.0) {
	fval = floor(fval/10.0);
	act++;
      }
      if (act == gval){ 
	return Cudd_addConst(gbm,pOne);
	
      }      
      else{
	return Cudd_addConst(gbm,pZero);
      }
    }
    return(NULL);

}

// take all vcurrent and action diag and merges them by picking the right action
DdNode *actionMerge(DdNode **valD, DdNode *actD, int numact)
{
  int i;
  DdNode *pickAct, *temp, *temp1, *temp2, *sum;
  Pair iPair;

  pickAct = Cudd_addConst(gbm,pOne);
  Cudd_Ref(pickAct);
  
  sum = Cudd_addConst(gbm,pZero);
  Cudd_Ref(sum);
  
  temp1 = Cudd_addConst(gbm,pOne);
  Cudd_Ref(temp1);

  for (i=0; i<numact; i++) {
    
    iPair.set((double) i+1);
    Cudd_RecursiveDeref(gbm,pickAct);
    pickAct = Cudd_addConst(gbm,&iPair);
    Cudd_Ref(pickAct);

    //temp = Cudd_addApply(gbm, PickOneAction, actD, pickAct);
    temp = Cudd_addApply(gbm, addAddExact, actD, pickAct);
    Cudd_Ref(temp);
    
    Cudd_RecursiveDeref(gbm, temp1);
    temp2 = valD[i];
    Cudd_Ref(temp2);

    temp1 = Cudd_addApply(gbm, Cudd_addTimes, temp,temp2);
    Cudd_Ref(temp1);
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm, temp);

    /** JH FEB 4 **/
    
     temp = Cudd_addApply(gbm, Cudd_addPlus, temp1, sum);
    //temp = Cudd_addApply(gbm, MyPlus, temp1, sum);
    Cudd_Ref(temp);
    //Cudd_RecursiveDeref(gbm,temp1);
    Cudd_RecursiveDeref(gbm,sum);
    sum = temp;
    
    //Cudd_Ref(sum);
    //Cudd_RecursiveDeref(gbm, temp);
  }

  Cudd_RecursiveDeref(gbm,temp1);
  Cudd_RecursiveDeref(gbm,pickAct);
    
  return sum;
}

/*
Loop like policy Iteration....
Backs up the value function valD using the
policy specified by actD
*/
DdNode *evalPolicy(DdNode *valD, DdNode *actD)
{
  int i,test,counter;
  DdNode **Vcurrent,*temp3;
  DdNode *Vpast,*temp,*temp1, *temp2;
  DdNode *Merged_Add;
  DdNode *tempA,*tempB;
  
  Vcurrent = (DdNode **)malloc(numactions*(sizeof(DdNode*)));
  
  test = 0;
    
  for(i=0;i<numactions;i++) {
    Vcurrent[i] = valD;
    Cudd_Ref(Vcurrent[i]);
  }  
  
  Vpast = valD;
  Cudd_Ref(Vpast);
  
  Merged_Add = valD;
  Cudd_Ref(Merged_Add);
  
  counter = 0;
  while ((horizon == -1 && test==0) || (horizon >= 0 && counter < horizon)) // && (counter < 50))
  { 
    // KEEP TRACK OF THE NUMBER OF ITERATIONS 
      counter = counter + 1;
      //fprintf(stderr,"iteration %d\n",counter);
      // KEEPING A COPY OF THE LAST ADD TO COMPARE FOR CONVERGENCE 
      Cudd_RecursiveDeref(gbm,Vpast);
      Vpast = Merged_Add;
      Cudd_Ref(Vpast);

      temp = Cudd_addSwapVariables(gbm,Merged_Add,Array1,Array2,numvars);
      Cudd_Ref(temp);
      //temp = Merged_Add.SwapVariables(ADDvector(numvars,&gbm,Array1),ADDvector(numvars,&gbm,Array2));
      for(i=0;i<numactions;i++){
#ifdef OPTIMIZED
	temp1 = primeReward(temp,&(Allprime[i][0]),prime_vars,0,&(allPlist[i][0]),numvars,i);
#else
	temp1 = newPrimeReward(temp,NewPrime[i]);
#endif
          
	temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp1,discount);
	Cudd_Ref(temp2);
	Cudd_RecursiveDeref(gbm,temp1);
	  
	temp1 = Cudd_addApply(gbm,Cudd_addPlus,temp2,RewardD);
	Cudd_Ref(temp1);
	Cudd_RecursiveDeref(gbm,temp2);
	
	Cudd_RecursiveDeref(gbm,Vcurrent[i]);
	Vcurrent[i] = temp1;
	
      }
      Cudd_RecursiveDeref(gbm,temp);
      Merged_Add = actionMerge(Vcurrent, actD, numactions);
      Cudd_Ref(Merged_Add);
      // COMPARE THE TWO MOST RECENT ADDs IN ORDER TO VERIFY THE CONVERGENCE 
      //test = Merged_Add.EqualSupNorm(Vpast,toleranceP,0);
      test= Cudd_EqualSupNorm(gbm,Merged_Add,Vpast,toleranceP,0);
    }
  
  for (i=0; i<numactions; i++)
    Cudd_RecursiveDeref(gbm, Vcurrent[i]);
  free(Vcurrent);
  
  return Merged_Add;
}

// Wrapper to get around some problems in Cudd
DdNode *Cudd_Else(DdNode *foo)
{
        if (Cudd_IsComplement(foo)) {
                return Cudd_Not(Cudd_E(foo));
        } else {
                return Cudd_E(foo);
        }
}

// Wrapper to get around some problems in Cudd
DdNode *Cudd_Then(DdNode *foo)
{
        if (Cudd_IsComplement(foo)) {
                return Cudd_Not(Cudd_T(foo));
        } else {
                return Cudd_T(foo);
        }
}
/*
Sums over the primed variables in dd from index first to index last
*/
DdNode *sumSubtrees(DdNode *dd, int first, int last,  rnum *prime_vars, int numvars)
{
  DdNode *tempT, *tempE, *temp, *cube;
  int i,j,index;

  j=0;
  //JH09/05/00  for (i=first-1; i>=last; i--) {
  for (i=first; i<last; i+=2) {
    ArraySum[j] = prime_vars[i/2].add_var;
    Cudd_Ref(ArraySum[j]);
    j++;
  }
  //JH09/05/00 cube = Cudd_addComputeCube(gbm,ArraySum,NULL,first-last);
  cube = Cudd_addComputeCube(gbm,ArraySum,NULL,(last-first)/2);
  Cudd_Ref(cube);

  j=0;
  for (i=first; i<last; i+=2) {
    Cudd_RecursiveDeref(gbm,ArraySum[j]);
    j++;
  }

  temp = Cudd_addExistAbstract(gbm,dd,cube);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,cube);
  return temp;
}
/*
Sums over the primed variables in dd from index first to index last
for multi-valued diagrams
*/
DdNode *newSumSubtrees(DdNode *dd, int first, int last)
{
  DdNode *tempT, *tempE, *temp, *cube;
  int i,j,index;

  j=0;
  for (i=first; i<last; i++) {
    ArraySum[j] = prime_vars[i].add_var;
    Cudd_Ref(ArraySum[j]);
    j++;
  }
  cube = Cudd_addComputeCube(gbm,ArraySum,NULL,(last-first));
  Cudd_Ref(cube);

  j=0;
  for (i=first; i<last; i+=2) {
    Cudd_RecursiveDeref(gbm,ArraySum[j]);
    j++;
  }

  temp = Cudd_addExistAbstract(gbm,dd,cube);
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,cube);
  return temp;
}

// old way
// mutliply reward by each dual action diagram and sum out subtrees
// only works for original-style (boolean) input
DdNode *newPrimeReward(DdNode *reward, DdNode **NewPrime) {
  
  int j;
  DdNode *temp, *temp1, *cube;
  int first, last;

  temp1 = reward;
  Cudd_Ref(temp1);
  for (j=0; j<numorigvars; j++) {
    temp = Cudd_addApply(gbm,Cudd_addTimes,NewPrime[j],temp1);
    Cudd_Ref(temp);
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;

    if (mvinput) {
      first = orig_vars[j].var1index;
      last = orig_vars[j].var1index+orig_vars[j].nbvars;      
      temp = newSumSubtrees(temp1,first,last);
    } else {
      cube = Cudd_addIte(gbm,prime_vars[j].add_var,One,Zero);
      Cudd_Ref(cube);
      
      temp = Cudd_addExistAbstract(gbm,temp1,cube);
      Cudd_Ref(temp);
      Cudd_RecursiveDeref(gbm,cube);
    }
    Cudd_RecursiveDeref(gbm,temp1);
    temp1 = temp;
  }
  return temp1;
}


/*
This does the multiplication of Pr(s,a,t) by V(t) summed over t
in the value iteration loop. Pr(s,a,t) is in AllPrime as an array of
the chunks of the complete action diagram. 
It does so by descending V(t) over the 'chunked'
levels to the lowest one, where it multiplies by the corresponding
chunk in Allprime, followed by a sum over the primed variables in 
the chunk. See section 4.2 in the Tech Report 
*/
DdNode  *primeReward(DdNode *reward, DdNode **AllPrime, rnum *prime_vars, 
                     int n, int *aplist, int numvars, int act)
{
  int i,j,count,test,dagCount,dagLeafCount,treeCount,treeLeafCount;
  DdNode *temp,*temp1,*tempvar,*prodT, *prodE, *tempT, *tempE, *totT, *totE, *tot, *R;
  DdNode *node;
  test = 0;

  /* Get the index of the root node of the reward tree */
  if ((Cudd_IsConstant(reward) == 1)) {
    temp = reward;
    Cudd_Ref(temp);
  } else {
    i = Cudd_NodeReadIndex(reward);
    
    //JH09/05/00if (i < aplist[n]) {
    // if we're at the end - do the multiplication
    if (i >= aplist[n]) {
      //JH09/05/00      if (aplist[n+1] == numvars) {
      if (aplist[n+1] == 2*numvars) {
        temp = reward;
        Cudd_Ref(temp);
      } else{
        temp = primeReward(reward, AllPrime, prime_vars, n+1, aplist, numvars,act);
      }
     
      temp1 = Cudd_addApply(gbm,Cudd_addTimes,AllPrime[n],temp);
      Cudd_Ref(temp1);
      
      Cudd_RecursiveDeref(gbm,temp);

      temp = sumSubtrees(temp1,aplist[n],aplist[n+1],prime_vars,numvars);
      Cudd_RecursiveDeref(gbm,temp1);
      
    } else {
      // descend to the next level
      tempT = primeReward(Cudd_Then(reward), AllPrime, prime_vars, n, aplist, numvars,act); 
      tempE = primeReward(Cudd_Else(reward), AllPrime, prime_vars, n, aplist, numvars,act); 
      //     temp = Cudd_addIte(gbm, prime_vars[2*numvars-1-i].add_var, tempT, tempE);
      // recombine the results
      temp = Cudd_addIte(gbm,prime_vars[i/2].add_var,tempT,tempE);
      Cudd_Ref(temp);
      /* deref these because I think that Cudd_addIte does its own reffing! */
      Cudd_RecursiveDeref(gbm, tempT);
      Cudd_RecursiveDeref(gbm, tempE);
    }
  }
  return temp;
}
/*
Gets the policy from the array of Vcurrents (value functions for each action) and the Vmerged (the
combined value function with the maximum values). The values at the terminal nodes of Vmerged correspond
to values at the terminal nodes of one or more of the Vcurrents. This correspondence is what we use
to figure out the policy (i.e. the value of which action was used to give the Vmerged when doing the
maximization). We do it as follows
1. subtract the kth Vcurrent from Vmerged - the nodes which match are 0 in the result
2. turn the resulting add into a BDD by thresholding at 0 - now the 1 nodes of the BDD are the 
   states for which the kth action was the one that contributed to Vmerged
3. turn that BDD back into a 0-1 ADD
4. offset the value by 1e+k (this is to keep track of ALL the possible actions)
5. maintain a running sum.
*/
DdNode *getActionAdd(DdNode **Vcurrent, DdNode *Vmerged, int numactions)
{
  double off;
  DdNode *actionAdd;
  actionAdd = Zero;
  Cudd_Ref(actionAdd);
  off = 1.0;
  getActionAdds(Vcurrent,Vmerged,numactions,&actionAdd,off);
  return actionAdd;
}
// this function takes the initial action ADD which in the general case (see getActionAdd above)
// is just Zero, but when we are maximizing a single action, it will be the current policy
// the initial offset is the current action index (usually 1.0) but will be the action number when we
// maximize a single action. All this is to avoid keeping the Q functions around during value iteration.
void getActionAdds(DdNode **Vcurrent, DdNode *Vmerged, int numactions,DdNode **actionAdd,double off)
{
  int i,j,k;
  DdNode *temp,*temp2,*offset;
  Pair offPair;
  /* loop over actions */
  for (k=0;k<numactions;k++) {

    temp = Cudd_addApply(gbm,Cudd_addMinus,Vcurrent[k],Vmerged);
    Cudd_Ref(temp);

    
    temp2 = Cudd_addBddThreshold(gbm,temp,pZero);
    Cudd_Ref(temp2);

 
    Cudd_RecursiveDeref(gbm,temp);
    temp = Cudd_BddToAdd(gbm,temp2);
    Cudd_Ref(temp);

   /* temp is a 0-1 ADD with all 1's are states where the
       kth action generated the value in Vmerged */

    Cudd_RecursiveDeref(gbm,temp2);
    offPair.set(off);
    offset = Cudd_addConst(gbm,&offPair);  
    Cudd_Ref(offset);
    
    temp2 = Cudd_addApply(gbm,Cudd_addTimes,temp,offset);
    Cudd_Ref(temp2);
    /*temp2 is the same structure as temp but has the
      number 1e+k instead of the 1 nodes*/


    Cudd_RecursiveDeref(gbm,temp);
    temp = Cudd_addApply(gbm,MyPlus,temp2,actionAdd[0]);
    Cudd_Ref(temp);

    
    Cudd_RecursiveDeref(gbm,temp2);
    Cudd_RecursiveDeref(gbm,actionAdd[0]);
    actionAdd[0] = temp;

    off = off+1.0;
  }
}
//Counts the number of internal nodes for the tree equivalent to the Ddnode *x
int count_internal_nodes(DdNode *x)
{
  DdNode *tempE, *tempT;

  int counter=0,leftcounter=0,rightcounter=0;
  
  if (!Cudd_IsConstant(x)) {
    tempE = Cudd_Else(x);
    tempT = Cudd_Then(x);
    counter = 1;
    leftcounter = count_internal_nodes(tempT);
    rightcounter = count_internal_nodes(tempE);
  }

  counter += leftcounter+rightcounter;
  return counter;
}

//Counts the number of leaves for a Tree similar to the Ddnode *x
int count_leaves(DdNode *x)
{
  DdNode *tempE, *tempT;

  int counter=0,leftcounter=0,rightcounter=0;
  
  if (Cudd_IsConstant(x)) {
    counter = 1;
  } else {
    tempE = Cudd_Else(x);
    tempT = Cudd_Then(x);

    leftcounter = count_leaves(tempT);
    rightcounter = count_leaves(tempE);
    counter = leftcounter+rightcounter;
  }
  return counter;
}
DdNode *SumSquareError(DdManager * gbm, DdNode ** f, DdNode ** g)
{
  DdNode *F, *G;
  DdNode *SSE;
  F = *f; G = *g;
  
  if (Cudd_IsConstant(F) && Cudd_IsConstant(G)) {
    numSumSquare++;
    SUMSQUARE = SUMSQUARE + ( (*Cudd_V(F)).get_max() * (*Cudd_V(G)).get_max());
    return(One);
  } 
  return(NULL);
} 
/* To dump the output as an HTML format*/
void dumpHtmlRew(FILE *log,char *outpath)
{
  char temp[2*MAXLEN];
 
  fprintf(log,"<html>\n\t<head>\n");
  fprintf(log,"<title>Reward</title>");
  fprintf(log,"</head>\n\n\t<body>\n\t<h1>Reward</h1>\n");
  strcpy(temp,"<img src=http://www.research.att.com/~north/cgi-bin/webdot.cgi/http://www.cs.ubc.ca/spider/jhoey/spudd/");
  strcat(temp,outpath);
  strcat(temp,"reward.dot.gif>\n");
  fprintf(log,temp);
  fprintf(log,"\n\n\t<hr>\n</body>\n</html>\n");
} 
void shuffleRandom() 
{
  /* generate random order */
  int i;
  double *infor = (double *) malloc(numvars*sizeof(double));
  int *reordList = (int *) malloc(2*numvars*sizeof(int));
  
  for (i=0; i<numvars; i++) {
    infor[i] = ((double) rand())/((double) RAND_MAX);
    //fprintf(stderr,"%f ",infor[i]);
  }
  //fprintf(stderr,"\n");
  
  infoSorter(infor,reordList);
  (void) Cudd_ShuffleHeap(gbm,reordList);    
  free(infor);
  free(reordList);
}

void testReorder(DdNode *OptimalValue, int reorderMeth, FILE *dataFile, FILE *dataFile2) 
{
  int optCount,reoCount;
  double temps0,temps1,temps;

  struct tms _t;
  long clk_tck = sysconf(_SC_CLK_TCK);

  fprintf(stderr,"\n********** EXPERIMENT 1 DATA ***********\n");
  optCount = Cudd_DagSize(OptimalValue);
  //  fprintf(dataFile,"Original Size: %d ",optCount);
  fprintf(dataFile2,"%d ",optCount);  
  shuffleRandom();
  optCount = Cudd_DagSize(OptimalValue);

  
  times(&_t);
  temps0 = 1.0*_t.tms_utime/clk_tck;
  
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
      if (!reorderMinSpan(OptimalValue))
	fprintf(stderr,"******ERROR: Reorder minSpan failed \n");
      break;
    case REORDER_EXACT:
      if (!Cudd_ReduceHeap(gbm, CUDD_REORDER_EXACT,0))
	fprintf(stderr,"******ERROR: Reorder exact failed \n");
      break;
    case REORDER_NONE:
      break;
    default:
      break;
    }
  /* We get the time of completion*/
  times(&_t);
  temps1 = 1.0*_t.tms_utime/clk_tck;
  
  temps = temps1- temps0;
  reoCount = Cudd_DagSize(OptimalValue);

  /*  fprintf(dataFile,"Size before Reordering: %d After: %d\n",optCount,reoCount);
  fprintf(dataFile,"Time to reorder: %f\n",temps);
  */
  fprintf(dataFile2,"%d %d ",optCount,reoCount);
  fprintf(dataFile2,"%f ",temps);
  fprintf(stderr,"****************************************\n\n");
} 
void writeVPDot(DdNode *valADD, DdNode *polADD,char *fileprefix) {
  FILE *value, *action;
  char fname[256];
  // write out data files
  sprintf(fname,"%svalue.dot",fileprefix);
  value = fopen(fname,"w");
  if (mvinput) {
    DumpDot(gbm, valADD, vars, orig_vars, NULL, value);
  } else {
    Array[0] = valADD;  
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars],NULL,value);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(value);
  
  sprintf(fname,"%spolicy.dot",fileprefix);
  
  action = fopen(fname,"w");
  
  if (mvinput) {
    DumpDot(gbm, polADD, vars, orig_vars, actionnames,action);
  } else {
    Array[0] = polADD;
    Cudd_Ref(Array[0]);
    My_DumpDot(gbm,1,Array,varnames,&lnames[numvars+1],actionnames,action);  
    Cudd_RecursiveDeref(gbm,Array[0]);
  }
  fclose(action);
}
