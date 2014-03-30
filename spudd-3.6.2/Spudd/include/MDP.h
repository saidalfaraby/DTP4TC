#ifndef __MDP
#define __MDP
 // MDP CLASS 
// By Jesse Hoey & Robert St. Aubin
// this class encapsulates an MDP and represents/does calculations on the
// components of the MDP using ADDs. This code was ripped from MVPSpudd in April 2003
// and esentially combines PSpudd.cc and JAPRICODD.cc into a single file
// certain things (like the parser) had to be mangled a bit - for instance
// it needs to have a global pointer the MDP object it is parsing.
// net2 library for socket communication
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <sys/times.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <map>

#include "util.h"
#include "st.h"
#include "cuddInt.h"

#include "cudd.h"
//#include "cuddObj.hh"
#include "dddmpInt.h"
#include "dddmp.h"
#include "net2.h"
#include "networkSettings.h"



// this is where the client will be running
#define SUPERVISOR_HOST_TMP "128.100.3.14"

//#define DEBUGP

#define MAXLINE 256	/* Maximum length of each input line read. */
#define MAXVARS 64
#define MAXVALS 1000
#define MAXLIN 1024
#define MAXACT 256
#define MAXLEN 64
#define DUMPIN 0
#define MAXMEM 700000000
#define MAXMEMHARD 2000000000 
#define BIGADD 10000
#define HREWARD 270
#define HACTION 271
#define HVALUE 272
#define HPOLICY 273
#define THRESHOLD  1.0     /* Maximum error allowed before prunning an ADD */
#define MAX_SIZE  10000
#define MAX_ERROR 0.00001
#define LIMIT_ERROR 371
#define ERROR_TOL 0.1
#define LIMIT_SIZE  372
#define SIZE_TOL  2
#define LARGEFLAG 1000
#define REORDERAPPLY_FIRSTFIVE 1500
#define REORDERAPPLY_ALL      1501
#define REORDER_SIFT 472
#define REORDER_RANDOM 473
#define REORDER_MINSPAN 474
#define REORDER_NONE 475
#define REORDER_EXACT 476
#define REORDER_SP_SIFT 477
#define APPROX_NONE 571
#define APPROX_ALLPAIRS 572
#define APPROX_ROUNDOFF 573
#define OPTIMAL_NONE 671
#define OPTIMAL_FILE 672
#define OPTIMAL_GENERATE 673
#define TOLERANCE_FIXED 772
#define TOLERANCE_SLIDING 773

//#define OPTIMIZED
//#define ALLACTIONS 1
//#define COSTSPRIMED 1

#define MAXNUMPOLICIES 20

using namespace std;
  struct ltstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    }
  };
  


int yyparse();
/* lex stuff */
extern FILE *yyin;
extern char *yytext;
extern int yydebug;
extern int yylex();
/* other structures for variables and actions */

typedef struct list_pair
{
  DdNode *add;
  struct list_pair *next;
} list_cst;

struct onum_ {
  char *name;
  int nvals;
  char *valname[MAXVALS];
  int nbvars;
  int var1index;
};
typedef struct onum_ onum;

struct rnum_ {
  char *name;
  int number;
  int orig_num;
  DdNode *add_var;
};

typedef struct rnum_ rnum;

struct actt_ {
  char *name;
  int index;
};

typedef struct actt_ actt;


DdNode *Cudd_Else(DdNode *);
DdNode *Cudd_Then(DdNode *);
DdNode *sumSubtrees(DdNode *, int, int,  rnum *, int);
DdNode  *primeReward(DdNode *, DdNode **, rnum *, int, int *, int, int);
DdNode  *newPrimeReward(DdNode *, DdNode **);

DdNode *myRoundOff(DdManager *, DdNode **, DdNode **);
DdNode *getTotalSpan(DdManager *, DdNode **, DdNode **);
DdNode *addAddExact(DdManager *dd, DdNode **, DdNode **);
DdNode * My_cuddAddApplyRecur(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode *My_addApply(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode * MyPlus(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode * PickAction(DdManager *dd, DdNode **f, DdNode **g);
DdNode * BinaryAction(DdManager *dd, DdNode **f, DdNode **g);

DdNode * getAction(DdManager *dd, DdNode **f, DdNode **g);

DdNode *getErrorAdd(DdManager *dd, DdNode **f, DdNode **x);
DdNode *My_addMaximum(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode * addMean(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *addBitwiseCompare(DdManager *dd , DdNode **f, DdNode **g);

DdNode * reducePrecision(DdManager *dd, DdNode **f, DdNode **g);


void pQuery(DdManager *dd, DdNode *value, DdNode *action, char * pfile, onum *orig_vars, rnum *vars, int, int);
void pdd(DdNode *b);

double  get_extent(DdNode *);
double  get_error(DdNode *);
Pair get_span(DdNode *);
int close_enough(double , double , double ); 
char* DumpDoth(DdManager *dd, DdNode *add, rnum *vars, int nv, onum *orig_vars, int nov, int rovar, char **,  FILE *fp, FILE *nfp, int *hmn);
char* DumpDoth_p(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int nv, onum *orig_vars, int nov, int rovar, char **,  FILE *fp, FILE *nfp, int *hmn);
void writeConstantNode(FILE *fp, char *nname, DdNode *add, char **lnames);
int DumpDot(DdManager *dd, DdNode *add, rnum *vars, int nv, onum *orig_vars, int nov, char **, FILE *fp);
int DumpDot_p(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int nv, onum *orig_vars, int nov, char **, FILE *fp);
int My_DumpDot(DdManager *,int  n , DdNode **, char **,char **,char **, FILE *);
int printDdNode(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp);
void printDdNode(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars,  FILE *fp, char *tabstop);
void printDdNode(DdManager *dd, DdNode **add, int numdds, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp);
int printDdNode(DdManager *dd, DdNode *add, rnum ***addlist, int & numadds, int & numexpadds, 
		rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars);
int storeAddInList(DdManager *dd, rnum ***addlist, int & numadds, int & numexpadds, rnum *add);
int findAddInList(rnum ***addlist, int & numadds, int & numexpadds, DdNode *add);
int printDdNode(DdManager *dd, rnum ***addlist, int & numadds, int & numexpadds, 
		DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp);


int aconvert(char * names, char **lnames, double lval, char *separator="\\n");

int pickAction(double lval);
int binaryAction(int lval);

void print_out(DdNode *, int, rnum *, DdNode *, DdNode *, int);

void dumpHtmlRew(FILE *log,char *outpath);
void infoSorter(double *infoList, int *reorderList, int numvars);
void removeElement(list_cst *list,list_cst *element);
DdNode *replacePairs(DdManager *dd, DdNode **f, DdNode **newC);

/*parser routines */
int parseVars(FILE *,int *, rnum *, rnum *);
int findVar(int *, int, char *, rnum *);
DdNode *parseADD(FILE *,int, rnum *,int);
DdNode *parseDiagram(FILE *errf, int, rnum*, int);
int parseActions(FILE *,int, int *, rnum *, actt *);
int pparse(FILE *,double *, double *, double *, int *, int *, rnum *, rnum *, actt *);

void addVar(onum *, int, char*);
int findOVar(onum *ovar, int nov, const char *findname);
int findOVal(onum *ovar, int novar, const char *findval);
void setActionCost(DdNode **ac, int nac, double val);
DdNode * buildCubeCPT(rnum *v, int coval, onum ovar, DdNode *add);
void removeDummy(DdNode **np, int);
void removeAllDummys(DdNode **np, int);
void removeDummyp(DdNode **np, int);
void removeAllDummysp(DdNode **np, int);
void addDummyStates(DdNode **np, int nv);

void buildGoodStateADDs(onum *, rnum *, rnum *, int);
DdNode *renormalizeCube(DdNode *, double);
bool member(int *, int, int);

void decodevvals(int *varvals, onum *orig_vars, int numorigvars, int state);


// non-MDP functions
// some of these are non-MDP only because they need access  to
// some of the global variables defined below. There is surely a more
// elegant way to do this.
DdNode *Convergence_Test(DdManager * gbm, DdNode ** f, DdNode ** g);
DdNode *Cudd_Else(DdNode *);
DdNode *Cudd_Then(DdNode *);
DdNode *sumSubtrees(DdNode *, int, int,  rnum *, int);


DdNode *sumOutPrime(DdNode *dd, onum *ovar, rnum *prime_vars);

// sums out primed variables primed_vars[indices[i]] for i=1..numindices from dd
DdNode *sumOutPrime(DdNode *dd, int *indices, int numindices, rnum *prime_vars, int numvars);
// sums out all primed variables in toSum
DdNode * sumPrimes(DdNode *toSum, rnum* , int numvars);


DdNode *newSumSubtrees(DdNode *dd, int first, int last, rnum *);
DdNode *actionMerge(DdNode **valD, DdNode *actD, int numact);
double length_span(DdNode *x,DdNode *y);
int count_internal_nodes_tree(DdNode *);
int count_leaves_tree(DdNode *);
void getAV(DdManager *dd, DdNode *val, DdNode *act, Pair & dval, Pair & aval, int *varvals,
	   rnum *, int, onum *, int);

// gets variable assignment for varvals
void getVarAss(int *varass, int *varvals, rnum *v, int nvars, onum *ov, int novars);
void getVarAssNoZero(int *varass, int *varvals, rnum *v, int nvars, onum *ov, int novars);



// non-MDP functions for comparing values of a variable in an ADD
int compareVarVals(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		   int ovar, double tol, int *compmat) ;
int compareVarVals(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		   int ovar, int ovarval1, int ovarval2, double tol);
DdNode *restrictVal(DdManager *dd, DdNode *val, rnum *v, int nvars, onum *ov, int novars,
		    int ovar, int ovarval);

//  other non-MDP functions used when reading in MDPs
int compareRnumStruct(rnum varsP[],rnum varsT[], int nv);
int compareOnumStruct(onum orig_varsP[],onum orig_varsT[], int nov);


// max of two ints 
int max(int, int);

// for starting a server using the net2 library (homer)
int connectToSupervisor     (Client **spudd_c, Server **spudd_s);

// gbm must be global
extern DdManager *gbm;
extern Pair tempOne,tempTwo;
extern double tolerance;
extern Pair *pZero, *pOne;
extern  DdNode *One, *Zero;

class MDP {
 protected:
  double horizon, discount_factor, prune, extR, bigadd, Max_Size, Max_Error;
  int numvars, numactions, numorigvars, counter;
  int toleranceMode,approxMeth,reorderMeth,reorderApplyMeth,shuffleFlag,limit_type;
  int ***ovarsupport;
  Pair * toleranceP;
  bool mvinput, delete_gbm;
  // for timing
  double temps;
  long int maxusememory;

  int allPlist[MAXACT][MAXVARS];
  rnum prime_vars[MAXVARS], vars[MAXVARS];
  onum orig_vars[MAXVARS];
  
  rnum p_vars[MAXVARS],p_prime_vars[MAXVARS];
  onum p_orig_vars[MAXVARS];
  
  actt actionlist[MAXACT];
  
  char **varnames, **actionnames, **lnames;
  char linebuf[MAXLINE],name[MAXLINE],primedname[MAXLINE];
  
  DdNode *VerySmall,*discount,*Half;
  DdNode **Array1,**Array2,**Array,**ArraySum;
  DdNode ***NewPrime,***Allprime, ***Actions;
  DdNode *RewardD,*RewardDNoDummy,*MergedMin;
  DdNode *actionCost[MAXACT], *actionCostNoDummy[MAXACT];
  DdNode *Qf[MAXACT];
  DdNode *Merged_Add;
  // from JAPRICODD
  DdNode *VMinpast;
  DdNode **VcurrentMin;
  DdNode *RewardMax;
  DdNode *nodeOne,*nodeTwo;


  // actually gets the value function
  virtual void getValueFunction();
  virtual void allocateMemory();
  virtual void buildAllPrime();

  // intializes some stuff - including the global ddManager (gbm)
  void init();

  // does the approximation
  void approximate();

  // helpers for approximate spudd
  void size_approx();

  // all pairs approximation method
  DdNode *roundOffMMA(DdNode *res);
  void allPairsApprox(DdNode **theADD);

  // reorder methods
  int reorderMinSpan(DdNode *addBase);
  int close_enough(double val, double limit, double tol) ;

  // calls parser
  virtual int parseInput();

  // policy server
  virtual void policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars);

  virtual void policyServe2(DdManager *gb, DdNode *act, DdNode *val, int numovars);
  
  // helpers for reading/writing dual ADD files
  int includedIn(DdManager *dd,int nsuppvars,char **Support);
  int fillRnumDdNode(DdManager *dd);
  int structCopy(rnum *varsIn,rnum *primeVarsIn,onum *origVarsIn);

  // helpers for value iteration /policy generation
  DdNode  *newPrimeReward(DdNode *, DdNode **);
  DdNode  *primeReward(DdNode *, DdNode **, rnum *, int, int *, int, int);

  void recoverPolicy(DdNode *, DdNode *, DdNode *, DdNode **, double);

  // replace newPrimeReward
  DdNode *multiplySumSet(DdNode *reward, DdNode **primes);
  DdNode *multiplySum(DdNode *reward, DdNode *prime, DdNode *varsum, onum *ovar);

  DdNode *sumOutAllOrigVars(DdNode *dd);


  
  //selects one branch probabilistically according
  //to the distribution over all possible states (which should sum to 1)
  //and returns the index of that branch in the original variable index ordering
  int selectBranchProb(DdNode *newVarDist, int ovar);

 public:
  // pointers to the policy and value functions
  DdNode *ApproximatePolicy, *ApproximateValue, *ApproximatePolicyValue, *ApproximateMidPointValue;
  DdNode *OptimalValue, *OptimalPolicy;

  // Default constructor
  MDP();

  // copy constructor
  MDP(MDP *mdp);

  // Constructs and reads in an ADD
  MDP(char *infile, double badd=BIGADD);

  //Desctructor
  ~MDP();
  

  // reads in an MDP 
  // if dualOptimal = false, this calls readMDPSpec (reads MDP specs using parser)
  // if dualOptimal = true - this reads a dual optimal MDP file using readDualOptimal (not implemented yet)
  int readMDP(char *infile, bool dualOptimal=false); 

  //data specification from file inpath/infile
  // opens the errf file int outpath/infile-error.dat for writing any error messages
  virtual int readMDPSpec(char *infile);

  // generates the policy as requested by all these flags


  // with no arguments, generatePolicy() generates the optimal policy using 'vanilla' Spudd
  // generatePolicy(rMeth) generates optimal policy with variable reordering according to rMeth
  // rAMeth can be set to REORDER_FIRSTFIVE to only reorder on the first five iterations
  // generatePolicy with at least 6 arguments generates approximate policies
  // sFlag is 1 if you want to shuffle the variable ordering before starting
  // ofile is where the dual optimal policy/value file gets written
  // outpath/infile-OPT/APP files are for dot representations, which can be
  // dumped if dFlag = true
  // for more information, see the usage directions in testmdp.cc, or see the web page
  virtual int generatePolicy(int rMeth = REORDER_NONE,
		     int rAMeth = REORDERAPPLY_ALL,
		     double tolerance=0.1,
		     int badd=BIGADD,
		     float prun=0.1,
		     int max_size=MAX_SIZE,
		     int aMeth = APPROX_NONE,
		     int ltype = LIMIT_SIZE,
		     int tMode = TOLERANCE_FIXED,
		     int sFlag = 0,
		     bool dFlag = false);

  // double-timing Spudd loop is to replace both
  // optimalSpudd and approximateSpudd in a single
  // super-efficient piece of code
  // optional argument specifies the value function to 
  // start with
  virtual int Spudd(DdNode *Reward=NULL);
  

  // accessor: returns the number of variables in this MDPs
  int getNumvars();
  int getNumOrigvars();
  int getNumActions();
  void setHorizon(int h);
  void setApproxMeth(int am) {approxMeth = am;};
  DdNode *getNewPrime(int act, int ovar);
  void setNewPrime(int act, int var, DdNode *newNewPrime);
  void flattenT(double ***, int &, int &);
  void decodevarvals(int *varvals, int state);
  DdNode *getOptValue() {return OptimalValue;}
  DdNode *getOptPolicy() {return OptimalPolicy;}
  void setOptValue(DdNode *val);
  
  onum *getOrigVars();
  rnum *getVars();
  char *getActionName(int i);
  DdNode *reward() {return RewardD;}
  int getHorizon() {return (int)horizon;}
  double getDiscount();

  // more accessors needed here...


  // prints the names and values of the original vars in state
  void printMDPSpudd(FILE *fp);
  void printState(int *state, FILE *fp);
  void printAction(int act, FILE *fp);
  void printVarName(int v, FILE *fp);
  void printVars(FILE *fp);
  void printReward(FILE *fp);
  void printNewPrime(FILE *fp);
  void printDiscount(FILE *fp);
  void printTolHor(FILE *fp);
  void writeVPDot(DdNode *valADD, DdNode *polADD,char *fileprefix);

  // simulates the MDP for actionTaken given state
 void simulate(int *state, int actionTaken, int *newstate, int *exclude=NULL);
 // returns reward achieved
 double simulater(int *state, int actionTaken, int *newstate, int *exclude=NULL);
  // simulates the MDP given the optimal action
  void simulate(int *state, int *newstate);

  // updateQ updates the Q functions (Vcurrents) based on current transition function
  // and reward function e.g. Dyna
  void updateQ(DdNode **Q, int acttaken, int *varvals, int k=6);
  
  // replaces dd1 by ddr for assignment of variables in varvals
  void replaceVal(DdNode **dd1, DdNode *ddr, int *varvals);

  // gets value of val and returns in rval
  void getVal(DdManager *dd, DdNode *val, Pair & rval, int *varvals, bool primedvars=false);
  void getVal(DdManager *dd, DdNode *val, Pair & rval, int *varvals, int *varvalsp);

  // normalizes funcition counts by dividing it by itself to yield (normalized) function f
  void normalizeFunction(DdNode *counts, DdNode **f, int j);
  void normalizeFunction(DdNode *counts, DdNode **f, bool primedvars=false);

  // extract a policy from Q functions
  DdNode * extractQPolicy(DdNode **Q, DdNode **ma);

  // gets action with maximal q value
  int getMaxQ(DdNode **Q, int *varvals);

  void classifySupport();
  DdNode *buildCPTCount(int opvar, int acttaken, int *opvarvals, int *ovarvals);

  void printDDSingleVariable(DdNode *dd, double *vv, int ovar);

  // performs variable reordering using reorderMeth according to ADD 'accordingTo'
  // (accordingTo only used for reorderMeth == REORDER_MINSPAN)
  void reorder(DdNode *accordingTo=NULL);

  // add a new ADD variable to gbm
  void newADDVar();
  
  // start policy Querying server
  void startServer();

  int consultOptimalPolicy(DdNode *pol, DdNode *val, int *varvals);

  int consultOptimalPolicy(int *varvals);

	   
  // looks at the OptimalPolicy and computes the 'agreement' between the policies for
  // all the pairs of values of original variable ovar
  // puts the results in mergeStates, which, if n = orig_vars[ovar].nvals
  // should have lenth n*(n-1)/2
  void findLargestSpan(int ovar, double *mergeStates);

  // counts the number of *states* spanned by the add != 0
  // over the novars original variables in listovars
  double countSpan(DdNode *add, int *listovars, int novars, double spanStates);

  // builds a 'stay-same' dd for original variable ovar
  DdNode *buildSameDD(int ovar);
  // builds a cube over orignal variables as in varvals
  // if primed = true - builds it over the primed variables
  // otherwise over the unprimed ones
  DdNode *buildCubeOrig(int *varvals, bool primedvars=false);
  // builds a cube over a single original variable value
  DdNode *buildOneCubeOrig(int ovar, int ovarval, bool primedvars=false);

  // writes policy and value functdion to Dual Optimal file
  // also writes to dot file if dotFlag = true
  void writePolicyToFile(char *, bool dotFlag = false);

  // backs up valD using policy given by actD
  DdNode *evalPolicy(DdNode *valD, DdNode *actD);

  // backs up valD using policy given by actD for numiters
  DdNode *modPolicyIteration(DdNode **valD, DdNode *actD, int numiters);
  DdNode *modifiedPolicyIteration(bool fromOpt=false);
  DdNode *evaluatePolicy(DdNode *actD, DdNode *valD);
  void pickPolicyAction();
  DdNode *improvePolicy(DdNode **valD);

  // randomly shuffle the order
  void shuffleRandom();

  int reinforcementLearn(MDP *simulator, int numiterations);

  // compares the ADDs corresponding to orig variable ovar's branches in the value function
  int compareValueVals(int ovar, double tol, int **compmat, int *nov);
  int compareSubPol(int ovar, int *mergeStates);

  void compareSubPol(int ovar);

  // writes out dot files for valDD and polADD to [fileprefix]value.dot and [fileprefix]action.dot (resp.)
  //void writeVPDot(DdNode *valADD, DdNode *polADD,char *fileprefix);
  
  // writes out dot files for dual action diagrams and reward function
  void writeMDPDot(char *fileprefix);
  void writeDdDot(DdNode *thedd, char *filename);

  // writes Optimal (Approximate) value and policy to file as in writeVPDot
  void writeOptimalVPDot(char *fileprefix);
  void writeApproximateVPDot(char *fileprefix);

  //Routine to store and upload ADDs
  int writeDualOptimal(char *binfilename,DdNode *action,DdNode *value,rnum vars[],rnum prime_vars[],onum orig_vars[]);
  int readDualOptimal(char *filename);
  int readDualOptimal(DdManager *dd,char *binfilename,DdNode ***DualOptimal);

  // for printing out statistics to outpath/infile-stats.dat
  void printStats(DdNode *, DdNode *, FILE*);
  void printOptimalStats(char *);
  void printApproximateStats(char *);

  void printAdd(DdNode *add, FILE *fp);
  void doshit(int *varvals, int sv, bool doq=false);
  void doshit2(int *varvals, int nsv, int *sv);

  // our friends the parsers - these have to be friends because they
  // are filling parts of the MDP pointed to by the global MDP pointer (themdp)
  friend int yyparse();
  friend int parseVars(FILE *,int *, rnum *, rnum *);
  friend int findVar(int *, int, char *, rnum *);
  friend DdNode *parseADD(FILE *,int, rnum *,int);
  friend DdNode *parseDiagram(FILE *errf, int, rnum*, int);
  friend int parseActions(FILE *,int, int *, rnum *, actt *);
  friend int pparse(FILE *,double *, double *, double *, int *, int *, rnum *, rnum *, actt *);
  
};
#endif
