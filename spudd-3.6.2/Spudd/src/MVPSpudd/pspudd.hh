#ifndef _SPUDD
#define _SPUDD

#include <sys/time.h>
#include <sys/times.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "cudd.h"
//#include "cuddObj.hh"
#include "dddmpInt.h"
#include "dddmp.h"
//#define OPTIMIZED

#define MAXLINE 256	/* Maximum length of each input line read. */
#define MAXVARS 64
#define MAXVALS 64
#define MAXLIN 1024
#define MAXACT 256
#define MAXLEN 64
#define DUMPIN 0
#define MAXMEM 700000000
#define MAXMEMHARD 800000000 
#define HREWARD 270
#define HACTION 271
#define HVALUE 272
#define HPOLICY 273
#define THRESHOLD  1.0     /* Maximum error allowed before prunning an ADD */
#define MAX_SIZE  10000
#define MAX_ERROR 0.00001
#define BIGADD 10000
#define LIMIT_ERROR 371
#define ERROR_TOL 0.1
#define LIMIT_SIZE  372
#define SIZE_TOL  2
#define LARGEFLAG 1000
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

/* Global Variables  */
extern DdNode *MergedMin_Add;
extern DdNode *VerySmall, *Zero,*One,*discount,*Half;
extern DdNode **Array1,**Array2,**Array,**ArraySum;
extern DdNode ***ActionsPrimed, ***NewPrime,***Allprime,***Actions;
extern DdNode *RewardD,*MergedMin;
extern DdNode *actionCost[MAXACT];
extern DdNode *ApproximatePolicy,*ApproximateValue, *ApproximatePolicyValue, *ApproximateMidPointValue;
extern DdNode *OptimalValue, *OptimalPolicy;

extern int allPlist[MAXACT][MAXVARS];
extern rnum prime_vars[MAXVARS], vars[MAXVARS];
extern onum orig_vars[MAXVARS];

extern actt actionlist[MAXACT];

extern char **varnames, **actionnames, **lnames;
extern char *outpath,*inpath,*infile;
extern char linebuf[MAXLINE],name[MAXLINE],primedname[MAXLINE];


extern double tolerance,off,discount_factor,horizon;
extern Pair avgError;
extern Pair *pZero;
extern Pair *pOne;
extern Pair*toleranceP;
extern int cntError,toleranceMode;

extern double prune; //pruning strength for the max-error mode
extern int Leaf_Counter;  //counter used in the approximation algorithm.....
extern int numorigvars,numvars,numactions;
extern double Max_Error;
extern int Max_Size;
extern int bigadd;
extern bool mvinput;
extern double extR;
extern double SUMSQUARE;
extern int numSumSquare;
/* Global BDD manager. */
//extern Cudd gbm;
extern DdManager *gbm; 

/* routines */
DdNode *Cudd_Else(DdNode *);
DdNode *Cudd_Then(DdNode *);
DdNode *sumSubtrees(DdNode *, int, int,  rnum *, int);
DdNode  *primeReward(DdNode *, DdNode **, rnum *, int, int *, int, int);
DdNode  *newPrimeReward(DdNode *, DdNode **);
DdNode *getActionAdd(DdNode **, DdNode *, int );
void getActionAdds(DdNode **, DdNode *, int , DdNode **, double);
DdNode *evalPolicy(DdNode *valD, DdNode *actD);

DdNode *Aspudd_Proc(int);
DdNode *NewAspudd_Proc(int,int,int,FILE*,FILE*);
DdNode *spudd_Proc(int,FILE*);
DdNode *addRoundOff(DdNode *,double *);
DdNode *myRoundOff(DdManager *, DdNode **, DdNode **);
DdNode *roundOffMMA(DdNode *);
DdNode *approxDd(DdManager *, DdNode **, DdNode **);
DdNode *getTotalSpan(DdManager *, DdNode **, DdNode **);
DdNode *addAddExact(DdManager *dd, DdNode **, DdNode **);
DdNode * My_cuddAddApplyRecur(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode *My_addApply(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode * MyPlus(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *getErrorAdd(DdManager *dd, DdNode **f, DdNode **x);
DdNode *My_addMaximum(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode * addMean(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *SumSquareError(DdManager * gbm, DdNode ** f, DdNode ** g);



void writeVPDot(DdNode *valADD, DdNode *polADD,char *fileprefix);


void pQuery(DdManager *dd, DdNode *value, DdNode *action, char * pfile, onum *orig_vars, rnum *vars, int, int);

void policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars);
void getAV(DdManager *dd, DdNode *val, DdNode *act, Pair & dval, Pair & aval, int *varvals);

//Routine to store and upload ADDs
int fillRnumDdNode(DdManager *dd, rnum *variables,rnum *prime_variables);
int compareRnumStruct(int num,rnum varsP[],rnum varsT[]);
int compareOnumStruct(onum orig_varsP[],onum orig_varsT[]);
int structCopy(rnum *varsIn,rnum *primeVarsIn,onum *origVarsIn);
int writeDualOptimal(char *binfilename,DdNode *action,DdNode *value,rnum vars[],rnum prime_vars[],onum orig_vars[]);
int readDualOptimal(DdManager *dd,char *binfilename,DdNode ***DualOptimal);




// parser routines
void addVar(onum *, int, char*);
void newADDVar(onum *ovar, int novar);
int findOVar(const char *findname);
int findOVal(int ovar, const char *findval);
void setActionCost(int nac, double val);
//DdNode *buildCPTHelp(rnum *, int *currv, int firstv, int lastv, DdNode **cpt, int cptsize);
//DdNode * buildCPT(rnum *, int cpov, DdNode **cpt);
DdNode * buildCubeCPT(rnum *v, int coval, onum ovar, DdNode *add);
void removeDummy(DdNode **np, int);
void removeAllDummys(DdNode **np);
void buildGoodStateADDs();

void shuffleRandom();
double  get_extent(DdNode *);
double  get_error(DdNode *);
Pair get_span(DdNode *);
int close_enough(double , double , double ); 
char* DumpDoth(DdManager *dd, DdNode *add, onum *orig_vars, int rovar,char **,  FILE *fp, FILE *nfp, int *hmn);
int DumpDot(DdManager *dd, DdNode *add, rnum *vars, onum *orig_vars,char **, FILE *fp);
int My_DumpDot(DdManager *,int  n , DdNode **, char **,char **,char **, FILE *);
int aconvert(char *, char **, double, char *separator="\\n");
void print_out(DdNode *, int, rnum *, DdNode *, DdNode *, int);
int count_internal_nodes(DdNode *);
int count_leaves(DdNode *);
int max(int, int);
void dumpHtmlRew(FILE *log,char *outpath);
int reorderMinSpan(DdNode *);
void infoSorter(double *infoList, int *reorderList);
void removeElement(list_cst *list,list_cst *element);
DdNode *replacePairs(DdManager *dd, DdNode **f, DdNode **newC);
void allPairsApprox();
void testReorder(DdNode *, int, FILE*, FILE*);
/*parser routines */
int parseVars(FILE *,int *, rnum *, rnum *);
int findVar(int *, int, char *, rnum *);
DdNode *parseADD(FILE *,int, rnum *,int);
DdNode *parseDiagram(FILE *errf, int, rnum*, int);
int parseActions(FILE *,int, int *, rnum *, actt *);
int pparse(FILE *,double *, double *, double *, int *, int *, rnum *, rnum *, actt *);


#endif
