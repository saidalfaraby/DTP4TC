#ifndef __MDPPSPUDD
#define __MDPPSPUDD 
#include <sys/time.h>
#include <sys/times.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "../../../CuddPP/include/cudd.h"
#include "../../../CuddPP/include/cuddObj.hh"
#include "../../../CuddPP/include/dddmpInt.h"
#include "../../../CuddPP/include/dddmp.h"
#define MAXLINE 256	/* Maximum length of each input line read. */
#define MAXVARS 64
#define MAXVALS 64
#define MAXLIN 1024
#define MAXACT 256
#define MAXLEN 64
#define DUMPIN 0
#define MAXMEM 700000000
#define MAXMEMHARD 800000000 
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

/* routines */
DdNode *Cudd_Else(DdNode *);
DdNode *Cudd_Then(DdNode *);
DdNode *sumSubtrees(DdNode *, int, int,  rnum *, int);
DdNode  *primeReward(DdNode *, DdNode **, rnum *, int, int *, int, int);
DdNode  *newPrimeReward(DdNode *, DdNode **);
DdNode *getActionAdd(DdNode **, DdNode *, int );
void getActionAdds(DdNode **, DdNode *, int , DdNode **, double);
DdNode *evalPolicy(DdNode *valD, DdNode *actD);

DdNode *myRoundOff(DdManager *, DdNode **, DdNode **);
DdNode *getTotalSpan(DdManager *, DdNode **, DdNode **);
DdNode *addAddExact(DdManager *dd, DdNode **, DdNode **);
DdNode * My_cuddAddApplyRecur(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode *My_addApply(DdNode * (*op)(DdManager *, DdNode **, DdNode **), DdNode * , DdNode * );
DdNode * MyPlus(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *getErrorAdd(DdManager *dd, DdNode **f, DdNode **x);
DdNode *My_addMaximum(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode * addMean(DdManager * dd, DdNode ** f, DdNode ** g);



void pQuery(DdManager *dd, DdNode *value, DdNode *action, char * pfile, onum *orig_vars, rnum *vars, int, int);

void policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars);

double  get_extent(DdNode *);
double  get_error(DdNode *);
Pair get_span(DdNode *);
int close_enough(double , double , double ); 
char* DumpDoth(DdManager *dd, DdNode *add, rnum *vars, int nv, onum *orig_vars, int nov, int rovar, char **,  FILE *fp, FILE *nfp, int *hmn);
int DumpDot(DdManager *dd, DdNode *add, rnum *vars, int nv, onum *orig_vars, int nov, char **, FILE *fp);
int My_DumpDot(DdManager *,int  n , DdNode **, char **,char **,char **, FILE *);

int aconvert(char *, char **, double);
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
int findOVal(onum *ovar, int ovar, const char *findval);
void setActionCost(DdNode **ac, int nac, double val);
DdNode * buildCubeCPT(rnum *v, int coval, onum ovar, DdNode *add);
void removeDummy(DdNode **np, int);
void removeAllDummys(DdNode **np, int);
void buildGoodStateADDs(onum *, rnum *, int);

bool member(int *, int, int);
#endif

