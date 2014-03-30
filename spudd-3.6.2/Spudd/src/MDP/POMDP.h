#include "MDP.h"
#define MAXNUMBELIEFS 200
class POMDP : public MDP {
 protected:
  // number of original observations
  int numorigobs;

  // total number of possible observation combinations
  int totalNumObs;

  // number of boolean observation ADD variables
  int numobs;

  // observation variable structure
  rnum obs[MAXVARS];
  onum orig_obs[MAXVARS];

  // factored observation function
  DdNode ***ObsFun;

  // product of all observation functions
  DdNode **totalObsFun;

  // initial, complete belief
  DdNode *initBelief;

  //factored belief
  DdNode **beliefState;

  // solution stuff
  int *bestactions;
  DdNode **alphas;
  int numalphas;
 public:
  // constructor 
  POMDP() : MDP() {
    numobs = 0;
    numorigobs = 0;
    totalNumObs = 0;
  }
  // constructor with parser read
  POMDP(char *infile, double badd=BIGADD);
  void getMostLikelyBelief(int *varvals);
  void computeTotalNumObs();
  void newADDObs();
  int generatePolicy(int rMeth = REORDER_NONE,
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
  int generatePolicy(bool factored_beliefs);
  
  DdNode * getBelief(int ovar, double *ovarvals, bool primedvars);
  double computeHorizon(double);
  void getInitialAlphas(DdNode **alphas, int numalphas, int typ);
  int *valueIteration(int horizon, int & numbeliefs, DdNode ***beliefs, bool factbelief);
  int * fastPBVI(int horizon, int & numbeliefs, DdNode ***beliefs, int & numalphas, DdNode ***palphas, bool factbelief);
  
  int getPolicy(int numalphas, DdNode **alphas, int *bestactions, DdNode *belief);
  int *getPolicy(int numalphas, DdNode **alphas, int *bestactions, int numbeliefs, DdNode **beliefs);
  int getPolicy(int numalphas, DdNode **alphas, int *bestactions, DdNode **belief);
  int *getPolicy(int numalphas, DdNode **alphas, int *bestactions, int numbeliefs, DdNode ***beliefs);

  double computeValue(int numalphas, int * maxi, DdNode **alphas, DdNode *belief);
  double computeValue(int numalphas, int * maxi, DdNode **alphas, DdNode **belief);
  DdNode *computeValue(DdNode *alpha, DdNode **belief);

  int expandBeliefs(int numbeliefs, DdNode ***pbeliefs, double horizon, double mthresh, bool allbeliefs = true);
  int expandBeliefs(int numbeliefs, DdNode ***beliefs, double mthresh, bool allbeliefs = true);

  int expandBeliefs(int numbeliefs, DdNode ****pbeliefs, double horizon, double mthresh, bool allbeliefs = true);
  int expandBeliefs(int numbeliefs, DdNode ****beliefs, double mthresh, bool allbeliefs = true);
  DdNode **beliefMarginal(DdNode *b);
  double supremumNormMar(DdNode **, DdNode **);

  int backupValue(int *policy, int numalphas, DdNode ***alpha, int numbeliefs,  DdNode **beliefs);
  DdNode *** step1Backup(int numalphas, DdNode **alpha);
  int step23Backup(int numalphas, DdNode **palpha, double *maxgamma, DdNode ***gamma_ao, DdNode *pbelief);
  int step23Backup(int numalphas, DdNode **palpha, double *maxgamma, DdNode ***gamma_ao, DdNode **pbelief);
  
  void printAlphas(DdNode **,int);
  void printBeliefs(DdNode **,int);
  DdNode *simulateAction(DdNode *b, int a);
  DdNode *simulateAction(DdNode **bb, int a);
  DdNode *simulateAction(DdNode *bb, int a, DdNode *obsprob);

  void simulate(double tol, DdNode *belief, int stages);
  void simulate(double tol, int stages);
  double simulateHandwashing(char *filename, char *ofile, int nub, char *lfile=NULL);


  DdNode * sampleObservation(DdNode *ofun, int o);
  DdNode * sample_Observation(DdNode *b, int o);

  DdNode * sampleBelief(DdNode *b);
  DdNode * sample_Belief(DdNode *b);
  DdNode * sample_Belief(DdNode **b);

  DdNode * bayesianUpdate(DdNode **b, int a, DdNode *obsprod);

  double L2dist(DdNode *, DdNode *);

  DdNode * buildOneCubeOrigObs(int ovar, int ovarval);
  DdNode* buildCubeOrigObs(int *obsvals);
  double evalOrigObs(DdNode *theadd, int ovar, int ovarval);

  void normalizeReduceFunction(DdNode *counts, DdNode **f, int n, bool primedvars);
  double checkNormalization(DdNode *b);
  void flattenT(double ***T, int & ns, int & na);

  DdNode *getSpan(DdNode *f, double extra, int primedvars);


  void printPolicy(int *policy, int numbeliefs, DdNode **beliefs);
  void printPolicy(int *policy, int numbeliefs, DdNode ***beliefs);
  void printAlphas(int numalphas, DdNode **alpha);
  void printBelief(DdNode *b);
  void printBelief(DdNode **b);
  void printBeliefComparetoActual(FILE *fid, DdNode *b, int *varvals);

  friend int yyparse();
  
};
bool converged(double *Vb, double *oVb, int numb, int counter, int horizon, double tolerance);
DdNode *My_addReducePrecision(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *My_addSample(DdManager * dd, DdNode ** f, DdNode ** g);
double supremumNorm(DdNode *dd1, DdNode *dd2);
int member(DdNode *belief, DdNode **beliefSet, int numbeliefs, double thresh);
DdNode *My_addAbsMinus(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *redistributeExtra(DdManager * dd, DdNode ** f, DdNode ** g);
