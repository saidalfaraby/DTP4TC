#ifndef __POMDP
#define __POMDP
#include "MDP.h"
#define MAXNUMBELIEFS 200
#define MAXMGVAL 270
#define MINMGVAL 0
#define MAXMGPVAL 1.8
#define MINMGPVAL 0.0
struct image_ {
  int nx;
  int ny;
  int *pixels;
  double *vals;
  int *action;
};

typedef struct image_ image;

class Sample {
 public:
  int index;
  double val;
  double *z;
  Sample() {
    index = 0;
    val = 0;
  };
  Sample (int n) {
    index = 0;
    val = 0;
    z = new double[n];
  };
};


class POMDP : public MDP {
 public:
  // number of original observation functions
  int numorigobs;

  // observation variable structure
  onum orig_obs[MAXVARS];

  // observation functions  one for each action-observation variable pair
  DdNode ***totalObsFun;
  
  // PDFs for samples
  DdNode ***obsSampleProb;
  MixGauss *****proposalsamples;
  DdNode ****pojk;
  DdNode ****fajk;

  // initial, complete belief
  DdNode *initBelief;

  //factored belief
  DdNode **initBeliefState;

  // forward propagate belief state
  DdNode ****baz;
  DdNode ***ba;
  // solution stuff
  int *bestactions;
  DdNode **alphas;
  int numalphas;

  DdNode ***beliefs;
  DdNode **ufbeliefs;
  int numbeliefs;

  bool dosampling;

  // constructor 
  POMDP() : MDP() {
    numorigobs = 0;
  }
  // constructor with parser read
  POMDP(char *infile, double badd=BIGADD);
  POMDP(char *infile, bool ds = false, double badd=BIGADD);
  void allocate2Dstructures(int nsamples);

  void getMostLikelyBelief(int *varvals);
  DdNode* getInitBelief();
  DdNode* getInitBeliefPrimed();
  void newADDObs();
  int generatePolicy(int nsamples, int mnb, char *, char *beliefFile=NULL);
  
  double computeHorizon(double);
  void getInitialAlphas(int typ);
  int *valueIteration(int horizon, double *beliefValues, char *);
  int * fastPBVI(int horizon, char *);
  int * fastPBVISampled(int horizon, char *);
  int * fastPBVISampled2(int horizon, char *);

  void partitionObservationSpace(int type, MixGauss **mgs, MixGauss **regions);
  void partitionObservationSpace(MixGauss **mgs, int & snr, double * & bounds, double* &maxmgs);


  DdNode *partitionIntegrateObsSpace(int pbelief, int a);

  void computeIntegrateRegions(DdNode **obsfun, DdNode **intObsFuns, int a, int b);
  DdNode *sampleBackup(int pbelief_ind, int a);
  DdNode *sampleBackup2(int pbelief_ind, int a, char *outpath);

  DdNode *backupAlphas(DdNode **obsfun, DdNode **transition);
  DdNode *regionIntegrateBackup(int b, int a);


  int *getPolicy(int *bestactions, int numbeliefs, DdNode ***beliefs, double *beliefValues);
  int getPolicy(int *bestactions, DdNode **belief, double & beliefValue);
  int getPolicy(int *bestactions, DdNode *belief, double & beliefValue);
  int getPolicy(DdNode **belief, double & beliefValue);
  int getPolicy(DdNode *belief, double & beliefValue);

  double computeValue(int & maxi, DdNode **belief);
  double computeValue(int & maxi, DdNode *belief);
  double computeValue(int numalphas, int & maxi, DdNode **alphas, DdNode **belief);
  double computeValue(int numalphs, int & maxi, DdNode **alphs, DdNode *belief);

  DdNode *computeValue(DdNode *alpha, DdNode **belief);
  DdNode *computeValue(DdNode *alpha, DdNode *belief);

  void expandBeliefs(double horizon, double mthresh, bool allbeliefs = true);
  void expandBeliefs(int mnumbeliefs, double mthresh, bool allbeliefs = true);
  void expandBeliefs(double mthresh, bool allbeliefs = true);

  DdNode *beliefProduct(DdNode **b);
  DdNode **beliefMarginal(DdNode *b);
  double supremumNormMar(DdNode **, DdNode **);

  int backup(DdNode **palpha, double *maxgamma, int pbelief_ind, DdNode **pbelief, char *outpath);
  
  void printAlphas();
  void printBeliefs(DdNode **,int);
  DdNode *simulateAction(DdNode **bb, int a);
  DdNode *simulateAction(DdNode **bb, int a, double **z);
  double simulator(DdNode **initb, int stages, char *outpath, bool verbose = false, bool funky=false);
  double simulator(int stages, char *outpath, bool verbose = false, bool funky = false);
  DdNode *simulateTransition(DdNode *stsamp, int a);

  DdNode * sample_Belief(DdNode *b);
  DdNode * sample_Belief(DdNode **b);
  void drawSampleFromBelief(DdNode *newb, int a, MixGauss **v) ;

  DdNode * getBelief(int ovar, double *ovarvals, bool primedvars);

  void addToReward(double addval);

  DdNode * bayesianUpdate(DdNode **b, int a);
  DdNode * bayesianUpdate(DdNode *belief, int a);

  DdNode * bayesianUpdate(DdNode **b, int a, DdNode *obsprod);
  DdNode * bayesianUpdate(DdNode *b, int a, DdNode *obsprod);
  DdNode * multBeliefObsProd(DdNode *belief, DdNode *obsprod);
  DdNode * computeSampleObsProd(DdNode *stsamp, int a, bool verbose = false);
  DdNode * computeSampleObsProdFunky(DdNode *stsamp, int a, bool verbose = false);

  double L2dist(DdNode *, DdNode *);

  void normalizeReduceFunction(DdNode *counts, DdNode **f, int n, bool primedvars);
  double checkNormalization(DdNode *b);
  void flattenT(double ***T, int & ns, int & na);

  DdNode *getSpan(DdNode *f, double extra, int primedvars);

  void getEdgeImage2D(int a, int o, char *edgeimage) ;

  void discretizeObservationFunction(int a, int o, image & oimage);

  MixGauss evaluateDdAtCube(DdNode *thedd, DdNode *thecube, bool primed);



  void printPolicy(int *policy, int numbeliefs, DdNode ***beliefs, double *beliefValues);
  void printPolicyHTML(char *fname, int *policy, int numbeliefs, DdNode ***beliefs, double *beliefValues, int hor, int ov);
  void printConditionalPlanRegions(char *fprefix, int *policy, int numbeliefs, DdNode ***beliefs, int o);
  void printConditionalPlanRegions2(char *fprefix, int *policy, int numbeliefs, DdNode ***beliefs, int o);
  void printConditionalPlanRegions2(char *fprefix, int *policy, int numbeliefs, int *pb, DdNode ***beliefs, int o);
  void printConditionalPlanRegions3(char *fprefix, int acttotake, DdNode *belief, int o);
  void printConditionalPlansForBeliefs(int *pb, int numpb, char *outpath);

  DdNode * computeObsProb(DdNode *lbaz, double *z);
  void printObsImage(char *fname, char *imname, int o, int a);

  void printAlphasFile(char *fname);
  void readAlphasFile(char *fname);
  void printBeliefsFile(char *fname);
  void readBeliefsFile(char *fname);
  void printBelief(DdNode **b, FILE *fd=stdout);
  void printBelief(DdNode *b, FILE *fd=stdout);
  void printBeliefComparetoActual(FILE *fid, DdNode *b, int *varvals);
  void writeObsFunction2D(int a);

  friend int yyparse();
  
};
void writeImagePPM(char *filename, unsigned char *im, int nx, int ny);
bool converged(double *Vb, double *oVb, int numb, int counter, int horizon, double tolerance);
DdNode *My_addReducePrecision(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *My_addSample(DdManager * dd, DdNode ** f, DdNode ** g);
double supremumNorm(DdNode *dd1, DdNode *dd2);
int member(DdNode *belief, DdNode **beliefSet, int numbeliefs, double thresh);
DdNode *My_addAbsMinus(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *redistributeExtra(DdManager * dd, DdNode ** f, DdNode ** g);
DdNode *combineObsFun(DdManager *dd, DdNode **f, DdNode **g);
DdNode *computePDF(DdManager *dd, DdNode **f, DdNode **g);
DdNode *computePDFType0(DdManager *dd, DdNode **f, DdNode **g);
DdNode *extractGaussians(DdManager *dd, DdNode **f, DdNode **g);
void writeObsImage(int ind, double *mgval, char *filename, char *eimage);
void braindeadEdgeDetect(int nx, int ny, int *iim, char *oim);
void writeEdgeImagePGM(char *filename, char *eimage);
void writeEdgeImagePPM(char *filename, char *imname, char *eimage);
void writeObsImageRaw(char *filename);
void writeObsImageMgs(char *filename);
DdNode *ddintegrate(DdManager *dd, DdNode **f, DdNode **g);
double integrate(MixGauss & integrand, MixGauss & region);
DdNode *drawSample(DdManager *dd, DdNode **f, DdNode **g);

void extractPrintSamples(FILE *fd, DdNode *dd);



void merge(const double *zc, const int nr, double **szc, int & snr);
#endif
