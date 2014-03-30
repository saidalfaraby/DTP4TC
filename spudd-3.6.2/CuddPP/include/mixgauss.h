#ifndef _MIXGAUSS
#define _MIXGAUSS
#include "gaussian.h"

//#define USE_MATHEMATICA
#ifdef USE_MATHEMATICA
#include "mathlink.h"
#endif
// a mixture of Gaussians is a linear combination of Gaussians with real-valued weights
// plus a constant. This constant is there to make it so we can also use the mixture of Gaussians 
// as a double, multiply a mixture of Gaussians by a double by multiplying two MixGauss, etc.
// if nmix = 0, then the MixGauss is just a double.
// if fvdim = 0, then the MixGauss is just an array of weights (the mixweights) + the double if needed.
//   this corresponds to a simple discrete probability function
class MixGauss {
 public:
  // initialises a 0-mixture Gaussian n the mixture with weight 1
  MixGauss();
  // 0-mixture gaussian with a value
  MixGauss(double);
  // intialises a mixture with nmix components, each with dimension fvdim
  // initial weights are all 1
  MixGauss(int fvdim, int nmix);
  // intializes all types of MixGauss
  MixGauss(int type, int numparams, double *);
  MixGauss(const MixGauss & mg);
  
  // destructor
  ~MixGauss();

  void copyFrom(const MixGauss & mg, int m=-1);
  void copyFrom(const MixGauss * mg);
  void getMem();
  void getMemNoGauss();
  void delMem();

  void setValue(double);
  void set(int type, int numparams, double *);
  void set(double);
  void set(int);

  void set_min(double);
  void set_max(double);
  void set_weights(double *);
  void setGaussians(double *, double **, double ***);
  void setGaussian(int, double, double *, double **);
  void setGaussian(int, double, double *, double *);

  double get_val();
  double get_min();
  double get_max();
  double pdf(double *obs) const;
  double pdf(MixGauss & mg) const;
  double pdf(int obs) const;
  double * drawSample() const;
  int drawSample(double * &) const;
  int drawMixtureSample() const;

  double integrate(int typ, double lo, double hi) const;

  int intersections(MixGauss & mg, double **zc);
  int oldintersections(MixGauss & mg, double **zc);

  MixGauss operator+(const MixGauss&) const;
  MixGauss operator*(const MixGauss&) const;
  MixGauss operator*(double) const;
  MixGauss operator/(MixGauss&) const;
  MixGauss operator-(const MixGauss&) const;
  MixGauss operator/(double) const;
  MixGauss operator-() const;
  MixGauss addorsubtract(const MixGauss & mg1, const MixGauss & mg2, bool) const;
  MixGauss multiply(const MixGauss & mg1, const MixGauss & mg2) const;

  bool operator>=(const MixGauss&);
  bool operator<=(const MixGauss&);
  bool operator>(const MixGauss&);
  bool operator<(const MixGauss&);
  bool operator!=(const MixGauss&); 
  bool operator==(const MixGauss&);

  void operator=(const MixGauss &);

  /*
  bool operator<(int);
  */

  void removeZeros();

  MixGauss absVal() const;

  double hashCode() const;

  char *toString() const;
  int readString(FILE *fp, char * buf);
  void parseString(char *buf);

  void print(FILE *) const;
  void printValMixWeights(FILE *fd) const;
  // variables
  // number of mixture components
  int nmix, fvdim;
  // and a constant
  double val;
  // array of pointers to Gaussians
  Gaussian **mixcomps;
  // array of weights
  double *mixweights;
  bool gotMem, gotMemNoGauss;

  friend class Gaussian;
};
MixGauss ceill(const MixGauss& source);
double erfc (double x);
double exparg ( int *l );
int ipmpar ( int *i );
void findroot(MixGauss & mg, double lo, double hi, double acc, double * theroots, int &numroots);
int ftbis(MixGauss  & mg, double a, double b, double acc, double & result);

#ifdef USE_MATHEMATICA
int getroots_mathematica(char *inputstring, double **theroots);
static void init_and_openlink( int argc, char* argv[]);
void initialize_mlink();
void terminate_mlink();
static void error( MLINK lp);
#endif



typedef MixGauss* PtrTerminal;
typedef MixGauss Terminal;
#endif
