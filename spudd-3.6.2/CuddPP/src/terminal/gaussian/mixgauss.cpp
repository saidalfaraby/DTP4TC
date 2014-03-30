#include "mixgauss.h"

#ifdef USE_MATHEMATICA
MLENV ep = (MLENV)0;
MLINK lp = (MLINK)0;
#endif

MixGauss::MixGauss()
{
  nmix = 0;
  fvdim = 0;
  val = 0;
  gotMem = false;
  gotMemNoGauss = false;
}
MixGauss::MixGauss(double v)
{
  nmix = 0;
  fvdim = 0;
  val = v;
  gotMem = false;
  gotMemNoGauss = false;
}
// constructs a new MixGauss of type 
// 0: fvdim = 0 is a discrete multinomial density in the mixweights
// N>0: fvdim = N is a N-D covariance Gaussian mixture
//       now the numparams should be either fvdim (type>0, in which case there is one mixture component)
//       or M*(fvdim+1) for M mixture components
MixGauss::MixGauss(int type, int numparams, double *v)
{
  set(type,numparams,v);
}

MixGauss::MixGauss(const MixGauss & mg)
{
  gotMem = false;
  gotMemNoGauss = false;
  copyFrom(mg);
}
MixGauss::~MixGauss()
{
  if (gotMem || gotMemNoGauss)
    delMem();
}
void MixGauss::copyFrom(const MixGauss *mg)
{
  copyFrom(*mg);
}
// copies the first m components of mg over to this
// if m < 0 [default], copies all components
void MixGauss::copyFrom(const MixGauss & mg, int m)
{
  if (m<0)
    m = mg.nmix;
  if (gotMem || gotMemNoGauss)
    delMem();
  val = mg.val;
  nmix = m;
  fvdim = mg.fvdim;
  getMemNoGauss();
  for (int k=0; k<nmix; k++) {
    mixweights[k] = mg.mixweights[k];
    // now, just set the pointers!
    if (fvdim > 0) 
      mixcomps[k] = mg.mixcomps[k];
  }
}
void MixGauss::getMem() {
  int k;
  if (nmix > 0) {
    mixweights = new double[nmix];
    for (k=0; k<nmix; k++) 
      mixweights[k] = 0.0;
    if (fvdim >  0) {
      mixcomps = new Gaussian*[nmix];
      for (k=0; k<nmix; k++) {
	mixcomps[k] =  new Gaussian(fvdim);
	mixcomps[k]->zero();
      }
    }
    gotMem = true;
    gotMemNoGauss = false;
  } else {
    gotMem = false;
    gotMemNoGauss = false;
  }
}
void MixGauss::getMemNoGauss() {
  int k;
  if (nmix > 0) {
    mixweights = new double[nmix];
    for (k=0; k<nmix; k++) 
      mixweights[k] = 0.0;
    if (fvdim >  0) {
      mixcomps = new Gaussian*[nmix];
    }
    gotMemNoGauss = true;
    gotMem = false;
  } else {
    gotMemNoGauss = false;
    gotMem = false;
  }
}
void MixGauss::delMem() 
{
  if (nmix > 0) {
    delete [] mixweights;
    if (fvdim > 0) {
      /*
      if (gotMem) {
	for (int k=0; k<nmix; k++) 
	  delete mixcomps[k];
      }
      */
      delete [] mixcomps;
    }
  }
  gotMem = false;
  gotMemNoGauss = false;
}

void MixGauss::setValue(double v)
{ 
  val = v;
}
void MixGauss::set(int type, int numparams, double *v)
{
  int i,k;
  fvdim = type;
  if (type == 0) {
    nmix = numparams;
    val = 0;
    getMem();
    if (v != NULL) 
      set_weights(v);
  } else if (type >= 1) {
    val = 0.0;
    // figure out nmix
    int paramspermixcomp = 1+fvdim+fvdim*(fvdim+1)/2;
    if (numparams < paramspermixcomp) 
      nmix = 1;
    else 
      nmix = numparams/paramspermixcomp;
    getMem();
    i=0;
    double *w = v;
    if (nmix == 1) {
      setGaussian(0,1.0,w,w+fvdim);
    } else {
      for (k=0; k<nmix; k++) {
	setGaussian(k,*w,w+1,w+1+fvdim);
	w += paramspermixcomp;
      }
    }
  }
}

void MixGauss::set(double v)
{
  val = v;
}
void MixGauss::set(int v)
{
  val = (double) v;
}

void MixGauss::set_min(double v)
{
  val = v;
}
void MixGauss::set_max(double v)
{
  val = v;
}
void MixGauss::set_weights(double *weight)
{
  for (int i=0; i<nmix; i++) 
    mixweights[i] = weight[i];
}
void MixGauss::setGaussians(double *weight, double **means, double ***cvs)
{
  for (int i=0; i<nmix; i++) {
    setGaussian(i,weight[i],means[i],cvs[i]);
  }
}
void MixGauss::setGaussian(int index, double weight, double *mean, double *cv)
{
  if (index < nmix) {
    mixweights[index] = weight;
    if (fvdim > 0) {
      mixcomps[index]->setMean(mean);
      mixcomps[index]->setCovariance(cv);
    }
  }
}
void MixGauss::setGaussian(int index, double weight, double *mean, double **cv)
{
  if (index < nmix) {
    mixweights[index] = weight;
    if (fvdim > 0) {
      mixcomps[index]->setMean(mean);
      mixcomps[index]->setCovariance(cv);
    }
  }
}
double MixGauss::pdf(MixGauss & mg) const
{
  if (fvdim == 0) {
    return pdf((int) (mg.val));
  } else {
    return pdf(mg.mixweights);
  }
}
double MixGauss::pdf(double *obs) const
{
  int k;
  double sum = val;
  if (fvdim > 0) {
    for (k=0; k<nmix; k++) {
      sum += mixweights[k]*mixcomps[k]->pdf(obs);
    }
  }
  return sum;
} 
// only for fvdim = 0;
double MixGauss::pdf(int obs) const
{
  double res = val;
  if (obs <nmix) 
    res += mixweights[obs];
  return res;
}
//only when fvdim = 0
// returns -1 if there are no mixture components
int MixGauss::drawMixtureSample() const
{
  int k;
  double rnd = ((double) rand())/((double) RAND_MAX+1.0);
  double sum = 0;
  k=0;
  while (sum < rnd && k < nmix) 
    sum += mixweights[k++];
  k--;
  return k;
}
double *MixGauss::drawSample() const
{
  int k,i;
  // draw a sample from each Gaussian
  if (nmix == 0 || fvdim == 0) {
    return NULL;
  }
  double *thesample= new double[fvdim];
  drawSample(thesample);
  return thesample;
}
int MixGauss::drawSample(double * &thesample) const
{
  int k,i;
  // draw a sample from each Gaussian
  if (nmix == 0 || fvdim == 0) {
    return 1;
  }
  // first, sample a mixture component
  // this assumes that the sampling weights are normalised!
  k = drawMixtureSample();
  if (k < 0) {
    fprintf(stderr,"---------------------- sampling failed -mixture weights not normalized\n");
    exit(0);
  }
  // draw from the kth mixture component
  thesample = mixcomps[k]->drawSample();
  return 1;
}
#ifdef USE_MATHEMATICA
int MixGauss::intersections(MixGauss & mg, double **zc) {
  // compute this-mg
  int i,j,minj,k;
  if ((*this) == mg) {
    return 0;
  }
  MixGauss tmpmg((*this)-mg);

  double lobound(1e+300), hibound(-1e+300), tmp;
  for (k=0; k<tmpmg.nmix; k++) {
    if ((tmp = tmpmg.mixcomps[k]->u[0]-tmpmg.mixcomps[k]->U[0][0]*10) < lobound) 
      lobound = tmp;
    if ((tmp = tmpmg.mixcomps[k]->u[0]+tmpmg.mixcomps[k]->U[0][0]*10) > hibound) 
      hibound = tmp;
  }

  // construct the string to evaluate
  char inputstring[4098];

  sprintf(inputstring,"IntervalBisection[");
  for (k=0; k<tmpmg.nmix; k++) {
    if (tmpmg.mixcomps[k]->u[0] < 0) {
      sprintf(inputstring,"%s+%.30fExp[-0.5(1.0/%.30f)(x+%.30f)^2]",inputstring,tmpmg.mixweights[k],tmpmg.mixcomps[k]->U[0][0],fabs(tmpmg.mixcomps[k]->u[0]));
    } else {
      sprintf(inputstring,"%s+%.30fExp[-0.5(1.0/%.30f)(x-%.30f)^2]",inputstring,tmpmg.mixweights[k],tmpmg.mixcomps[k]->U[0][0],tmpmg.mixcomps[k]->u[0]);
    }
  }
  double accuracy(0.01);
  sprintf(inputstring,"%s, x, Interval[{%f, %f}], %f, MaxRecursion -> 100]%c",inputstring,lobound,hibound,accuracy,10);
  //fprintf(stderr,"%s\n",tmpmg.toString());
  //fprintf(stderr,"inputstring is %s",inputstring);
  
  int nroots = getroots_mathematica(inputstring, zc);
  /*
  fprintf(stderr,"num roots found %d: ",nroots);
  for (k=0; k<nroots; k++) {
    fprintf(stderr,"%g ",(*zc)[k]);
  }
  fprintf(stderr,"\n");
  */
  return nroots;
}
#endif
// outputs an ordered  vector *zc with all the places where this MixGauss crosses mg
// returns the number of such crossings
// only for 1D functions
int MixGauss::oldintersections(MixGauss & mg, double **zc) {
  // compute this-mg
  int i,j,minj,k;
  if ((*this) == mg) {
    return 0;
  }
  MixGauss tmpmg((*this)-mg);
  double minmn, largestmn;
  double *mns = new double[tmpmg.nmix];
  int *mnind = new int[tmpmg.nmix];
  bool *picked = new bool[tmpmg.nmix];

  for (j=0; j<tmpmg.nmix; j++) {
    mnind[j] = -1;
    picked[j] = false;
  }
  largestmn = tmpmg.mixcomps[0]->u[0];
  // find largest mean
  for (j=1; j<tmpmg.nmix; j++) {
    if (tmpmg.mixcomps[j]->u[0] > largestmn) 
      largestmn = tmpmg.mixcomps[j]->u[0];
  }
  // sort the means
  for (k=0; k<tmpmg.nmix; k++) {
    minmn = largestmn;
    minj = 0;
    for (j=0; j<tmpmg.nmix; j++) {
      if (!picked[j] && tmpmg.mixcomps[j]->u[0] < minmn) {
	minmn = tmpmg.mixcomps[j]->u[0];
	minj = j;
      }
    }
    picked[minj] = true;
    mnind[k] = minj;
    mns[k] = minmn;
  }
  // find largest and smallest bounds on this 
  double lobound(1e+300), hibound(-1e+300), tmp;
  for (k=0; k<tmpmg.nmix; k++) {
    if ((tmp = tmpmg.mixcomps[k]->u[0]-tmpmg.mixcomps[k]->U[0][0]*3) < lobound) 
      lobound = tmp;
    if ((tmp = tmpmg.mixcomps[k]->u[0]+tmpmg.mixcomps[k]->U[0][0]*3) > hibound) 
      hibound = tmp;
  }
  double acc(0.01);
  double a(lobound),b;
  int foundone;
  int numroots(0), nnumroots(0);

  k = 0;

  // step through the function from mean to mean and find roots
  b = mns[0];
  *zc = new double[256];
  while (k <= tmpmg.nmix) {
    // find a root in [a,b]
    nnumroots = 0;
    // at end points, there can be only one root left
    findroot(tmpmg,a,b,acc, (*zc)+numroots,nnumroots);
    if (nnumroots > 0) {
      i=numroots;
      numroots += nnumroots;
      //while (i < numroots) 
      //fprintf(stderr,"found a root at %lf\n",zc[i++]);
    }
    a = b;
    k++;
    if (k == tmpmg.nmix) 
      b = hibound;
    else
      b = mns[k];
  }

  // cleanup
  delete [] picked;
  delete [] mnind;
  delete [] mns;

  return numroots;
}
void findroot(MixGauss & mg, double lo, double hi, double acc, double * theroots, int &numroots) 
{
  int k;
  int foundone(0);
  // if there is only one mixture component then we're done
  if (mg.nmix <= 1) {
    return;
  }
  double vlo = mg.pdf(&lo);
  double vhi = mg.pdf(&hi);
  double vmid;
  double mid, aroot;
  mid = (hi-lo)/2.0+lo;
  if (hi-lo > acc) {
    vmid = mg.pdf(&mid);
    if (vlo*vhi >= 0.0) {
      // same sign
      //fprintf(stderr,"same sign...");
      // check midpoint split
      if (vmid*vlo >= 0) {
	// its the same sign
	//recurse on both halves
	findroot(mg,lo,mid,acc,theroots+numroots,numroots);
	findroot(mg,mid,hi,acc,theroots+numroots,numroots);
      } else {
	// different sign, so each half contains exactly one root
	// this is the base case
	// this is not true - there can be infinitely many roots in here :(
	foundone = ftbis(mg,lo,mid,acc,aroot);
	theroots[numroots++] = aroot;
	foundone = foundone && ftbis(mg,mid,hi,acc,aroot);
	theroots[numroots++] = aroot;
	if (!foundone) {
	  fprintf(stderr,"!aaye aye aaye ! ya des problemes ici!\n");
	  exit(0);
	} 
      }
    } else {
      // different sign, we have to look in both halves
      findroot(mg,lo,mid,acc,theroots+numroots,numroots);
      findroot(mg,mid,hi,acc,theroots+numroots,numroots);
    }
  } else {
    // if the signs are still different, then midpoint is the zero
    if (vlo*vhi < 0) 
      theroots[numroots++] = mid;
  }
}
#define MAXBIT 100
// using bisection, find the root of a mixture mg known to lie between a and b with accurary acc
// performs a maximum of MAXBIT iterations
// copied (and improved in doing so) from numerical recipies in C sect 9.1
int ftbis(MixGauss  & mg, double a, double b, double acc, double & result)
{
  int j;
  double dx,f,fmid,xmid,rtb;
  f = mg.pdf(&a);
  fmid = mg.pdf(&b);
  
  if (f*fmid >= 0.0)
    return 0;
  rtb = f < 0.0 ? (dx=b-a,a) : (dx=a-b,b);
  for (j=0; j<MAXBIT; j++) {
    xmid=rtb+(dx *= 0.5);
    fmid = mg.pdf(&xmid);
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < acc || fmid == 0.0) {
      result = rtb;
      return 1;
    }
  } 
  return 0;
}

// only works for 1D functions right now
// if typ = 0 then lo and hi are both non-infinite
// if typ = 1 then lo is taken as -inf (regardless of the input argument 'lo')
// if type = 2 then hi is taken as +inf (regardless of the input argument 'hi')
// if typ = 3 then lo is -inf and hi is +inf
// lo is assumed less than hi
double MixGauss::integrate(int typ, double lo, double hi) const
{
  int k;
  double tmp, sum, rt2,rt2sig;
  if (fvdim > 1) 
    return 0.0;
  // infinite integral - return 0 instad
  if (typ > 0 && fabs(val) > 0) 
    return 0.0;
  sum = val*(hi-lo);
  rt2 = sqrt(2);
  for (k=0; k<nmix; k++) {
    rt2sig = rt2*sqrt(mixcomps[k]->U[0][0]);
    if (typ != 1 && typ != 3) {
      sum += 0.5*mixweights[k]*erfc((lo-mixcomps[k]->u[0])/rt2sig);
    } else {
      // erfc(-inf) = 2
      sum += mixweights[k];
    }
    if (typ != 2 && typ != 3) {
      sum -= 0.5*mixweights[k]*erfc((hi-mixcomps[k]->u[0])/rt2sig);
    }
  }
  return sum;
}
MixGauss MixGauss::operator*(double g) const
{
  MixGauss newMG;
  newMG.val = g;
  return multiply((*this),newMG);
}
MixGauss  MixGauss::operator*(const MixGauss & mg)  const
{
  return multiply(*this,mg);
}
MixGauss  MixGauss::multiply(const MixGauss & mg1, const MixGauss & mg2) const
{
  int k;

  MixGauss newMG, newMG1;
  newMG.fvdim = mg1.fvdim;
  if (mg1.nmix == 0 && mg2.nmix ==0) {
    newMG.nmix = 0;
    newMG.val = mg1.val*mg2.val;
    return newMG;
  }
  if (mg1.nmix == 0) {
    newMG.copyFrom(mg2);
    //    newMG.print(stderr);
    //mg.print(stderr);
    for (k=0; k<mg2.nmix; k++) 
      newMG.mixweights[k] *= mg1.val;
    newMG.val = mg1.val*mg2.val;
    //    newMG.print(stderr);
    return newMG;
  }
  if (mg2.nmix == 0) {
    newMG.copyFrom(mg1);
    for (k=0; k<mg1.nmix; k++) 
      newMG.mixweights[k] *= mg2.val;
    newMG.val = mg1.val*mg2.val;
    return newMG;
  }
  if (mg1.fvdim == 0 && mg2.fvdim == 0 && mg2.nmix == mg1.nmix) {
    // we have a mixture (regular probability vector)
    // just mutliply each component
    newMG.copyFrom(mg1);
    for (k=0; k<mg1.nmix; k++) 
      newMG.mixweights[k] *= mg2.mixweights[k];
    newMG.val = mg1.val*mg2.val;
    return newMG;
  }
  fprintf(stderr,"SHOULD NOT GET HERE EVER - THIS DOESN'T WORK YET!");
}
MixGauss  MixGauss::operator+(const MixGauss & mg) const
{
  int i,k;
  return addorsubtract(*this,mg,true);
}
MixGauss MixGauss::operator-(const MixGauss& mg) const
{
  return addorsubtract(*this,mg,false);
}

MixGauss  MixGauss::addorsubtract(const MixGauss & mg1, const MixGauss & mg2, bool addition) const
{
  int i,k;
  MixGauss newMG;
  newMG.fvdim = mg1.fvdim;
  // check for any duplicate Gaussians
  int *isdup, ndup(0);
  if (addition) {
    newMG.val = mg1.val + mg2.val;
  } else {
    newMG.val = mg1.val - mg2.val;
  }
  if (mg2.nmix == 0) {
    newMG.nmix = 0;
    return newMG;
  }
  if (mg1.fvdim != mg2.fvdim) {
    fprintf(stderr,"can't add two MixGauss' with different fvdims\n");
    exit(0);
  }
    
  if (mg1.fvdim == 0 && mg2.fvdim == 0) {
    // if fvdim = 0 this is just a discrete belief function
    // so we just add up the mixture components
    newMG.nmix = mg1.nmix;
    newMG.getMem();
    if (mg1.nmix != mg2.nmix) {
      fprintf(stderr,"can't add two MixGauss' with fvdim = 0 with different numbers of nmix\n");
      exit(0);
    }
    for (k=0; k<mg1.nmix; k++) {
      if (addition) {
	newMG.mixweights[k] = mg1.mixweights[k] + mg2.mixweights[k];
      } else {
	newMG.mixweights[k] = mg1.mixweights[k] - mg2.mixweights[k];
      }
    }	
    return newMG;
  }
  // now, we know we have mixture components too
  // check how many duplicates there are
  isdup = new int[mg2.nmix];
  for (i=0; i<mg2.nmix; i++) 
    isdup[i] = -1;
  for (k=0; k<mg1.nmix; k++) {
    for (i=0; i<mg2.nmix; i++) {
      if (mg1.mixcomps[k] == mg2.mixcomps[i]) {
	isdup[i] = k;
	ndup++;
      }
    }
  }
  newMG.nmix = mg1.nmix+mg2.nmix-ndup;
  if (newMG.nmix < mg1.nmix) {
    fprintf(stderr,"************************ INTEGRITY COMPROMISED *********************\n");
    exit(0);
  }
  newMG.getMemNoGauss();

  for (k=0; k<mg1.nmix; k++) {
    newMG.mixweights[k] = mg1.mixweights[k];
    if (fvdim > 0)
      newMG.mixcomps[k] = mg1.mixcomps[k];
  }
  // now, copy all the new mixtures over, or add their weights
  // if its a duplicate
  i = mg1.nmix;
  for (k=0; k<mg2.nmix; k++) {
    if (isdup[k] >= 0) {
      // its a duplicate, so just add the value
      if (addition) {
	newMG.mixweights[isdup[k]] += mg2.mixweights[k];
      } else {
	newMG.mixweights[isdup[k]] -= mg2.mixweights[k];
      }
    } else {
      newMG.mixweights[i] = mg2.mixweights[k];
      newMG.mixcomps[i] = mg2.mixcomps[k];
      i++;
    }
  }
  if (mg2.nmix > 0) 
    delete [] isdup;
  newMG.removeZeros();
  return newMG;
}

// removes all mixture components with zero weghts
void MixGauss::removeZeros()
{
  int k,i;
  int numzeros = 0;
  // not necessary if fvdim is 0
  if (fvdim == 0) 
    return;
  int *iszero = new int[nmix];
  for (k=0; k<nmix; k++) {
    iszero[k] = ((mixweights[k] == 0.0) ? 1 : 0);
    numzeros = numzeros+iszero[k];
  }
  if (numzeros) {
    // make a copy
    MixGauss tmpMG(*this);
    if (gotMem || gotMemNoGauss)
      delMem();
    nmix = tmpMG.nmix-numzeros;
    getMem();
    i=0;
    for (int k=0; k<tmpMG.nmix; k++) {
      if (!iszero[k]) {
	mixweights[i] = tmpMG.mixweights[k];
	mixcomps[i++] = tmpMG.mixcomps[k];
      }
    }
  }
  delete [] iszero;
  
}
void MixGauss::operator=(const MixGauss & mg) 
{
  copyFrom(mg);
}
bool MixGauss::operator==(const MixGauss & mg) 
{
  int i;
  bool res = (mg.nmix == nmix) && (mg.fvdim == fvdim);
  res = res && (fabs(mg.val-val)<1.0e-12);
  if (!res || nmix ==  0) {
    return res;
  }
  if (fvdim == 0) {
    for (i=0; res && i<nmix; i++) 
      res = res && (fabs(mixweights[i]-(mg.mixweights[i])) < 1e-12);
    return res;
  }
  // new method - subtract the two and see if there are no mixture 
  // components left
  if (res) {
    MixGauss tmpmg = addorsubtract(*this,mg,false);
    res = res && (tmpmg.nmix == 0);
  }
  return res;
    
}
MixGauss MixGauss::operator-() const
{
  MixGauss newMG(-1.0);
  newMG = newMG*(*this);
  return newMG;
}
MixGauss MixGauss::operator/(MixGauss& mg) const
{
  int k;
  MixGauss newMG;
  newMG.copyFrom(*this);
  for (k=0; k<nmix; k++) 
    newMG.mixweights[k] /= mg.val;
  newMG.val = val/mg.val;
  return newMG;
}
bool MixGauss::operator>=(const MixGauss& mg) {
  bool v = (val >= mg.val);
  return v;
}
bool MixGauss::operator<=(const MixGauss& mg) {
  bool v = (val <= mg.val);
  return v;
}
bool MixGauss::operator>(const MixGauss& mg) {
  bool v = ((*this)<=mg);
  return !v;
}
bool MixGauss::operator<(const MixGauss& mg) {
  bool v = ((*this)>=mg);
  return !v;
}
bool MixGauss::operator!=(const MixGauss& mg) {
  bool v = ((*this)==mg);
  return !v;
}
double MixGauss::get_val() {
  return val;
}
double MixGauss::get_min() {
  return val;
}
double MixGauss::get_max() {
  return val;
}
MixGauss MixGauss::absVal() const {
  MixGauss newMG(*this);
  int k;
  newMG.val = fabs(newMG.val);
  for (k=0; k<nmix; k++) {
    newMG.mixweights[k] = fabs(newMG.mixweights[k]);
  }
  return newMG;
}

MixGauss MixGauss::operator/(double num) const
 {
  MixGauss newMG;
  newMG.val = val/num;
  return newMG;
}
void MixGauss::print(FILE *fd) const
{
  fprintf(fd,"%d %d\n",nmix,fvdim);
  fprintf(fd,"%f\n",val);
  for (int i=0; i<nmix; i++) {
    fprintf(fd,"%f\n",mixweights[i]);
    if (fvdim > 0)
      mixcomps[i]->print(fd);
  }
}
void MixGauss::printValMixWeights(FILE *fd) const
{
  fprintf(fd,"%f ",val);
  for (int i=0; i<nmix; i++) {
    fprintf(fd,"%f ",mixweights[i]);
  }
  fprintf(fd,"\n");
}
char *MixGauss::toString() const
{
  int numlen = 4;
  // only print first ten nmixs
  int pnmix = (10 < nmix) ? 10 : nmix;
  int strlen = numlen+pnmix*(numlen+10+fvdim*2*numlen)+10;
  //fprintf(stderr,"numlen: %d nmix %d pnmix %d fvdim %d allocating string length %d\n",numlen,nmix,pnmix,fvdim,strlen);
  //fprintf(stderr,"allocating string length %d\n",strlen);
  char *newstring = new char[strlen+128];
  int pos = 0;
  pos += sprintf(newstring,"%10.9f ",val);
  for (int i=0; i<pnmix; i++) {
    if (i >0)
      pos += sprintf(newstring+pos,"+");
    if (fvdim == 0) {
      pos += sprintf(newstring+pos,"%3.2f ",mixweights[i]);
    } else {
      pos += sprintf(newstring+pos,"%3.2f*N([",mixweights[i]);
      for (int j=0; j<fvdim; j++) {
	if (j>0)
	  pos += sprintf(newstring+pos,",");
	pos += sprintf(newstring+pos,"%3.2f",mixcomps[i]->u[j]);
      }
      pos += sprintf(newstring+pos,"],[");
      for (int j=0; j<fvdim; j++) {
	if (j>0)
	  pos += sprintf(newstring+pos,",");
	pos += sprintf(newstring+pos,"%3.2f",mixcomps[i]->U[j][j]);
      }
      pos += sprintf(newstring+pos,"])");
    }
  }
  if (pnmix < nmix)
    sprintf(newstring+pos," ... ");
  return newstring;
}
int MixGauss::readString(FILE *fp, char * buf) {
  int retval;
  retval = fscanf(fp,"%s\n",buf);
  return retval;
}
void MixGauss::parseString(char *buf) {
} 



double MixGauss::hashCode() const
{
  double tmp = 0.0;
  double ctmp;
  double ttmp = 1.0;
  tmp += val;
  for (int i=0; i<nmix; i++) {
    tmp += ttmp*mixweights[i];
    if (fvdim > 0) {
      ctmp = (double) ((long) mixcomps[i]);
      tmp += ttmp*ctmp;
    }
    //ttmp *= HPI;

  }
  return tmp;
}

MixGauss ceill(const MixGauss& source){
  
  MixGauss newMixGauss;

  newMixGauss.val = ceil(source.val);
  return newMixGauss;
}

// complementary error fuction. Taken from
// http://www.csit.fsu.edu/~burkardt/cpp_src/dcdflib/dcdflib.html
// (dcdflib)
double erfc (double x)

//****************************************************************************
//
//  Purpose:
// 
//    ERFC1 evaluates the complementary error function.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int *IND, chooses the scaling.
//    If IND is nonzero, then the value returned has been multiplied by
//    EXP(X*X).
//
//    Input, double *X, the argument of the function.
//
//    Output, double ERFC1, the value of the complementary 
//    error function.
//
{
  static double c = .564189583547756e0;
  static double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  static double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  static double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  static double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  static double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  static double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  static int K1 = 1;
  static double erfc1,ax,bot,e,t,top,w;

//
//                     ABS(X) .LE. 0.5
//
    ax = fabs(x);
    if(ax > 0.5e0) goto S10;
    t = x*x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erfc1 = 0.5e0+(0.5e0-x*(top/bot));
    //if(*ind != 0) erfc1 = exp(t)*erfc1;
    return erfc1;
S10:
//
//                  0.5 .LT. ABS(X) .LE. 4
//
    if(ax > 4.0e0) goto S20;
    top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
      7];
    bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
      7];
    erfc1 = top/bot;
    goto S40;
S20:
//
//                      ABS(X) .GT. 4
//
    if(x <= -5.6e0) goto S60;
    //if(*ind != 0) goto S30;
    if(x > 100.0e0) goto S70;
    if(x*x > -exparg(&K1)) goto S70;
S30:
    t = pow(1.0e0/ x,2.0);
    top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
    bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
    erfc1 = (c-t*top/bot)/ax;
S40:
//
//                      FINAL ASSEMBLY
//
    goto S50; //if(*ind == 0) goto S50;
    if(x < 0.0e0) erfc1 = 2.0e0*exp(x*x)-erfc1;
    return erfc1;
S50:
    w = x*x;
    t = w;
    e = w-t;
    erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
    if(x < 0.0e0) erfc1 = 2.0e0-erfc1;
    return erfc1;
S60:
//
//             LIMIT VALUE FOR LARGE NEGATIVE X
//
    erfc1 = 2.0e0;
    //if(*ind != 0) erfc1 = 2.0e0*exp(x*x);
    return erfc1;
S70:
//
//             LIMIT VALUE FOR LARGE POSITIVE X
//                       WHEN IND = 0
//
    erfc1 = 0.0e0;
    return erfc1;
}
double exparg ( int *l )

//****************************************************************************
//
//  Purpose:
// 
//    EXPARG returns the largest or smallest legal argument for EXP.
//
//  Discussion:
//
//    Only an approximate limit for the argument of EXP is desired.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int *L, indicates which limit is desired.
//    If L = 0, then the largest positive argument for EXP is desired.
//    Otherwise, the largest negative argument for EXP for which the
//    result is nonzero is desired.
//
//    Output, double EXPARG, the desired value.
//
{
  static int K1 = 4;
  static int K2 = 9;
  static int K3 = 10;
  static double exparg,lnb;
  static int b,m;

    b = ipmpar(&K1);
    if(b != 2) goto S10;
    lnb = .69314718055995e0;
    goto S40;
S10:
    if(b != 8) goto S20;
    lnb = 2.0794415416798e0;
    goto S40;
S20:
    if(b != 16) goto S30;
    lnb = 2.7725887222398e0;
    goto S40;
S30:
    lnb = log((double)b);
S40:
    if(*l == 0) goto S50;
    m = ipmpar(&K2)-1;
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
S50:
    m = ipmpar(&K3);
    exparg = 0.99999e0*((double)m*lnb);
    return exparg;
}
int ipmpar ( int *i )

//****************************************************************************
//
//  Purpose:
//  
//    IPMPAR returns integer machine constants. 
//
//  Discussion:
//
//    Input arguments 1 through 3 are queries about integer arithmetic.
//    We assume integers are represented in the N-digit, base-A form
//
//      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
//
//    where 0 <= X(0:N-1) < A.
//
//    Then:
//
//      IPMPAR(1) = A, the base of integer arithmetic;
//      IPMPAR(2) = N, the number of base A digits;
//      IPMPAR(3) = A**N - 1, the largest magnitude.
//
//    It is assumed that the single and double precision floating
//    point arithmetics have the same base, say B, and that the
//    nonzero numbers are represented in the form
//
//      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
//
//    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
//    EMIN <= E <= EMAX.
//
//    Input argument 4 is a query about the base of real arithmetic:
//
//      IPMPAR(4) = B, the base of single and double precision arithmetic.
//
//    Input arguments 5 through 7 are queries about single precision
//    floating point arithmetic:
//
//     IPMPAR(5) = M, the number of base B digits for single precision.
//     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
//     IPMPAR(7) = EMAX, the largest exponent E for single precision.
//
//    Input arguments 8 through 10 are queries about double precision
//    floating point arithmetic:
//
//     IPMPAR(8) = M, the number of base B digits for double precision.
//     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
//     IPMPAR(10) = EMAX, the largest exponent E for double precision.
//
//  Reference:
//
//    Fox, Hall, and Schryer,
//    Algorithm 528,
//    Framework for a Portable FORTRAN Subroutine Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, 1978, pages 176-188.
//
//  Parameters:
//
//    Input, int *I, the index of the desired constant.
//
//    Output, int IPMPAR, the value of the desired constant.
//
{
  static int imach[11];
  static int ipmpar;
//
//     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
//       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
//       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). 

   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
    ipmpar = imach[*i];
    return ipmpar;
}


#ifdef USE_MATHEMATICA

// gets roots using mathematica - 
// you have to call initialize() before calling this to 
// initialize the link
int getroots_mathematica(char *inputstring, double **theroots)
{	
  int pkt, n, prime, expt;
  long len, lenin, lenp, k, j;
  double lobound, hibound;
  
  MLPutFunction( lp, "EvaluatePacket", 1);
  MLPutFunction(lp,"ToExpression",1);
  MLPutString(lp,inputstring);
  MLEndPacket( lp);
  
  while( (pkt = MLNextPacket( lp)) && pkt != RETURNPKT && pkt != RETURNTEXTPKT) {
    MLNewPacket( lp);
  }
  if (!pkt) {
    fprintf(stderr,"error in there!\n");
  }
  if ( ! MLCheckFunction( lp, "Interval", &lenin)) error(lp);
  //fprintf(stderr,"got Interval back with %d intervals\n",lenin);
  *theroots = new double[lenin];
  for (j = 1; j <= lenin; j++) {
    if ( ! MLCheckFunction( lp, "List", &len)) error(lp);
    //fprintf(stderr,"interval %d has %d elements\n",j,len);
    if (MLGetDouble( lp, &lobound)  &&  MLGetDouble( lp, &hibound)){
      //fprintf(stderr,"%f %f %f\n", lobound, hibound,(hibound+lobound)/2.0);
      (*theroots)[j-1] = (hibound+lobound)/2.0;
    }else{
      error(lp);
    }
  }
  return lenin;
}
// MATHEMATICA STUFF

static void error( MLINK lp)
{
  if( MLError( lp)){
    fprintf( stderr, "Error detected by MathLink: %s.\n",
	     MLErrorMessage(lp));
  }else{
    fprintf( stderr, "Error detected by this program.\n");
  }
  exit(3);
}


static void deinit( void)
{
  if( ep) MLDeinitialize( ep);
}


static void closelink( void)
{
  if( lp) MLClose( lp);
}
void terminate_mlink() {
  MLPutFunction( lp, "Exit", 0);
  //closelink();
}

static void init_and_openlink( int argc, char* argv[])
{
  long err;

  ep =  MLInitialize( (MLParametersPointer)0);
  if( ep == (MLENV)0) exit(1);
  atexit( deinit);

  lp = MLOpenArgv( ep, argv, argv + argc, &err);
  if(lp == (MLINK)0) exit(2);
  atexit( closelink);
	
}
void initialize_mlink() {
  int pkt, n, prime, expt;
  int argc = 5;
  char *argv[argc];
  argv[0] = "gspudd";
  argv[1] = "-linkmode";
  argv[2] = "launch";
  argv[3] = "-linkname";
  argv[4] = "/h/23/jhoey/Mathematica/bin/math -mathlink";

  init_and_openlink( argc, argv);
  
  MLPutFunction( lp, "EvaluatePacket", 1);
  MLPutFunction(lp,"Get",1);
  MLPutString(lp,"NumericalMath`IntervalRoots`");
  //MLPutFunction(lp,"ToExpression",1);
  //MLPutString(lp,"Get[NumericalMath`IntervalRoots`]");
  MLEndPacket( lp);
  while ( (pkt= MLNextPacket(lp))  && pkt != RETURNPKT) {
    MLNewPacket(lp);
  }
  // throw away the Null return packet
  MLNewPacket(lp);
  
  fprintf(stderr,"loaded package numerical math\n");
  
}

#endif
