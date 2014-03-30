#include "gaussian.h"
Gaussian::Gaussian()
{
  fvdim = 0;
  Usr = NULL;
}
Gaussian::Gaussian(int fvd)
{
  fvdim = fvd;
  getMem();
}
Gaussian::Gaussian(Gaussian *g)
{
  fvdim = g->fvdim;
  getMem();
  copyFrom(g);
}
Gaussian::~Gaussian() {
  int k;
  for (k=0; k<fvdim; k++) {
    delete [] UI[k];
    delete [] U[k];
  }

  delete [] u;
  delete [] UI;
  delete [] U;
}
void Gaussian::copyFrom(Gaussian *g)
{
  int i,j;
  for (j=0; j<fvdim; j++) {
    u[j] = g->u[j];
    for (i=0; i<fvdim; i++) {
      U[j][i] = g->U[j][i];
    }
  }
  Usr = NULL;
  invert_covar();
}
void Gaussian::zero()
{
  int i,j;
  det = 0.0;
  for (j=0; j<fvdim; j++) {
    u[j] = 0.0;
    for (i=0; i<fvdim; i++) {
      U[j][i] = 0.0;
      UI[j][i] = 0.0;
    }
  }
}
void Gaussian::getMem() {
  u = new double[fvdim];
  U = new double*[fvdim];
  UI = new double*[fvdim];
  for (int k=0; k<fvdim; k++) {
    U[k] = new double[fvdim];
    UI[k] = new double[fvdim];
  }
  Usr = NULL;
}
const double Gaussian::hashCode()
{
  double tmp = 0.0;
  double ttmp = 1.0;
  for (int i=0; i<fvdim; i++) {
    tmp += ttmp*u[i];
    for (int j=0; j<fvdim; j++) {
      tmp += ttmp*U[i][j];
    }
    ttmp *= HPI;
  }
  //fprintf(stderr,"hash code is %g\n",tmp);
  return tmp;
}
Gaussian Gaussian::multiply(const Gaussian & g, double & fact)
{
  Gaussian newGauss;
  double newsig, tmp;
  double *newmu = new double[fvdim];
  double *tmp1 = new double[fvdim];
  double *tmp2 = new double[fvdim];
  double *tmp3 = new double[fvdim];
  newGauss.fvdim = fvdim;
  newGauss.getMem();
  addMatrix(fvdim,UI, g.UI, newGauss.UI);
  newGauss.det = invertMatrix(fvdim, newGauss.UI, newGauss.U);
  multiplyMatrixVector(fvdim, UI, u, tmp1);
  multiplyMatrixVector(fvdim, g.UI, g.u, tmp2);
  addVector(fvdim,tmp1,tmp2,tmp3);
  multiplyMatrixVector(fvdim, newGauss.U, tmp3, newGauss.u);
  tmp = multiplyMatrix2Vector(fvdim, UI, u);
  tmp += multiplyMatrix2Vector(fvdim, g.UI, g.u);
  tmp -= multiplyMatrix2Vector(fvdim, newGauss.UI, newGauss.u);
  fact = (sqrt(newGauss.det/(2*HPI*det*g.det)))*exp(-0.5*tmp);
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3;
  return newGauss;
}
bool Gaussian::operator==(const Gaussian & g)
{
  int i,j;
  bool res = true;
  for (i=0; res && i<fvdim; i++) {
    res = res && (fabs(u[i]-(g.u[i])) < 1e-12);
    for (j=0; res && j<fvdim; j++) 
      res = res && (fabs(U[i][j]-(g.U[i][j])) < 1e-12);
  }
  return res;
}
void Gaussian::setMean(double *mn) {
  int i;
  for (i=0; i<fvdim; i++)
    u[i] = mn[i];
}
void Gaussian::setCovariance(double **cv) {
  int i,j;
  for (i=0; i<fvdim; i++)
    for (j=0; j<fvdim; j++)
      U[i][j] = cv[i][j];
  invert_covar();
}
// here cv is a fvdim*(fvdim+1)/2 array specifying the 
// upper (or lower) triangular elements of the covariance matrix 
// c00 c01 c02 ... c0fvdim c11 c12 c13...c1fvdim c22 c23 ... etc
void Gaussian::setCovariance(double *cv) {
  int i,j,k(0);
  // if there is only one dimension, we expect this
  // to be the variance, not the covariance, so we square
  if (fvdim == 1) {
    cv[0] = cv[0]*cv[0];
  }
  for (i=0; i<fvdim; i++) {
    for (j=i; j<fvdim; j++) {
      U[i][j] = *(cv+k);
      k++;
    }
  }
  // copy over other diagonals
  for (i=0; i<fvdim; i++) 
    for (j=0; j<i; j++) 
      U[i][j] = U[j][i];
  invert_covar();
}
void Gaussian::invert_covar()
{
  det = invertMatrix(fvdim, U, UI);
  // compute multiplicative factor for Gaussians
  oPFg = 1.0/(sqrt(absv(det))*pow((2.0*HPI),((double) fvdim)/2.0));
}
// evaluates the P(obs | gaussian in thetaZ)
double Gaussian::pdf(double *obs)
{
  int k;
  double sum=0.0;
  for (int i=0; i<fvdim; i++) 
    for (k=0; k<fvdim; k++) 
      sum += (obs[k]-u[k])*UI[k][i]*(obs[i]-u[i]);
  
  sum = -0.5*sum;
  
  return oPFg*exp(sum);
}

double *Gaussian::drawSample()
{
  int k,l;
  double *sample = new double[fvdim];
  double z;
  if (Usr == NULL)
    computeSqrtCovar();
  // draw zero-mean variance 1 samples using boxmuller for each dimension
  for (k=0; k<fvdim; k++) 
    sample[k] = 0.0;
  for (k=0; k<fvdim; k++) {
    z = drawGaussianSample();
    for (l=0; l<fvdim; l++) 
      sample[l] += z*Usr[k][l];
  }
  // add the mean
  for (k=0; k<fvdim; k++) 
    sample[k] += u[k];
  return sample;
}
 // compute square root of covariance matrix --> Usr
void Gaussian::computeSqrtCovar() 
{
  int k,i,j;
  Usr = new double*[fvdim];
  for (k=0; k<fvdim; k++) 
    Usr[k] = new double[fvdim];
  if (fvdim < 2) {
    Usr[0][0] = sqrt(U[0][0]);
    return;
  }
  int m=fvdim;
  double *tempe = dvectorNR(1,m);
  double **temp = dmatrix(1,m,1,m);
  double **temp2 = dmatrix(1,m,1,m);
  double **eigvals = new double*[m];
  double **eigvecs = new double*[m];
  double **eigvecsp = new double*[m];
  for (i=0; i<m; i++) {
    eigvals[i] = new double[m];
    eigvecs[i] = new double[m];
    eigvecsp[i] = new double[m];
  }

  int nrot;

  //fprintf(stderr,"computing sqrt of \n");
  //fprintf(stderr,"\n%f %f\n%f %f\n",U[0][0],U[0][1],U[1][0],U[1][1]);

  for (i=1; i<=m; i++) 
    for (j=1; j<=m; j++) 
      temp[i][j] = U[i-1][j-1];
  jacobi(temp, m, tempe, temp2, &nrot);
  // temp2 are the eigenvectors, eigvs are the eigenvalues
  for (i=0; i<m; i++) {
    for (j=0; j<m; j++) {
      if (i==j) {
	eigvals[i][j] = sqrt(tempe[i+1]);
      } else {
	eigvals[i][j] = 0.0;
      }
      eigvecs[i][j] = temp2[i+1][j+1];
    }
  }
  // now perform the multiplications
  transposeMatrix(m,eigvecs,eigvecsp);
  multiplyMatrix(m,eigvals,eigvecsp,Usr);
  multiplyMatrix(m,eigvecs,Usr,eigvecsp);
  for (i=0; i<m; i++) 
    for (j=0; j<m; j++) 
      Usr[i][j] = eigvecsp[i][j];
  
  free_dmatrix(temp,1,m,1,m);
  free_dmatrix(temp2,1,m,1,m);
  free_dvectorNR(tempe,1,m);

  for (i=0; i<m; i++) {
    delete [] eigvals[i];
    delete [] eigvecs[i];
    delete [] eigvecsp[i];
  }
  delete [] eigvals;
  delete [] eigvecs;
  delete [] eigvecsp;

  /*
  double l1,l2,den1,den2;
  double e11,e12,e21,e22;
  double a = U[0][0];
  double b = U[1][1];
  double c = U[0][1];
  l1 = 0.5*((a+b)+sqrt((a+b)*(a+b)-4*(a*b-c*c)));
  l2 = 0.5*((a+b)-sqrt((a+b)*(a+b)-4*(a*b-c*c)));
  fprintf(stderr,"\n%f %f\n",l1,l2);
  den1 = (a-l1)*(a-l1)+c*c;
  den2 = (a-l2)*(a-l2)+c*c;
  e11 = sqrt(c*c/den1);
  e12 = sqrt((a-l1)*(a-l1)/den1);
  e21 = sqrt(c*c/den2);
  // one has to be negative
  e22 = -sqrt((a-l2)*(a-l2)/den2);
  fprintf(stderr,"\n%f %f\n%f %f\n",e11,e21,e12,122);
  l1 = sqrt(l1);
  l2 = sqrt(l2);

  Usr[0][0] = l1*e11*e11+l2*e21*e21;
  Usr[0][1] = l1*e11*e12+l2*e22*e21;
  Usr[1][0] = l1*e11*e12+l2*e21*e22;
  Usr[1][1] = l1*e12*e12+l2*e22*e22;
  fprintf(stderr,"\n%f %f\n%f %f\n",Usr[0][0],Usr[0][1],Usr[1][0],Usr[1][1]);
  */
}

// only draws from diagonal covariance Gaussians
/*
double *Gaussian::drawSample() 
{
  double *sample = new double[fvdim];
  for (int i=0; i<fvdim; i++) {
    sample[i] = u[i]+drawGaussianSample(U[i][i]);
  }
  return sample;
}
*/
double Gaussian::logpdf(double *obs)
{
  int i,k;
  double sum=mahalDist(obs);
  return log2(oPFg)-0.5*sum;
}
double Gaussian::mahalDist(double *obs)
{
  int i,k;
  double sum=0.0;
  for (i=0; i<fvdim; i++) 
    for (k=0; k<fvdim; k++) 
      sum += (obs[k]-u[k])*UI[k][i]*(obs[i]-u[i]);
  return sum;
}
void Gaussian::print(FILE *fd)
{
  int i,k;
  double sum=0.0;
  for (i=0; i<fvdim; i++)
    fprintf(fd,"%f ",u[i]);
  fprintf(fd,"\n");
  for (i=0; i<fvdim; i++) {
    for (k=0; k<fvdim; k++) 
      fprintf(fd,"%f ",U[i][k]);
    fprintf(fd,"\n");
  }
}
void multiplyMatrix(int m, double **mat1, double **mat2, double **mat)
{
  int i,j,k;
  for (i=0; i<m; i++) {
    for (j=0; j<m; j++) {
      mat[i][j] = 0.0;
      for (k=0; k<m; k++) 
	mat[i][j] += mat1[i][k]*mat2[k][j];
    }
  }
}
void transposeMatrix(int m, double **mat1, double **mat)
{
  int i,j;
  for (i=0; i<m; i++) 
    for (j=0; j<m; j++) 
      mat[i][j] = mat1[j][i];
}
void multiplyMatrixVector(int m, double **mat1, double *vec1, double *vec)
{
  int i,j,k;
  for (i=0; i<m; i++) {
    vec[i] = 0.0;
    for (k=0; k<m; k++) 
      vec[i] += mat1[i][k]*vec1[k];
  }
}
double multiplyMatrix2Vector(int m, double **mat1, double *vec1)
{
  int i,j,k;
  double val;
  double *tmp = new double[m];
  multiplyMatrixVector(m,mat1,vec1,tmp);
  val = multiplyVector(m,vec1,tmp);
  delete [] tmp;
  return val;
}
double multiplyVector(int m, double *vec1, double *vec2)
{
  int i;
  double val = 0;
  for (i=0; i<m; i++)
    val += vec1[i]*vec2[i];
  return val;
}

void addVector(int m, double *vec1, double *vec2, double *vec)
{
  int i,j;
  for (i=0; i<m; i++) 
    vec[i] = vec1[i] + vec2[i];
}
void addMatrix(int m, double **mat1, double **mat2, double **mat)
{
  int i,j;
  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      mat[i][j] = mat2[i][j] + mat1[i][j];
}

// inverts the m x m matrix in mat, puts the result in mati, and
// returns the determinant
double invertMatrix(int m, double **mat, double **mati)
{
  int i,j,k,ts;
  double **temp, d, *col;
  int *indx;

  if (m == 1) {
    mati[0][0] = 1.0/mat[0][0];
    return mat[0][0];
  }
  col = dvectorNR(1,m);
  indx = ivectorNR(1,m);
  temp = dmatrix(1,m,1,m);
  
  for (i=1; i<=m; i++)   
    for (j=1; j<=m; j++) 
      temp[i][j] = mat[i-1][j-1];
  
  dludcmp(temp, m, indx, &d);
  
  for (j=1; j <= m; j++) {
    d *= temp[j][j];
    for (i=1; i<= m; i++) 
      col[i] = 0.0;
    col[j] = 1.0;
    dlubksb(temp,m,indx,col);
    for (i=1; i<=m; i++) 
      mati[i-1][j-1] = (double) (col[i]);
  } 
  
  if (d < 0 )
    fprintf(stderr," what??? determinant < 0?? \n");
  
  free_dmatrix(temp,1,m,1,m);
  free_ivectorNR(indx,1,m);
  free_dvectorNR(col,1,m);
  return d;
}
// draws a gaussian sample using the box Muller method
double drawGaussianSample(double sig) {
  return sig*drawGaussianSample();
}
// draws gaussian sample with variance 1
double drawGaussianSample() {
   double rn1 = rand()/(RAND_MAX+1.0);
   double rn2 = rand()/(RAND_MAX+1.0);
   return (sqrt(-2*log(rn1))*sin(2*3.141592*rn2));
}
// draws a gaussian sample from a zero mean sig = 1.0 
/*
double drawGaussianSample(double sig) {
  double rn = rand()/(RAND_MAX+1.0);
  int numbins(100);
  
  double sum = 0.0;
  double x = -5.0*sqrt(sig);
  double dx = 10.0*sqrt(sig)/numbins;
  double *gg = new double[numbins];
  double factor = 1.0/sqrt(2.0*3.1415*sig);
  int i=0;
  for (i=0; i<numbins; i++) {
    gg[i] = factor*exp(-0.5*x*x/sig);
    sum += gg[i];
    x += dx;
  }
  for (i=0; i<numbins; i++) {
    gg[i] /= sum;
  }
  sum = 0.0;
  i=0;
  x = -5.0*sqrt(sig);
  while (i < numbins && rn > sum) {
    sum += gg[i];
    x += dx;
    i++;
  }
  delete [] gg;
  return x;
}
*/
