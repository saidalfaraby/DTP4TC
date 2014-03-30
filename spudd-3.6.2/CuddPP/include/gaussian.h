#ifndef _GAUSSIAN
#define _GAUSSIAN

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

#define absv(x) (((x)<0) ? -(x) : (x))
#define HPI 3.1415926535897
class Gaussian {
 public:
  Gaussian();
  Gaussian(int fvdim);
  Gaussian(Gaussian *g);
  ~Gaussian();
  void getMem();
  const double hashCode();
  void copyFrom(Gaussian *g);
  // invert the covariance in U --> UI
  void invert_covar();
  virtual void zero();
  
  void setMean(double *);
  void setCovariance(double **);
  void setCovariance(double *);
  
  // compute P(obs | M)
  double pdf(double *obs);
  double logpdf(double *obs);

  double* drawSample();

  double mahalDist(double *obs);
  // accessors
  double *get_mean();
  double ** get_covar();
  double ** get_covar_inv();
  double get_det();

  void computeSqrtCovar();

  bool operator==(const Gaussian&);
  Gaussian multiply(const Gaussian &, double &);
  void print(FILE *);
  // variables
  // dimension
  int fvdim;
  // mean
  double *u;
  // covariance and inverse
  double **U, **UI; 
  // square root of covariance
  double **Usr;
  // determinant of covariance
  double det;
  // multiplicative factor
  double oPFg;
  // integer index for each unique Gaussian
  //int gaussianUniqueIndex;
};
// helper functions for matrix manipulation
void multiplyMatrix(int m, double **mat1, double **mat2, double **mat);
void multiplyMatrixVector(int m, double **mat1, double *vec1, double *vec);
double multiplyMatrix2Vector(int m, double **mat1, double *vec1);
void addVector(int m, double *vec1, double *vec2, double *vec);
void addMatrix(int m, double **mat1, double **mat2, double **mat);
double multiplyVector(int m, double *vec1, double *vec2);
double invertMatrix(int m, double **mat, double **mati);
void transposeMatrix(int m, double **mat1, double **mat);
double drawGaussianSample(double sig);
double drawGaussianSample();

#endif
