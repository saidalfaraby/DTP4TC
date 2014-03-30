#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include "mixgauss.h"

int main(int argc, char *argv[]) 
{
  int i, j, k, nmix(2);
  double weight, *mn, **cv, rnd, res;
  int fvdim = 2;
  mn = new double[fvdim];
  cv = new double*[fvdim];
  for (i=0; i<fvdim; i++) {
    cv[i] = new double[fvdim]; 
    for (k=0; k<fvdim; k++) 
      cv[i][k] = 0.0;
  }
  struct timeval sd;
  gettimeofday(&sd,NULL);
  int rseed = (int) sd.tv_sec;
  srand(rseed);

  Gaussian gg(2);
  mn[0] = 0.0;
  mn[1] = 0.0;
  cv[0][0] = 6.563;
  cv[0][1] = -0.211;
  cv[1][0] = -0.211;
  cv[1][1] = 0.723;
  gg.setMean(mn);
  gg.setCovariance(cv);
  double *sample;
  sample = gg.drawSample();
  //fprintf(stderr,"%f %f\n",sample[0],sample[1]);

  /*
  MixGauss mgg(nmix,fvdim);
  for (i=0; i<nmix; i++) {
    weight = ((double) rand())/((double) RAND_MAX+1.0);
    for (j=0; j<fvdim; j++) {
      mn[j] =  20*(((double) rand())/((double) RAND_MAX+1.0) - 0.5);
      cv[j][j] = 10*((double) rand())/((double) RAND_MAX+1.0);
    }
    mgg.setGaussian(i,weight,mn,cv);
  }
  for (j=0; j<fvdim; j++) {
    mn[j] =  20*(((double) rand())/((double) RAND_MAX+1.0) - 0.5);
  }
  res = mgg.pdf(mn);
  fprintf(stderr,"hash to %lf     -------- result is %lf\n",mgg.hashCode(),res);
  // compare two gaussians
  MixGauss nmgg(nmix,fvdim);
  for (i=0; i<nmix; i++) {
    weight = ((double) rand())/((double) RAND_MAX+1.0);
    for (j=0; j<fvdim; j++) {
      mn[j] =  20*(((double) rand())/((double) RAND_MAX+1.0) - 0.5);
      cv[j][j] = 10*((double) rand())/((double) RAND_MAX+1.0);
    }
    nmgg.setGaussian(i,weight,mn,cv);
  }
  mgg.print(stderr);
  nmgg.print(stderr);
  if (mgg == nmgg) 
    fprintf(stderr,"two different gaussians equal?\n");
  else 
    fprintf(stderr,"two different gaussians are not equal\n");
  MixGauss mgg;
  mgg = nmgg;
  mgg.print(stderr);
  nmgg.print(stderr);
  if (mgg == nmgg) 
    fprintf(stderr,"two same gaussians are equal\n");
  else 
    fprintf(stderr,"two same gaussians not equal?\n");
  MixGauss nnmgg(0,fvdim);
  nnmgg = mgg+nmgg;
  nnmgg.print(stderr);
  MixGauss n2mgg(0,fvdim);
  mgg.setValue(1.2);
  n2mgg.setValue(10.0);

  
  if (nnmgg.absVal() < mgg.absVal()) {
    fprintf(stderr,"its less!\n");
  } else {
    fprintf(stderr,"its more!\n");
  }
  */
}




