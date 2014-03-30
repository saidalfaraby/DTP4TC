#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
void fillinimage(char *im, int nx, int ny, int stepsize, char *);

double * readObsImageRaw(char *filename, int nx, int ny)
{
  // read im as raw doubles
  double *vv = new double[nx*ny];
  int fd = open(filename,O_RDONLY);
  read(fd,vv,nx*ny*sizeof(double));
  return vv;
}

void readObsImageMgs(char *filename, int nx, int ny, int *pixels)
{

  // read im as raw ints
  int fd = open(filename,O_RDONLY);
  read(fd,pixels,nx*ny*sizeof(int));
}
void readEdgeImagePGM(char *filename, int nx, int ny, char *eimage)
{
  FILE *stream = fopen(filename, "r" );
  if ( !stream )
    return;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(eimage, sizeof(char), ny*nx, stream);
  fclose(stream);
}

void writeImagePPM(char *filename, char *im, int nx, int ny)
{
  FILE *stream = fopen(filename, "w" );
  if ( !stream )
    return;
  fprintf( stream, "P6\n%d %d\n255\n", nx, ny);
  fwrite(im, sizeof(char), 3*ny*nx, stream);
  fclose(stream);
}

// usage: 
// rescaleims  fileprefix (e.g. obsreg34) numalphas numbeliefs numactions
int main(int argc, char *argv[])
{
  char *fileprefix, *edgeimprefix;
  char filename[256];
  // prefixes
  fileprefix = *++argv;

  int nmgs = atoi(*++argv);

  int nb = atoi(*++argv);
  int na = atoi(*++argv);
  int ov = atoi(*++argv);
  int i,j, index,a,b, k;
  double ***dvals;
  double maxval, minval;
  char *im;
  int nx, ny;
  nx = 160;
  ny = 120;
  // gather all raw images
  dvals = new double**[nb];
  for (b=0; b<nb; b++) {
    dvals[b] = new double*[na];
    for (a=0; a<na; a++) {
      //sprintf(filename,"%s%d_%d_%d.dat",fileprefix,b,a,ov);
      sprintf(filename,"%s%d_%d.dat",fileprefix,b,ov);
      dvals[b][a] = readObsImageRaw(filename,nx,ny);
    }
  }
  // read in edge images
  char **eims = new char*[na];
  for (a=0; a<na; a++) {
    eims[a] = new char[nx*ny];
    sprintf(filename,"edges%d.pgm",a);
    readEdgeImagePGM(filename, nx, ny, eims[a]);
  }

  
  im = new char[nx*ny*3];

  // find max and min vals across all images
  maxval = dvals[0][0][0];
  minval = dvals[0][0][0];
  int stepsize = 1;
  int modu = 0;
  for (b=0; b<nb; b++) {
    for (a=0; a<na; a++) {
      modu = 0;
      for (i=0; i<ny; i+=stepsize) {
	for (j=modu; j<nx; j+=stepsize) {
	  index = i*nx+j;
	  if (dvals[b][a][index] > maxval) 
	    maxval = dvals[b][a][index];
	  if (dvals[b][a][index] < minval) 
	    minval = dvals[b][a][index];
	}
	modu = (modu+1)%stepsize;
      }
    }
  }
  fprintf(stderr,"maxval %g minval %g\n",maxval,minval);

  int colmap[3];
  int numcolorshades = 4;
  int colorshades[4] = {255,192,128,64};
  // now, rescale
  double basevalue = 64;
  double pfac = (maxval-minval)/(256.0-basevalue);
  int *pixels  = new int[nx*ny];
  double ppfac(1.0/255.0);
  double pval;
  for (b=0; b<nb; b++) {
    for (a=0; a<na; a++) {
      // reset colormap
      for (k=0; k<3; k++)
	colmap[k] = 0;
      // read in mg pixel vals
      //sprintf(filename,"%s%d_%d_%d.mgs",fileprefix,b,a,ov);
      sprintf(filename,"%s%d_%d.mgs",fileprefix,b,ov);
      readObsImageMgs(filename,nx,ny,pixels);
      for (k=0; k<nmgs; k++) {
	modu = 0;
	for (i=0; i<ny; i+=stepsize) {
	  for (j=modu; j<nx; j+=stepsize) {
	    index = i*nx+j;
	    if (pixels[index] == k) {
	      pval  = MAX(basevalue,MIN(255.0,floor((dvals[b][a][index]-minval)/pfac + basevalue)));
	      //pval  = MAX(0.0,MIN(255.0,floor(255.0*(1.0-(dvals[b][a][index]-minval)/(maxval-minval)))));
	      im[index*3] = (char) (ppfac*colorshades[colmap[0]]*pval);
	      im[index*3+1] = (char) (ppfac*colorshades[colmap[1]]*pval);
	      im[index*3+2] = (char) (ppfac*colorshades[colmap[2]]*pval);
	    }
	  }
	  modu = (modu+1)%stepsize;
	}
	fillinimage(im,nx,ny,stepsize,eims[a]);
	//sprintf(filename,"%s%d_%d_%d.ppm",fileprefix,b,a,ov);
	sprintf(filename,"%s%d_%d.ppm",fileprefix,b,ov);
	writeImagePPM(filename,im,nx,ny);
	
	colmap[0]++;
	if (colmap[0] >= numcolorshades) {
	  colmap[0] = 0;
	  colmap[1]++;
	  if (colmap[1] >= numcolorshades) {
	    colmap[1] = 0;
	    colmap[2]++;
	    if (colmap[2] >= numcolorshades) 
	      colmap[2] = 0;
	  }
	}
      }
    }
  }

}
// fill in missing pixels from image
void fillinimage(char *im, int nx, int ny, int stepsize, char *eims)
{
  int i,j, k, ii, jj, iindex, index;
  int modu = 0, st;
  for (i=0; i<ny; i+=stepsize) {
    for (ii=0; ii<stepsize; ii++) {
      for (j=modu; j<nx; j+=stepsize) {
	index = i*nx+j;
	if (j==modu) 
	  st = -modu;
	else
	  st = 0;
	for (jj=st; jj<stepsize; jj++) {
	  iindex = (i+ii)*nx+(j+jj);
	  // original image (correct color) is at im[index]
	  // spread to new iindex
	  for (k=0; k<3; k++) 
	    im[iindex*3+k] = im[index*3+k];
	}
      }
    }
    modu = (modu+1)%stepsize;
  }
  // add edge image
  for (i=0; i<ny; i++) {
    for (j=0; j<nx; j++) {
      index= i*nx+j;
      if (eims[index] != 0) {
	im[index*3] = im[index*3+1] = im[index*3+2] = 255;
      }
    }
  }
}
