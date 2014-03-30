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
void addedges(unsigned char *im, int nx, int ny, unsigned char *eims);

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
void readEdgeImagePGM(char *filename, int nx, int ny, unsigned char *eimage)
{
  FILE *stream = fopen(filename, "r" );
  if ( !stream )
    return;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(eimage, sizeof(unsigned char), ny*nx, stream);
  fclose(stream);
}

void writeImagePPM(char *filename, unsigned char *im, int nx, int ny)
{
  FILE *stream = fopen(filename, "w" );
  if ( !stream )
    return;
  fprintf( stream, "P6\n%d %d\n255\n", nx, ny);
  fwrite(im, sizeof(unsigned char), 3*ny*nx, stream);
  fclose(stream);
}

// usage: 
// rescaleims  fileprefix (e.g. obsreg34) numalphas numbeliefs numactions
int main(int argc, char *argv[])
{
  char *fileprefix, *edgeimprefix;
  char filename[256];
  char *uimname;
  // prefixes
  uimname = *++argv;

  int nmgs = atoi(*++argv);
  int nt = atoi(*++argv);
  int na = atoi(*++argv);
  int ov = atoi(*++argv);

  int stepsize = 1;
  int modu = 0;

  int i,j, index,a,b,t, k;
  double ***dvals;
  double maxval, minval;
  unsigned char *im;
  int nx, ny;
  nx = 160;
  ny = 120;

  // read in edge images
  unsigned char **eims = new unsigned char*[na];
  for (a=0; a<na; a++) {
    eims[a] = new unsigned char[nx*ny];
    sprintf(filename,"edges%d.pgm",a);
    readEdgeImagePGM(filename, nx, ny, eims[a]);
  }

  
  im = new unsigned char[nx*ny*3];
  unsigned char *uim = new unsigned char[nx*ny*3];
  FILE *stream =  fopen(uimname,"r");
  if ( !stream )
    return 0;
  fscanf( stream, "%*s\n%*d %*d\n%*d\n");
  fread(uim, sizeof(unsigned char), 3*ny*nx, stream);
  fclose(stream);

  int colmap[3];
  int numcolorshades = 4;
  int colorshades[4] = {255,192,128,64};
  // now, rescale
  int *pixels  = new int[nx*ny];
  for (t=0; t<nt; t++) {
    // reset colormap
    for (k=0; k<3; k++)
      colmap[k] = 0;
    // read in mg pixel vals
    sprintf(filename,"sim%dobsreg0_%d.mgs",t,ov);
    fprintf(stderr,"reading %s\n",filename);
    readObsImageMgs(filename,nx,ny,pixels);
    for (k=0; k<nmgs; k++) {
      modu = 0;
      for (i=0; i<ny; i+=stepsize) {
	for (j=modu; j<nx; j+=stepsize) {
	  index = i*nx+j;
	  if (pixels[index] == k) {
	    im[index*3] = (unsigned char) (colorshades[colmap[0]]);
	    im[index*3+1] = (unsigned char) (colorshades[colmap[1]]);
	    im[index*3+2] = (unsigned char) (colorshades[colmap[2]]);
	  }
	}
	modu = (modu+1)%stepsize;
      }
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
    // assumes they're always the same
    addedges(im,nx,ny,eims[0]);
    sprintf(filename,"simobsreg%d.ppm",t);
    writeImagePPM(filename,im,nx,ny);
  }
}
void addedges(unsigned char *im, int nx, int ny, unsigned char *eims)
{
  int i,j,index;
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
