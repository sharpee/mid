#include <fstream>
#include <iostream>
using namespace std;
#include "nrutil.h"
#include "patch.h"
#include "myutil.h"


void runChecks( FILE * inf, unsigned long Ntrials, int dimx, int dimy, int x0, int y0, int Nh, int Nv, int cx, int cy, int bytesPerData) 
{
	if(Nh%cx !=0)
	  myerror("x downsample does not divide evenly into sta width.");
	if(Nv%cy !=0)
	  myerror("y downsample does not divide evenly into sta height.");
	if(Nh + x0 > dimx+1) 
	  myerror("x offset plus sta width should be less than or equal to width + 1.");
	if(Nv + y0 > dimy +1) 
	  myerror("y offset plus sta height should be less than or equal to height + 1.");
	if(x0 < 1 || x0 > Nh)
	  myerror("x offset should be between 1 and Nh.");
	if(y0 < 1 || y0 > Nh)
	  myerror("y offset should be between 1 and Nv.");
	long expectedFileSize = Ntrials*dimx*dimy*bytesPerData;
	
    fseek(inf, 0L, SEEK_END);
    long fileSize = ftell(inf);
    fseek(inf, 0L, SEEK_SET);
	 
	if(expectedFileSize != fileSize) {
	   cout << "Nh: " << Nh << " Nv:" << Nv << " Ntrials: " << 
	     Ntrials << " expectedFileSize: " << expectedFileSize << 
		 " fileSize: " << fileSize <<endl ;
	   myerror("expectedFileSize does not equal fileSize");
	}
}

//reads as unsigned bytes, rescales (stim - 128)/255, and outputs as double
void readinframes_to_double(char *movieFile,double *stimuli, unsigned long Ntrials,int dimx, int dimy, int x0, int y0, int Nh, int Nv, int cx, int cy)
{
  int nmax,i,j,Ntcurrent=0;
  unsigned long k;
  int fsize=(int)((Nv*Nh)/(cx*cy));
  double sum;
 
  unsigned char *b;
  float *bfloat;
  unsigned char *tstpatch;
  char buffer[256];
 
  tstpatch=cvector(0,fsize-1);
  b=cvector(0,dimx*dimy-1);
  bfloat=vector(0,dimx*dimy-1);
 
  FILE *inf;

  Patch ptc(dimx, dimy, x0-1, y0-1, Nh, Nv, cx, cy);
 

  inf=fopen(movieFile,"rb");
	runChecks(inf, Ntrials, dimx, dimy, x0, y0, Nh, Nv, cx, cy, sizeof(char));
    if (inf ==NULL)      erroropen(movieFile);  
    for(i=1;i<=Ntrials;i++){
      fread(b,sizeof(unsigned char),dimx*dimy,inf);
      if (ferror(inf)) myerror("error reading from movie file.");
      ptc.Convert(b,tstpatch);		 
      for(j=1;j<=fsize;j++){
	stimuli[(i-1)*fsize+j]= (double) (tstpatch[j-1]-128)/255.0;
      }
    }
    fclose(inf); 

  free_vector(bfloat,0,dimx*dimy-1); 
  free_cvector(b,0,dimx*dimy-1); 
  free_cvector(tstpatch,0,fsize-1);
}

void readinframes_as_double(char *movieFile,double *stimuli, unsigned long Ntrials,int dimx, int dimy, int x0, int y0, int Nh, int Nv, int cx, int cy)
{
  int nmax,i,j,Ntcurrent=0;
  unsigned long k;
  int fsize=(int)((Nv*Nh)/(cx*cy));
  double sum;
  void erroropen(char *);

  double *b;
  double *tstpatch;
 
  tstpatch=dvector(0,fsize-1);
  b=dvector(0,dimx*dimy-1);
 
  FILE *inf;

  Patch ptc(dimx, dimy, x0-1, y0-1, Nh, Nv, cx, cy);
 
    inf=fopen(movieFile,"rb");
	
	runChecks(inf, Ntrials, dimx, dimy, x0, y0, Nh, Nv, cx, cy, sizeof(double));
 
    if (inf ==NULL)  erroropen(movieFile);  
    for(i=1;i<=Ntrials;i++){
      fread(b,sizeof(double),dimx*dimy,inf);
      if (ferror(inf)) myerror("error reading from movie file.");
      ptc.dConvert(b,tstpatch);		 
      for(j=1;j<=fsize;j++){
	    stimuli[(i-1)*fsize+j]= tstpatch[j-1];
      }
    }
    fclose(inf); 
  free_dvector(b,0,dimx*dimy-1); 
  free_dvector(tstpatch,0,fsize-1);
}

