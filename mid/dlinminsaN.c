#include <math.h>
#include "myutil.h"
#include "prototypes.h"

double (*nrfuncdN)(int Nvec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials);
void (*nrdfundN)(int Nvec, double a, double [], double [], double [], double *stimuli, unsigned long n, int *locator, unsigned long Nspikes, unsigned long Ntrials);


void NdlinminsaN_double(int Nvec, double p[], double xi[], unsigned long n, int type, double *stimuli, int *locator, unsigned long Nspikes, unsigned long Ntrials, double *fret, double *pT, double (*func)(int Nvec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials), void (*dfunc)(int Nvec, double a, double *v, double *vN, double *dv, double *stimuli, unsigned long n, int *locator, unsigned long Nspikes, unsigned long Ntrials)) 
{
  unsigned long j;
  double xx, xmin, fx, fb, fa, bx, ax;
  double temp;

  nrfuncdN=func;
  nrdfundN=dfunc;

  ax = 0.0;
  xx = length(p+(Nvec-1)*n,n);
  xx /= length(xi,n);
  temp = xx;
  xx *= -2.0;

  NmnbrakN(Nvec,&ax,&xx,&bx,&fa,&fx,&fb,Nf1dimN,p,xi,n,stimuli,locator,Nspikes, Ntrials);
  *fret=NdbrentsaN_double(Nvec,ax,xx,bx,pT,Nf1dimN,Ndf1dimN,2.0e-4,&xmin,p,xi,n,stimuli,locator,Nspikes,Ntrials);
  if(xmin==0) (*pT)=0;
  for(j = 1; j <= n; j++) {
    xi[j] *= xmin;
    p[j+(Nvec-1)*n] += xi[j];
  }
}

