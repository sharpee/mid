#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "myutil.h"
#include "prototypes.h"
#define EPS 1.0e-10
extern unsigned long fsize;
extern int Nbins;
extern long idum;
extern double *sx,*sxt,*px,*pxt,*pstim1,ispike;
extern double *c;//*fx;
extern int *max_ind;
extern int Ngrid;
extern int Ngrid_linear;
extern int *transl;
extern int nlags;

extern double *xminN, *xmaxN, *stepN;
extern int *indN, *sizeN;
extern double *pstim3,*pstim4;
extern int size;


/* Calculates dInfo(v1,...,vN+a*dv)/dvN
 * Input - number of vectors, alpha, v, vN+a*dv,output derivative, stimulus, size of each vector, spike train, number of spikes, number of trials
 * Note: pstim1 is size Ntrials*Nvec and contains vector projections
 * Note: pstim3 is size Ntrials and contains dv projections
 * Note: pstim4 is size Ntrials and will contain projections of vN + a * dv
 * Note: xminN, xmaxN, stepN, and indN are size Nvec
 * Note: sizeN is size Nvec + 1 and contains 1,Nbins,...,Nbins^Nvec */
void Nderivative_double(int Nvec, double a, double *v, double *vN, double *dv, double *stimuli, unsigned long Nn, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  unsigned long j, k, l;
  unsigned long ind;
  char found_ok;
  double temp;

  for(j = 1; j <= Ntrials; j++) pstim4[j] = pstim1[(Nvec-1)*Ntrials +j] + a * pstim3[j];

  for(j = 1; j < Nvec; j++) {
    xminN[j] = min(pstim1+(j-1)*Ntrials,Ntrials);
    xmaxN[j] = max(pstim1+(j-1)*Ntrials,Ntrials);
    stepN[j] = (xmaxN[j] - xminN[j])/Nbins;
  }

  xminN[Nvec] = min(pstim4,Ntrials);
  xmaxN[Nvec] = max(pstim4,Ntrials);
  stepN[Nvec] = (xmaxN[Nvec] - xminN[Nvec])/Nbins;

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    px[j] = 0.0;
    pxt[j] = 0.0;
    for(k = 1; k <= Nn; k++) {
      sx[(j-1)*Nn + k] = 0.0;
      sxt[(j-1)*Nn + k] = 0.0;
    }
  }

  for(j = 1; j <= Ntrials; j++) {
    for(k = 1; k < Nvec; k++) {
      indN[k] = (pstim1[j + Ntrials*(k-1)] - xminN[k])/stepN[k];
      if(indN[k] == Nbins) indN[k]--;
      if((indN[k] < 0)||(indN[k] >= Nbins)) myerror("index out of bounds in Nderivative_double");
    }
    indN[Nvec] = (pstim4[j] - xminN[Nvec])/stepN[Nvec];
    if(indN[Nvec] == Nbins) indN[Nvec]--;
    if((indN[Nvec] < 0)||(indN[Nvec] >= Nbins)) myerror("index out of bounds in Nderivative_double");
    ind = 0;
    for(k = 1; k <= Nvec; k++) ind += sizeN[k] * indN[k];
    px[ind + 1] += 1.0;
    pxt[ind + 1] += locator[j];
    for(k = 1; k <= Nn; k++) {
      sx[ind*Nn+k] += stimuli[(j-1)*fsize+k];
      sxt[ind*Nn+k] += stimuli[(j-1)*fsize+k]*locator[j];
    }
  }

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    if(px[j] > 0) {
      pxt[j] = (pxt[j]/Nspikes)/(px[j]/Ntrials);
      if(pxt[j] > 0) {
	for(k = 1; k <= Nn; k++) {
	  sx[(j-1)*Nn+k] /= Ntrials;
	  sx[(j-1)*Nn+k] -= sxt[(j-1)*Nn+k]/Nspikes/pxt[j];
	}
      }
      if(pxt[j] == 0) {
	for(k = 1; k <= Nn; k++) {
	  sx[(j-1)*Nn+k] /= Ntrials;
	}
      }
    }
  }

  for(k = 1; k <= Nn; k++) dv[k] = 0.0;

  found_ok = 0;
  for(j = 1; j <= sizeN[Nvec]; j++) {
    for(k = 2; k < Nbins - 2; k++) {
      if((px[j+(k-2)*sizeN[Nvec]]>0)&&(px[j+(k-1)*sizeN[Nvec]]>0)&&(px[j+k*sizeN[Nvec]]>0)&&(px[j+(k+1)*sizeN[Nvec]]>0)&&(px[j+(k+2)*sizeN[Nvec]]>0)) {
	temp=(c[1]*pxt[j+(k-2)*sizeN[Nvec]]+c[2]*pxt[j+(k-1)*sizeN[Nvec]]+c[4]*pxt[j+(k+1)*sizeN[Nvec]]+c[5]*pxt[j+(k+2)*sizeN[Nvec]])/stepN[Nvec];
	if(!found_ok) found_ok = 1;
	if(temp != 0 ) {
	  for(l = 1; l <= Nn; l++) dv[l] += sx[(k*sizeN[Nvec]+j-1)*Nn+l]*temp;
	}
      }
    }
  }

  if(!found_ok) {
    for(j = 1; j <= sizeN[Nvec]; j++) {
      for(k = 2; k < Nbins - 2; k++) {
	if(pxt[j+k*sizeN[Nvec]]>0) {
	  temp=(c[1]*pxt[j+(k-2)*sizeN[Nvec]]+c[2]*pxt[j+(k-1)*sizeN[Nvec]]+c[4]*pxt[j+(k+1)*sizeN[Nvec]]+c[5]*pxt[j+(k+2)*sizeN[Nvec]])/stepN[Nvec];
	  if(temp != 0) {
	    for(l = 1; l <= Nn; l++) dv[l] += sx[(k*sizeN[Nvec]+j-1)*Nn+l]*temp;
	  }
	}
      }
    }
  }

  temp = length(dv,Nn);
  if(temp == 0) {cout<<"zero derivative length"<<endl;
  for(k = 1; k <= Nn; k ++) dv[k] = 0.1*length(vN,Nn)*(ran2(&idum)-0.5);
  }
  /*
  for(j = 0; j < Nvec - 1; j++) {
    temp = dot(dv,v+Nn*j,Nn)/dot(v+Nn*j,v+Nn*j,Nn);
    for(k = 1; k <= Nn; k++) dv[k] -= v[k+j*Nn] * temp;
  }
  */

  temp = dot(dv,vN,Nn)/dot(vN,vN,Nn);
  for(k = 1; k <= Nn; k++) dv[k] -= vN[k] * temp;
}

/* Calculates dInfo(vN+a*dv)/dvN and makes it normal to v
 * Input - number of vectors, alpha, v, vN+a*dv,output derivative, stimulus, size of each vector, spike train, number of spikes, number of trials
 * Note: pstim1 is size Ntrials*Nvec and contains vector projections
 * Note: pstim3 is size Ntrials and contains dv projections
 * Note: pstim4 is size Ntrials and will contain projections of vN + a * dv
 * Note: xminN, xmaxN, stepN, and indN are size Nvec
 * Note: sizeN is size Nvec + 1 and contains 1,Nbins,...,Nbins^Nvec */
void Nderivative1d_double(int Nvec, double a, double *v, double *vN, double *dv, double *stimuli, unsigned long Nn, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  unsigned long j, k, l;
  unsigned long ind;
  char found_ok;
  double temp;

  if(a==0.0) {
    for(j = 1; j <= Ntrials; j++) pstim4[j] = pstim1[(Nvec-1)*Ntrials + j];
  } else {
    for(j = 1; j <= Ntrials; j++) pstim4[j] = pstim1[(Nvec-1)*Ntrials + j] + a * pstim3[j];
  }

  xminN[Nvec] = min(pstim4,Ntrials);
  xmaxN[Nvec] = max(pstim4,Ntrials);
  stepN[Nvec] = (xmaxN[Nvec] - xminN[Nvec])/Nbins;

  for(j = 1; j <= Nbins; j++) {
    px[j] = 0.0;
    pxt[j] = 0.0;
    for(k = 1; k <= Nn; k++) {
      sx[(j-1)*Nn + k] = 0.0;
      sxt[(j-1)*Nn + k] = 0.0;
    }
  }

  for(j = 1; j <= Ntrials; j++) {
    indN[Nvec] = (pstim4[j] - xminN[Nvec])/stepN[Nvec];
    if(indN[Nvec] == Nbins) indN[Nvec]--;
    if((indN[Nvec] < 0)||(indN[Nvec] >= Nbins)) myerror("index out of bounds in Nderivative_double");
    ind = indN[Nvec];
    px[ind + 1] += 1.0;
    pxt[ind + 1] += locator[j];
    for(k = 1; k <= Nn; k++) {
      sx[ind*Nn+k] += stimuli[(j-1)*fsize+k];
      sxt[ind*Nn+k] += stimuli[(j-1)*fsize+k]*locator[j];
    }
  }

  for(j = 1; j <= Nbins; j++) {
    if(px[j] > 0) {
      pxt[j] = (pxt[j]/Nspikes)/(px[j]/Ntrials);
      if(pxt[j] > 0) {
	for(k = 1; k <= Nn; k++) {
	  sx[(j-1)*Nn+k] /= Ntrials;
	  sx[(j-1)*Nn+k] -= sxt[(j-1)*Nn+k]/Nspikes/pxt[j];
	}
      }
      if(pxt[j] == 0) {
	for(k = 1; k <= Nn; k++) {
	  sx[(j-1)*Nn+k] /= Ntrials;
	}
      }
    }
  }

  for(k = 1; k <= Nn; k++) dv[k] = 0.0;

  found_ok = 0;
  for(j = 3; j < Nbins - 1; j++) {
    if((px[j-2]>0)&&(px[j-1]>0)&&(px[j]>0)&&(px[j+1]>0)&&(px[j+2]>0)) {
      temp=(c[1]*pxt[j-2]+c[2]*pxt[j-1]+c[4]*pxt[j+1]+c[5]*pxt[j+2])/stepN[Nvec];
      if(!found_ok) found_ok = 1;
      if(temp != 0 ) {
	for(l = 1; l <= Nn; l++) dv[l] += sx[(j-1)*Nn+l]*temp;
      }
    }
  }

  if(!found_ok) {
    for(j = 3; j < Nbins - 1; j++) {
      if(pxt[j]>0) {
	temp=(c[1]*pxt[j-2]+c[2]*pxt[j-1]+c[4]*pxt[j+1]+c[5]*pxt[j+2])/stepN[Nvec];
	if(temp != 0) {
	  for(l = 1; l <= Nn; l++) dv[l] += sx[(j-1)*Nn+l]*temp;
	}
      }
    }
  }

  temp = length(dv,Nn);
  if(temp == 0) for(k = 1; k <= Nn; k ++) dv[k] = 0.1*length(vN,Nn)*(ran2(&idum)-0.5);

  for(j = 0; j < Nvec - 1; j++) {
    temp = dot(dv,v+Nn*j,Nn)/dot(v+Nn*j,v+Nn*j,Nn);
    for(k = 1; k <= Nn; k++) dv[k] -= v[k+j*Nn] * temp;
  }

  temp = dot(dv,vN,Nn)/dot(vN,vN,Nn);
  for(k = 1; k <= Nn; k++) dv[k] -= vN[k] * temp;
}

