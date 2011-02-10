#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "myutil.h"
#include "prototypes.h"

extern int Nbins;
extern double *px,*pxt,*pstim1;
extern double ispike;
extern long idum;
extern int Ngrid;
extern int *max_ind;
extern int *max_ind_test;
extern double *MASTER_S;

extern double *pstim2,*pstim3,*pstim4;
extern double *xminN, *xmaxN, *stepN;
extern int *indN, *sizeN;

extern int size;

/* Calculates Info(vN + a * dvN) using fast method
 * Input - number of vectors, alpha, spike train, number of spikes, number of trials
 * Output - negative information
 * Note: pstim1 is size Ntrials * Nvec and contains vector projections
 * Note: pstim3 is size Ntrials and contains dvN projections
 * Note: pstim4 is size Ntrials and will contain projections of vN + a * dvN
 * Note: xminN, xmaxN, stepN, and indN are size Nvec
 * Note: sizeN is size Nvec + 1 and contains 1,Nbins,...,Nbins^Nvec */
double Nneginform1d_double(int Nvec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  unsigned long j, k;
  unsigned long ind;
  double neginform;

  for(j = 1; j <= Ntrials; j++) pstim4[j] = pstim1[(Nvec - 1)*Ntrials + j] + a * pstim3[j];

  xminN[Nvec] = min(pstim4,Ntrials);
  // cout<<"xminN="<<xminN[Nvec]<<endl;
  xmaxN[Nvec] = max(pstim4,Ntrials);
  stepN[Nvec] = (xmaxN[Nvec]-xminN[Nvec])/Nbins;

  for(j = 1; j <= Nbins; j++) {
    px[j] = 0.0;
    pxt[j] = 0.0;
  }

  for(j = 1; j <= Ntrials; j++) {
    indN[Nvec] = (pstim4[j] - xminN[Nvec])/stepN[Nvec];
    if(indN[Nvec] == Nbins) indN[Nvec]--;
    if((indN[Nvec] < 0)||(indN[Nvec] >= Nbins)) {
      cout<<j<<" "<<Nvec<<" "<<indN[Nvec]<<endl; 
      myerror("index out of bounds in Nneginform1d_double");
    }
    ind = indN[Nvec] + 1;
    px[ind] += 1.0;
    pxt[ind] += locator[j];
  }

  for(j = 1; j <= Nbins; j++) {
    px[j] /= Ntrials;
    pxt[j] /= Nspikes;
  }

  neginform = 0.0;

  for(j = 1; j <= Nbins; j++) {
    if(pxt[j] > 0.0) neginform -= pxt[j]*log(pxt[j]/px[j]);
  }

  return neginform / ispike;

}

/* Calculates information on test set
 * Input - number of vectors, vectors, test stimulus, size of vectors, test spike train, number of spikes, number of trials, projection vector
 * Output - negative information */
double Nneginform_double_test(int Nvec, double  *v,double *stimuli,unsigned long Nn, int *locator,unsigned long Nspikes, unsigned long Ntrials,double *pstim1_test)
{
  unsigned long i, j;
  unsigned long ind;
  double neginform;

  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1_test+i*Ntrials,v+i*Nn,stimuli,Nn,Ntrials);

  for(i = 1; i <= Nvec; i++) {
    xminN[i] = min(pstim1_test+Ntrials*(i-1),Ntrials);
    xmaxN[i] = max(pstim1_test+Ntrials*(i-1),Ntrials);
    stepN[i] = (xmaxN[i] - xminN[i])/Nbins;
  }

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    px[j] = 0.0;
    pxt[j] = 0.0;
  }

  for(i = 1; i <= Ntrials; i++) {
    for(j = 1; j <= Nvec; j++) {
      indN[j] = (pstim1_test[i+(j-1)*Ntrials]-xminN[j])/stepN[j];
      if(indN[j]==Nbins) indN[j]--;
    }
    ind = 1;
    for(j = 1; j <= Nvec; j++) ind += sizeN[j]*indN[j];
    px[ind] += 1.0;
    pxt[ind] += locator[i];
  }

  for(i = 1; i <= sizeN[Nvec+1]; i++) {
    px[i] /= Ntrials;
    pxt[i] /= Nspikes;
  }

  neginform = 0.0;

  for(i = 1; i <= sizeN[Nvec+1]; i++) {
    if(pxt[i] > 0.0) neginform -= pxt[i]*log(pxt[i]/px[i]);
  }

  return neginform / ispike;
}



/* Calculates Info(v1,...,vN + a * dvN) using fast method
 * Input - number of vectors, alpha, spike train, number of spikes, number of trials
 * Output - negative information
 * Note: pstim1 is size Ntrials * Nvec and contains vector projections
 * Note: pstim3 is size Ntrials and contains dvN projections
 * Note: pstim4 is size Ntrials and will contain projections of vN + a * dvN
 * Note: xminN, xmaxN, stepN, and indN are size Nvec
 * Note: sizeN is size Nvec + 1 and contains 1,Nbins,...,Nbins^Nvec */

double Nneginform_double(int Nvec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  unsigned long j, k;
  unsigned long ind;
  double neginform;

  for(j = 1; j <= Ntrials; j++) pstim4[j] = pstim1[(Nvec - 1)*Ntrials + j] + a * pstim3[j];

  for(j = 1; j < Nvec; j++) {
    xminN[j] = min(pstim1+(j-1)*Ntrials,Ntrials);
    xmaxN[j] = max(pstim1+(j-1)*Ntrials,Ntrials);
    stepN[j] = (xmaxN[j]-xminN[j])/Nbins;
  }

  xminN[Nvec] = min(pstim4,Ntrials);
  xmaxN[Nvec] = max(pstim4,Ntrials);
  stepN[Nvec] = (xmaxN[Nvec]-xminN[Nvec])/Nbins;

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    px[j] = 0.0;
    pxt[j] = 0.0;
  }

  for(j = 1; j <= Ntrials; j++) {
    for(k = 1; k < Nvec; k++) {
      indN[k] = (pstim1[j + Ntrials*(k-1)] - xminN[k])/stepN[k];
      if(indN[k] == Nbins) indN[k]--;
      if((indN[k] < 0)||(indN[k] >= Nbins)) {cout<<k<<" "<<indN[k]<<endl; myerror("index out of bounds in Nneginform_double");}
    }
    indN[Nvec] = (pstim4[j] - xminN[Nvec])/stepN[Nvec];
    if(indN[Nvec] == Nbins) indN[Nvec]--;
    if((indN[Nvec] < 0)||(indN[Nvec] >= Nbins)) {
      cout<<j<<" Nvec="<<Nvec<<" indN="<<indN[Nvec]<<" xmaxN="<<xmaxN[Nvec]<<" xminN="<<xminN[Nvec]<<endl; 
      myerror("index Nvec out of bounds in Nneginform_double");
    }
    ind = 1;
    for(k = 1; k <= Nvec; k++) ind += sizeN[k] * indN[k];
    px[ind] += 1.0;
    pxt[ind] += locator[j];
  }

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    px[j] /= Ntrials;
    pxt[j] /= Nspikes;
  }

  neginform = 0.0;

  for(j = 1; j <= sizeN[Nvec+1]; j++) {
    if(pxt[j] > 0.0) neginform -= pxt[j]*log(pxt[j]/px[j]);
  }

  return neginform / ispike;
}




double neginform_double(double  *v,double *stimuli,unsigned long Nn, int *locator,unsigned long Nspikes, unsigned long Ntrials)
{
  int bin_ind;
  unsigned long j;
  double stepp,xmin;
  double xmax,neginfrm=0;

  Proj1d_double(pstim1,v,stimuli,Nn,Ntrials);
  //cout<<"elements inside the neginform_double"<<v[1]<<" "<<v[2]<<endl;
  xmax=max(pstim1,Ntrials);
  xmin=min(pstim1,Ntrials);

  // stepp=(xmax-xmin)/(Nbins-1);
  //  cout<<"inside negin: "<<xmin<<" "<<xmax<<endl;
  if (xmin==xmax) {	
    cout<<"vector inside neginform_double "<<v[1]<<" "<<v[Nn]<<" "<<length(v,Nn)<<endl;
    cout<<"stimuli inside neginform_double "<<stimuli[1]<<" "<<stimuli[Nn]<<" "<<length(stimuli,Nn)<<endl;
    cout<<"xmin="<<xmin<<" "<<xmax<<endl; 
    myerror("zero width range in neginform_double");
  }
  stepp=(xmax-xmin)/Nbins;
  // step1=stepp;
  //xmin1=xmin;
  for(j=1;j<=Nbins;j++)  { px[j]=0; pxt[j]=0;}
  for(j=1;j<=Ntrials;j++){
    pstim1[j]-=xmin;
    pstim1[j]/=stepp;
    bin_ind=(int)(pstim1[j])+1;
    if (bin_ind==Nbins+1) bin_ind--;
    if ((bin_ind <1) ||(bin_ind >Nbins)) { 
      cout<<"index is out of bounds in neginform"<<endl;
      cout<<bin_ind<<" "<<pstim1[j]<<" "<<xmin<<" "<<xmax<<" "<<stepp<<" "<<(pstim1[j]-xmin)/stepp;
      cout<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[Nn]<<endl;
      //      myerror(" index out of bounds in neginform_double");
      cout<<" index out of bounds in neginform_double"<<endl;
      return(10);
    }
    px[bin_ind]+=1.0/(double)Ntrials;  
    pxt[bin_ind]+=(double)locator[j]/(double)Nspikes;
  }
 
  for(j=1;j<=Nbins;j++){
    if (pxt[j]>0)  neginfrm-=pxt[j]*log(pxt[j]/px[j]);
  }
 
  return(neginfrm/ispike);
}

