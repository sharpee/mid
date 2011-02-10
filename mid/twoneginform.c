using namespace std;
#include <iostream>
#include <math.h>
#include <cstdlib>
#include "myutil.h"
#include "prototypes.h"

extern int size,Nbins;
extern double *pstim2,*pstim1,*px,*pxt;
extern double ispike;
extern int Ngrid;
extern int *max_ind;
extern int *max_ind_test;


double twoneginform_double(double  *v,double *stimuli,unsigned long n, int *locator,unsigned long Nspikes,unsigned long Ntrials)
{
  int ind1,ind2;
  unsigned long j;
  double step1,xmin1,xmax1,step2,xmin2,xmax2,twoneginfrm=0;  
  unsigned long Nn=n>>1;

  Proj1d_double(pstim1,v,stimuli,Nn,Ntrials);
  xmax1=max(pstim1,Ntrials);
  xmin1=min(pstim1,Ntrials);
  step1=(xmax1-xmin1)/Nbins;

  Proj1d_double(pstim2,v+Nn,stimuli,Nn,Ntrials); 
  xmax2=max(pstim2,Ntrials);  
  xmin2=min(pstim2,Ntrials);
  step2=(xmax2-xmin2)/Nbins;

  for(j=1;j<=size;j++){ px[j]=0; pxt[j]=0;}
  
  for(j=1;j<=Ntrials;j++){
    pstim1[j]-=xmin1;
    pstim1[j]/=step1;
    ind1=(int)(pstim1[j])+1;
    if (ind1==Nbins+1) ind1--;
    if  ((ind1 <1) || (ind1> Nbins)) {
      cout<<ind1<<" xmin1="<<xmin1<<" j="<<j<<" pstim1[j]"<<pstim1[j]<<" step1 "<<step1<<endl;
      myerror("index1 <1 in twoneginform ");
    }
    
    ind2=(int)((pstim2[j]-xmin2)/step2);
    if (ind2==Nbins) ind2--;
    if ((ind2 <0) || (ind2 >=Nbins)) { 
      cout<<ind2<<" xmin2="<<xmin2<<endl; myerror("index2 out of bounds twoneginform_double ");
    }
    
    px[ind2*Nbins+ind1]+=1.0/(double)Ntrials; 
    pxt[ind2*Nbins+ind1]+=(double)locator[j]/(double)Nspikes; 
  }
  
  for(j=1;j<=size;j++){
    if (pxt[j]>0) twoneginfrm-=pxt[j]*log(pxt[j]/px[j]);
  }
  return(twoneginfrm/ispike);
}

