using namespace std;
#include <iostream>
#include <fstream>
#include <math.h>
#include "nrutil.h"
#include "myutil.h"
#include "prototypes.h"


extern double *pstim1, *pstim3;

double chooseN_double(int Nvec, double *v1,double *v2,double *stimuli,unsigned long n, int *locator,unsigned long Nspikes,unsigned long Ntrials,double (*func)(int Nvec, double a, int *locator,unsigned long Nspikes,unsigned long Ntrials))
{
  double temp,temp1,temp2;
  double *temp_vector;
  double pstim_temp;
  temp_vector=dvector(1,n);
  unsigned long i,j;
  double imax;
  unsigned long nV = (Nvec-1)*n;
  unsigned long nT = (Nvec-1)*Ntrials;

 //initial normalization 
  temp2=length(v2+nV,n); 
  temp1=length(v1+nV,n); 
  for(i=1+nV;i<=n*Nvec;i++) {
    v2[i]/=temp2;
    v1[i]/=temp1;
  }

  temp=dot(v1+nV,v2+nV,n);

  if (fabs(temp)==1) {
    cout<<"vectors are identical, no rotation is needed"<<endl;
    return(imax);
  }

  for(i=1+nV;i<=n*Nvec;i++) v2[i] -= temp*v1[i];
  temp=length(v2+nV,n); 
  cout<<"temp "<<temp<<endl;
  for(i=1+nV;i<=n*Nvec;i++) v2[i]/=temp;

  Proj1d_double(pstim1+nT,v1+nV,stimuli,n,Ntrials);
  Proj1d_double(pstim3,v2+nV,stimuli,n,Ntrials);

  cout<<"before imax"<<endl;
  cout<<"v1 "<<v1[1+nV]<<" v2 "<<v2[1+nV]<<endl;
  cout<<"ps1 "<<pstim1[1+nT]<<" ps3 "<<pstim3[1]<<endl;

  imax=-(*func)(Nvec,0.0,locator,Nspikes,Ntrials);

  for(i=0;i<=19;i++){
    if (i<=10) temp1=-(*func)(Nvec,-1.0+(double)i/5.0,locator,Nspikes,Ntrials);
    else temp1=-(*func)(Nvec,3.0-(double)i/5.0,locator,Nspikes,Ntrials);

   //switch v1 and v2 at i=10
    if(i==10) {
      for(j = 1; j <= Ntrials; j++) {
	pstim_temp = pstim1[j+nT];
	pstim1[j+nT] = pstim3[j];
	pstim3[j] = pstim_temp;
      }
    }
    if (temp1>=imax) {
      imax=temp1;
      temp2=i;
    }
  }

 //switch v1 and v2 back
  for(j = 1; j <= Ntrials; j++) {
    pstim_temp = pstim1[j+nT];
    pstim1[j+nT] = pstim3[j];
    pstim3[j] = pstim_temp;
  }

  for(j=1;j<=n;j++) {
    if (temp2<=10){
      v1[j+nV]+=(-1.0+temp2/5.0)*v2[j+nV];
    }
    else{
      v1[j+nV]=(3.0-temp2/5.0)*v1[j+nV]+v2[j+nV];
    }
  }

 //final normalization
  temp2=length(v2+nV,n); 
  temp1=length(v1+nV,n); 
  for(i=1;i<=n;i++) {
    v2[i+nV]/=temp2;
    v1[i+nV]/= temp1;
  }
  free_dvector(temp_vector,1,n);
  return(imax); 
}

