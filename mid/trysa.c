#include <iostream>
#include <fstream>
using namespace std;
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "myutil.h"
#include "prototypes.h"

extern double *pstim1, *pstim2, *pstim3;


void Ntrysa_double(int Nvec,double p[], double pbest[], double xi[], unsigned long n, int type, double *stimuli, int *locator, unsigned long Nspikes, unsigned long Ntrials, double ftol, double *pT, int ITMAX, int *count, double UPDATEFACTOR, double *negmin, char *filename1, double Tstart, double tdown, double tup, double tstopup, int *iter, double *fret, double *ave_training, double (*func)(int Nvec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials), void (*dfunc)(int Nvec, double a, double *v, double *vN, double *dv, double *stimuli, unsigned long Nn, int *locator, unsigned long Nspikes, unsigned long Ntrials)) 
{
  int its;
  double fp, temp;
  unsigned long k;
  *ave_training = 0.0;
  FILE *fpointer;
  char localbuffer[256];
  int nV;
  unsigned long nT;

  nV = (Nvec - 1) * n;
  nT = (Nvec - 1) * Ntrials;

  fp=(*func)(Nvec,0.0,locator,Nspikes,Ntrials);
  if(fp==10) {
    *fret=0;
    return;
  }

  for(its = 1; its <= ITMAX; its++) {
    *iter = its;
    temp = length(p+nV,n);
    if (temp > 50000) {
      for(k = 1; k <= n; k++) p[k+nV] /= temp;
      for(k=1;k<=Ntrials;k++) pstim1[k+nT]/=temp;
    }
    if(temp==0) myerror("zero vector in Ntrysa_double");

    //Calculate derivative with respect to vN
    (*dfunc)(Nvec,0.0,p,p+(Nvec-1)*n,xi,stimuli,n,locator,Nspikes,Ntrials);
    Proj1d_double(pstim3,xi,stimuli,n,Ntrials);

    //Find optimal step size
    NdlinminsaN_double(Nvec,p,xi,n,type,stimuli,locator,Nspikes,Ntrials,fret,pT,func,dfunc);

    if(*fret==10) {
      *fret=its;
      return;
    }

    (*pT) *= tdown;
    if(*pT < tstopup) (*pT) = tstopup;
    if(*fret < *negmin) {
      *negmin = *fret;
      *count = 0;
      temp = length(p+nV,n);
      for(k = 1 + nV; k <= Nvec*n; k++) pbest[k]  = p[k] / temp;
    }
    else {
      (*count)++;
      //if vector is very bad, change to a linear combination of p and pbest
      if(*fret - *negmin > UPDATEFACTOR * (*pT)) {
	cout<<"attempting to change vectors in trysa"<<endl;
	chooseN_double(Nvec,pbest,p,stimuli,n,locator,Nspikes,Ntrials,func);
	for(k = 1 + nV; k <= n*Nvec; k++) p[k] = pbest[k];
	Proj1d_double(pstim1+nT,p+nV,stimuli,n,Ntrials);
	fpointer=fopen(filename1,"a");
	if(fpointer==NULL) myerror(filename1);
	fprintf(fpointer,"#changing vector to the best-so-far one from %f, and I thought it to = %f\n",-(*fret),-(*negmin));
	*negmin=*fret=(*func)(Nvec,0.0,locator,Nspikes,Ntrials);
	if(*fret==10){
	  *fret=ITMAX+1;
	  return;
	}
	fprintf(fpointer,"and it is = %f (could be larger because of choose_double)\n",-(*fret));
	fclose(fpointer);
	cout<<"successfully changed vector inside trysa"<<endl;
      }
    }
    if(2.0*fabs((*fret)-fp)<=ftol*(fabs(*fret)+fabs(fp))) {
      cout<<"getting stuck inside trysa"<<endl;
      Proj1d_double(pstim1+nT,p+nV,stimuli,n,Ntrials);
      (*dfunc)(Nvec,0.0,p,p+nV,xi,stimuli,n,locator,Nspikes,Ntrials);
      temp=-(*pT)*length(p+nV,n)/length(xi,n);
      for(k = 1; k <= n; k++) {
	xi[k] *= temp;
	p[k+nV]+= xi[k];
      }
      Proj1d_double(pstim1+nT,p+nV,stimuli,n,Ntrials);
      *fret=(*func)(Nvec,0.0,locator,Nspikes,Ntrials);
      if(*fret==10) {
	*fret=ITMAX+2;
	return;
      }
      cout<<"correction was necessary: new xmin"<<-(*pT)<<" and -I="<<(*fret)<<endl;
      //increating temperature
      if(tup > 1.0) {
	if(*pT < Tstart) (*pT) *= tup;
	if(*pT > Tstart) (*pT) = Tstart;
      }
    }
    fp = (*fret);
    (*ave_training) -= fp;
  }
  (*ave_training) /= (*iter);
  return;
}

