#include <math.h>
#include <iostream>
extern unsigned long fsize;
extern int Ngrid;
extern int nlags;
extern int LAG_SHIFT;
extern double xmin1;
extern double step1;
extern double *px_special,*pxt_special;
extern long idum;
extern int Nbins;
extern int *locator;


void Proj1d_double(double *p1,double *v,double *stimuli,unsigned long Nn, unsigned long Ntrials)
{ 
   int  n;
   unsigned long i;
   double dSum;
   double * pdVector;
   double * pStimulus;
   int chunk=Ntrials/2;
   
   for(n=0; n < Ntrials; n++){
       dSum = 0;
       pdVector = & v[1];
       pStimulus = & stimuli[n*fsize + 1];
       for(i=1; i <= Nn; i++){
		   dSum += *pdVector++ * *pStimulus++;
       }
       p1[n+1] = dSum;
   }
}

