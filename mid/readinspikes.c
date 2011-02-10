using namespace std;
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include "nrutil.h"

int countSpikeFrames(char * spikefile) 
{
  ifstream inf;
  FILE *fp;
 
  char buffer[256];
  fp=fopen(spikefile,"r");
  if (fp==NULL) erroropen(spikefile);
  cout<<"reading spikes from "<<spikefile<<endl;
  
  unsigned long count = 0;
  while(!ferror(fp) && (fgets(buffer,30,fp) != NULL)) count++;
  return count;
  fclose(fp);
}

void readinspikes(char * spikefile, int *repetition,int *locator,int Ntrials,int nlags)
{
  unsigned long i,j;
  unsigned long sum(int *,unsigned long);
  ifstream inf;
  FILE *fp;
  int k;
  int N_violators=0;
  char buffer[256];

  fp=fopen(spikefile,"r");
  if (fp==NULL) erroropen(spikefile);

  for(i=1;i<=Ntrials;i++){
     if (! ferror(fp)) {
       k=atoi(fgets(buffer,30,fp));
       if (i>=nlags) locator[i-nlags+1]+=k;
     }
   }

   fclose(fp);
}

