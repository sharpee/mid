using namespace std;
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include "nrutil.h"
#include "myutil.h"
#include "prototypes.h"

extern double *pstim1,*pstim2;

/* my standard error handler */
void myerror(const char *error_text)
{
  fprintf(stderr,"my run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/* my standard error handler */
void erroropen(char error_text[])
{
  fprintf(stderr,"error opening file...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

void putzeros(char * buffer)
{
  int i=0;
  while((buffer[i]!='\0')) {
    if (buffer[i]==' ') buffer[i]='0';
    i++;
  }
}

unsigned long sum(int v[], unsigned long n)
{
  unsigned long s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=v[i];};
  return(s);
}

unsigned long csum(unsigned char v[], unsigned long n)
{
  unsigned long s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=(int)v[i];};
  return(s);
}

double dsum(double v[], unsigned long n)
{
  double s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=v[i];};
  return(s);
}

double mean(double v[],unsigned long n)
{
  double s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=v[i];}
  return(s/(double)n);
}

double infosmall(int *locator, unsigned long ntrials)
{
  unsigned long nspikes=sum(locator,ntrials);
  unsigned long i;
  double s1;
  if (nspikes>0){
    s1= -log((double)nspikes/(double)ntrials)/log(2);
    for(i=1;i<=ntrials;i++) {
      if (locator[i]>0) s1+=(double)locator[i]*log((double)locator[i])/log(2)/(double)nspikes;
    }
  }
  else s1=-100;
  return(s1); 
}

double varsmall(int *locator, unsigned long ntrials)
{
  unsigned long nspikes=sum(locator,ntrials);
  unsigned long i;
  double s1;
  if (nspikes>0){
    s1= 0;
    for(i=1;i<=ntrials;i++){
      if (locator[i]>0) s1+=pow((double)locator[i]/(double)nspikes,2);}
    s1*=ntrials;
    s1-=1.0;
  }

  else s1=-100;
  return(s1); 
}

double infosmall_intracell(int *locator1, int *locator2, unsigned long ntrials)
{
  unsigned long nspikes1=sum(locator1,ntrials);
  unsigned long nspikes2=sum(locator2,ntrials);
  unsigned long i;

  unsigned long nspikes12=0;
  unsigned long  nspikes1not2=0;
  unsigned long nspikesnot1not2=0;
  unsigned long  nspikes2not1=0;
  //cout<<nspikes1<<" "<<nspikes2<<" ";
  for(i=1;i<=ntrials;i++){
    if ((!locator1[i])&&(!locator2[i])) nspikesnot1not2++;
    else {
      if ((locator1[i])&&(!locator2[i])) nspikes1not2++;
      if ((locator1[i])&&(locator2[i])) nspikes12++;
      if ((!locator1[i])&&(locator2[i])) nspikes2not1++;
    }
  }
  nspikes1=nspikes12+nspikes1not2;
  nspikes2=nspikes12+nspikes2not1;
  //  cout<<"trials ="<<ntrials<<" "<<nspikes1not2+nspikes2not1+nspikesnot1not2+nspikes12<<endl;

  double s1; 
 
  if (nspikes12>0) s1=(double)nspikes12/(double)ntrials*(log((double)(nspikes12)/(double)nspikes1)+log((double)ntrials/(double)nspikes2));
  else s1=0;
  //cout<<s1<<endl; 
  if (nspikes1not2>0) s1+=(double)nspikes1not2/(double)ntrials*(log((double)nspikes1not2/(double)nspikes1)+log((double)ntrials/(double)(ntrials-nspikes2)));
  // cout<<s1<<endl; 
  if (nspikes2not1>0)  s1+=(double)nspikes2not1/(double)ntrials*(log((double)nspikes2not1/(double)nspikes2)+log((double)ntrials/(double)(ntrials-nspikes1)));
  // cout<<s1<<endl; 
  if (nspikesnot1not2>0) s1+=(double)nspikesnot1not2/(double)ntrials*(log((double)nspikesnot1not2/(double)(ntrials-nspikes2))+log((double)ntrials/(double)(ntrials-nspikes1)));
  // cout<<s1<<endl; 
  //  exit(1);
  s1/=log(2);
  //    cout<<s1<<" "<<"1 not2:"<<nspikes1not2<<" 2not1:"<<nspikes2not1<<" not1no2: "<<nspikesnot1not2<<" 12: "<<nspikes12<<endl;
  return(s1); 
}


double mean(double **v,unsigned long nv,unsigned long nh,int offset1,int offset2)
{
  double s=0;
  unsigned long i,j;
  for(i=offset1;i<(nv+offset1);i++){
    for(j=offset2;j<(nh+offset2);j++){
      s+=v[i][j];
    }
  }
  return (s/(double)(nv*nh));
}

double mean(int v[], unsigned long n)
{
  double s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=v[i];};
  return((double)s/(double)n);
}

unsigned short mean(unsigned short v[], unsigned long n)
{
  unsigned short s=0;
  unsigned long i;
  for(i=1;i<=n;i++){s+=v[i];};
  return(s/n);
}

double length(double v[],unsigned long n)
{
  double s=0;
  unsigned long i;
  for(i=1;i<=n;i++) s += v[i]*v[i];
  return(sqrt(s));
}

double dot(double *v1,double *v2,unsigned long n)
{
  double s=0;
  unsigned long i;
  for(i=1;i<=n;i++) s+=v1[i]*v2[i];
  return(s);
}

double det(double *f1,double *f2,double *v1,double *v2, unsigned long n)
{
  return( (dot(f1,v1,n)*dot(f2,v2,n)-dot(f1,v2,n)*dot(f2,v1,n))/sqrt(pow(length(v1,n)*length(v2,n),2)-pow(dot(v1,v2,n),2))/sqrt(pow(length(f1,n)*length(f2,n),2)-pow(dot(f1,f2,n),2)));
}

int factrl(int n)
{ 
  if (n==0) return 1;
  else return n*factrl(n-1);
} 


void print(double **m,unsigned long nv,unsigned long nh)     
{
  unsigned long i,j;
  for(j=1;j<=nv;j++){
    for(i=1;i<=nh;i++){cout<<"  "<<m[j][i]<<"  ";}
    cout<<endl;}
}

void print(unsigned long n,unsigned long m[])     
{
  unsigned long i;
  for(i=1;i<=n;i++){cout<<"  "<<m[i]<<"  ";}
  cout<<endl;
  cout<<endl;
}

void print(unsigned long n,double *m)     
{
  unsigned long i;
  for(i=1;i<=n;i++){cout<<m[i]<<"  ";}
  cout<<endl;
  cout<<endl;
}

double max(double ra[],unsigned long n)
{
  double a;
  unsigned long i;
  a=ra[1];
  for(i=2;i<=n;i++){if (a<ra[i]) {a=ra[i];}}
  return(a);
}

double min(double ra[],unsigned long n)
{
  double a;
  unsigned long i;
  a=ra[1];
  for(i=2;i<=n;i++){if (a>ra[i]) {a=ra[i];}}
  return(a);
}

int imax(int ra[],unsigned long n)
{
  int a;
  unsigned long i;
  a=ra[1];
  for(i=2;i<=n;i++){if (a<ra[i]) {a=ra[i];}}
  return(a);
}

int imin(int ra[],unsigned long n)
{
  int a;
  unsigned long i;
  a=ra[1];
  for(i=2;i<=n;i++){if (a>ra[i]) {a=ra[i];}}
  return(a);
}

double var(double *v,unsigned long size)
{
  double s2=0;
  double s1;
  unsigned long i;
  s1=mean(v,size);
  if (size == 1) myerror("unable to calculate variance");
  for(i=1;i<=size;i++){s2+=pow(v[i]-s1,2)/(double)(size-1);}
  return(s2);
}

double excess(double *v,unsigned long size)
{
  double s4=0;
  double s1,s2;
  unsigned long i;
  s1=mean(v,size);
  s2=var(v,size);
  for(i=1;i<=size;i++){s4+=pow(v[i]-s1,4)/(s2*s2*(double)size)-3;}
  return(s4);
}

void meanstim(double *sta,signed char *mi,unsigned long n,unsigned long fsize,  unsigned long Ntrials)
{
  unsigned long i;
  unsigned long k;

  for(i=1;i<=n;i++){ 
    sta[i]=0;
    for(k=1;k<=Ntrials;k++){
      sta[i]+=(double)mi[(k-1)*fsize+i]/(double)(Ntrials*255);}
  }
}

void meanstim(double *sta,double *mi,unsigned long n,unsigned long fsize,  unsigned long Ntrials)
{
  unsigned long i;
  unsigned long k;

  for(i=1;i<=n;i++){ 
    sta[i]=0;
    for(k=1;k<=Ntrials;k++){
      sta[i]+=mi[(k-1)*fsize+i];
    }
    sta[i]/=(double)(Ntrials);
  }
}

void meanresp(double *sta,signed char *mi,unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials)
{
  unsigned long i;
  unsigned long k;
  unsigned long sum(int *, unsigned long);
  for(i=1;i<=Nn;i++)   sta[i]=0;
  double rbar;
  Nspikes=sum(locator,Ntrials);
  rbar=(double)Nspikes/(double)Ntrials;
  for(k=1;k<=Ntrials;k++){
    for(i=1;i<=Nn;i++){
      //     sta[i]+=(double)mi[(k-1)*fsize+i]/255.0*((double)locator[k]/(double)Nspikes-(double)Nspikes/pow((double)Ntrials,2)/Nrepeat);
      //       sta[i]+=(double)mi[(k-1)*fsize+i]/255.0*((double)locator[k]-rbar)/(double)Ntrials);
       sta[i]+=(double)mi[(k-1)*fsize+i]*((double)locator[k]-rbar)/(double)Ntrials;
       //sta[i]+=(double)mi[(k-1)*fsize+i]*((double)locator[k])/(double)Ntrials;
      // sta[i]+=(double)mi[(k-1)*fsize+i]/255.0/(double)Ntrials;
    }
  }
}


void  meanresp(double *sta,double *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials)
{
  unsigned long i,k;
  double temp;
  for(i=1;i<=Nn;i++)   sta[i]=0;
  double rbar;
  rbar=(double)Nspikes/(double)Ntrials;
  for(k=1;k<=Ntrials;k++){
    temp= ((double)locator[k]-rbar)/(double)Ntrials;
    for(i=1;i<=Nn;i++){
      sta[i]+=stimuli[(k-1)*fsize+i]*temp;
    }
  }
}

void  meanresp_char(double *sta,unsigned char *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials)
{
  unsigned long i,k;
  double temp;
  for(i=1;i<=Nn;i++)   sta[i]=0;
  double rbar;
  rbar=(double)Nspikes/(double)Ntrials;
  for(k=1;k<=Ntrials;k++){
    temp= ((double)locator[k]-rbar)/(double)Ntrials;
    for(i=1;i<=Nn;i++){
      sta[i]+=((double)stimuli[(k-1)*fsize+i]-128)*temp;
    }
  }
}


void  meanresportho(double *sta2,double *sta,double *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials)
{
  unsigned long i,k;
  double temp,tempdot;
  int j;
  temp=length(sta,Nn); 
  for(i=1;i<=Nn;i++)   sta[i]/= temp;
  for(i=1;i<=Nn;i++)   sta2[i]=0;
  double *tempvector;
  tempvector=dvector(1,Nn);
  double rbar;
  Nspikes=sum(locator,Ntrials);
  rbar=(double)Nspikes/(double)Ntrials;
  for(k=1;k<=Ntrials;k++){
    temp= ((double)locator[k]-rbar)/(double)Ntrials;
    tempdot=dot(sta,stimuli+(k-1)*fsize,Nn);

    for(i=1;i<=Nn;i++){
      tempvector[i]=(stimuli[(k-1)*fsize+i]-sta[i]*tempdot);
      sta2[i]+=(stimuli[(k-1)*fsize+i]-sta[i]*tempdot)*temp;
    }
    //    cout<<"projection indi "<<dot(tempvector,sta,Nn)/length(tempvector,Nn)<<endl;
    //cin>>j;
  }
  temp=length(sta2,Nn); 
  for(i=1;i<=Nn;i++)   sta2[i]/= temp;
  temp=dot(sta2,sta,Nn); 
  cout<<"projection between 1st and 2nd sta before is "<<dot(sta,sta2,Nn)<<endl;
  for(i=1;i<=Nn;i++)   sta2[i]-= sta[i]*temp;
  temp=length(sta2,Nn); 
  for(i=1;i<=Nn;i++)   sta2[i]/= temp;
  cout<<"projection between 1st and 2nd sta is "<<dot(sta,sta2,Nn)<<endl;
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch) 
{
  free((char*) (m[nrl]+ncl-1));
  free((char*) (m+nrl-1));
}

char **cmatrix(long nrl, long nrh, long ncl, long nch) 
{
  long i, nrow, ncol;
  char **m;

  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;

  m = (char **)malloc((size_t)((nrow+1)*sizeof(char*)));
  if(!m) nrerror("allocation failure 1 in cmatrix");
  m += 1;
  m -= nrl;

  m[nrl]=(char *)malloc((size_t)((nrow*ncol+1)*sizeof(char)));
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}

