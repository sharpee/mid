#include <iostream>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#include "myutil.h"
#include "prototypes.h"
#include "ConfigFile.h"

#define NUM_THREADS 2

enum DataType {BYTE = 1, DOUBLE = 2};

int size,Nbins;
long idum=-138;

double *px,*pxt,*sx,*sxt,*pstim1,*pstim2,*pstim3,*pstim4,ispike,*xt,*df;
double *c;
int Ngrid;
int *max_ind;
int *max_ind_test;

int nlags = -1;
int LAG_SHIFT = 0;

unsigned long fsize;

double *xminN, *xmaxN, *stepN;
int *indN, *sizeN;


int main(int argc, char **argv)
{
  char errmsg[256];
  unsigned long i,k,j,n;
  
  if(argc != 4) {
    cout << "Usage: " << argv[0] << " [# of directions] [config.xml] [testrep]" << endl;
    exit(1);
  }
  cout << "Loading configuration file ..."<<endl;
  ConfigFile *cf = new ConfigFile(argv[2]);

  int spikeFileCount = 0;
  ParameterGroup *spikeGroup = cf->getParameterGroup("Spike Parameters");
  static char ** spikeFileNames = spikeGroup->getStringListValue("spike files", &spikeFileCount);  
  int testrep = atoi(argv[3]);
  int Nparts=spikeGroup->getIntValue("number of parts");
  unsigned long Ntrialsin = spikeGroup->getIntValue("number of trials");

  ParameterGroup *annealingGroup = cf->getParameterGroup("Annealing Parameters");
  int itmax = annealingGroup->getIntValue("max annealing iterations");
  double tStart =  annealingGroup->getDoubleValue("start temperature");
  double tStopTup =  annealingGroup->getDoubleValue("stop temperature");
  double tDown = annealingGroup->getDoubleValue("down temperature factor");
  double tUp = annealingGroup->getDoubleValue("up temperature factor");
  double ftol =  annealingGroup->getDoubleValue("function tolerance");
  int updateFactor = annealingGroup->getIntValue("updateFactor"); 

  ParameterGroup *movieGroup = cf->getParameterGroup("Movie Parameters"); 
  
  int movieFileCount = 0;
  static char ** movieFileNames = movieGroup->getStringListValue("movie files", &movieFileCount);
  
  if(spikeFileCount != movieFileCount) {
    myerror("There must be the same number of spike files as movie files.");
  }
  int dataType = movieGroup->getEnumValue("data type");
  //the size of receptive field in pixels
  int dimx= movieGroup->getIntValue("width");
  int dimy= movieGroup->getIntValue("height");
  
  // downsample by a factor of 4 or 2
  int cx= movieGroup->getIntValue("x downsample");
  int cy= movieGroup->getIntValue("y downsample");
  int x0 = movieGroup->getIntValue("x offset");
  int Nh = movieGroup->getIntValue("sta width");
  int y0 = movieGroup->getIntValue("y offset");
  int Nv = movieGroup->getIntValue("sta height");
  nlags = movieGroup->getIntValue("sta duration");
  LAG_SHIFT = movieGroup->getIntValue("skipped sta frames");
  int *binstep, *biniter;
  int nsteps, nstepiter;
  binstep = movieGroup->getIntegerListValue("number of bins",&nsteps);
  biniter = movieGroup->getIntegerListValue("number of iterations",&nstepiter);
  if(nsteps!=nstepiter) {
    sprintf(errmsg,"Number of Nbins values (%u) does not match number of iteration values (%u)",nsteps,nstepiter);
    myerror(errmsg);
  }
  Nbins = binstep[1];
  for(i=1;i<=nsteps;i++) cout<<"Nbins("<<i<<")="<<binstep[i]<<endl;
  
  ParameterGroup *outputGroup = cf->getParameterGroup("Output Parameters");
  static const char *prefix = outputGroup->getStringValue("prefix");
 
  cout << " done.\n";

  int Nvec = atoi(argv[1]);
  
  unsigned long Nn;
  unsigned long Ntrials;

  unsigned long Nspikes;
  unsigned long Nspikes_overall;
  double rbar;
  double *x1;
  unsigned long Ntrials_test;
  unsigned long Ntrials_test_short;
  unsigned long Ntrials_short;
  unsigned long Nspikes_test;


  ifstream inf,inf1;
  ofstream outf;
  FILE *fp;
  int *locator_test;
  double *stimuli_test;
  double *pstim1_test;
  double *stimuli; 
  double info1;
  double temp,temp1,temp2;
  double temptr,fret=0.0;
  double error_margin;
  double *sta;
  double *x;
  double fret1=0.0;

  int *repetition;
  int count,overall_count;
  int max_count;
  double *v,*dv,*vbest,*vtest;
  double std,xmin,xmax;
  double info_ave;
  double sigma_ispike;
  double ista;
  double info_test_max,info_test_ave;
  double twoinfo_test_max,twoinfo_test_ave,twoinfo_test_current;
  double achieved_value,achieved_value1;

  int iter=0;
  int nskip=0;

  char filenameLog[256]; 
  char **filenamebest;
  char **filenametest;
  char filenamePX[256];
  char buffer[256];

  cout<<"test set is "<<testrep<<endl;
  cout<<"x0="<<x0<<endl;
  cout<<"y0="<<y0<<endl;
  cout<<"Nv="<<Nv<<endl;
  cout<<"Nh="<<Nh<<endl;
  cout<<"nlags="<<nlags<<endl;
  cout<<"cx="<<cx<<endl;
  cout<<"cy="<<cy<<endl;
  cout<<"Nvec="<<Nvec<<endl;

  sizeN=ivector(1,Nvec+1);
  indN=ivector(1,Nvec);
  xminN=dvector(1,Nvec);
  xmaxN=dvector(1,Nvec);
  stepN=dvector(1,Nvec);

  sizeN[1] = 1;
  for(i = 1; i <= Nvec; i++) sizeN[i+1] = sizeN[i] * Nbins;

  fsize=(int)((Nv*Nh)/(cx*cy));
  Nn=fsize*nlags; cout<<"Nn="<<Nn<<" fsize="<<fsize<<" nlags="<<nlags<<endl;
  size=sizeN[Nvec+1];
  x1=dvector(1,binstep[nsteps]); 
  v=dvector(1,Nn*Nvec); 
  dv=dvector(1,Nn*Nvec); 
  c=dvector(1,5);
  c[3]=0;
  c[2]=-0.1; c[4]=0.1;
  c[1]=-0.2; c[5]=0.2;
  vbest=dvector(1,Nn*Nvec);
  vtest=dvector(1,Nn*Nvec);

  int Movie_length = 0;

  int *spikeFramesPerFile = ivector (1, spikeFileCount);
  spikeFramesPerFile[0] = 0;
  for(int i=1; i<=spikeFileCount; i++) {
    spikeFramesPerFile[i] = countSpikeFrames(spikeFileNames[i]);
	Movie_length += spikeFramesPerFile[i];
  }
  int * locator=ivector(1,Movie_length);
  
  for(int i=1; i<=Movie_length; i++) locator[i] = 0;
  
  Movie_length = 0;
  for(int i=1; i<=spikeFileCount; i++) {
    readinspikes(spikeFileNames[i], repetition,locator + Movie_length
	,spikeFramesPerFile[i],nlags);
	Movie_length += spikeFramesPerFile[i];
    Ntrials=Movie_length;
  } 
  
  if((Ntrialsin>=1)&&(Ntrialsin<=Movie_length)) {
    Ntrials=Ntrialsin;
  } else {
    Ntrials=Movie_length;
  }

  Ntrials_test=Ntrials/Nparts;
  Ntrials-=Ntrials_test;
  Ntrials_short=Ntrials-nlags+1-LAG_SHIFT;
  Ntrials_test_short=Ntrials_test-nlags+1-LAG_SHIFT;

  stimuli=dvector(1,Movie_length*fsize);  

  ispike=sigma_ispike=error_margin=0;
  
  if (error_margin<0.25) error_margin=0.25;
  pstim1=dvector(1,Ntrials_short*Nvec);
  pstim1_test=dvector(1,Ntrials_test_short*Nvec);
  pstim2=dvector(1,Ntrials_short);
  pstim3=dvector(1,Ntrials_short);
  pstim4=dvector(1,Ntrials_short);
  px=dvector(1,size*NUM_THREADS);
  pxt=dvector(1,size*NUM_THREADS);
  sx=dvector(1,size*Nn*NUM_THREADS);
  sxt=dvector(1,size*Nn*NUM_THREADS);   
  xt=dvector(1,Nn*Nvec);
  df=dvector(1,Nn*Nvec);
  stimuli_test=dvector(1,Ntrials_test*fsize);
  locator_test=ivector(1,Ntrials_test);


  //readin stimuli
  Movie_length = 0;
  for(int i=1; i<=spikeFileCount; i++) {
    if(dataType == BYTE)
      readinframes_to_double(movieFileNames[i],stimuli + Movie_length*fsize,
	    spikeFramesPerFile[i],dimx,dimy,x0,y0,Nh,Nv,cx,cy);
  else if(dataType == DOUBLE) {
    readinframes_as_double(movieFileNames[i],stimuli + Movie_length*fsize,
	  spikeFramesPerFile[i],dimx,dimy,x0,y0,Nh,Nv,cx,cy);
  } else {
	myerror("Data type of movie file is not recognized.\n");
  }
   Movie_length += spikeFramesPerFile[i];
 }

  cout<<"overall number of spikes ="<<sum(locator,Ntrials+Ntrials_test)<<endl;
  

  Nspikes_overall=sum(locator,Ntrials+Ntrials_test);
  rbar=(double)Nspikes_overall/(double)(Ntrials+Ntrials_test)/imax(locator,Ntrials_short+Ntrials_test_short);
  cout<<Nspikes_overall<<endl;
  temp= -log(rbar)/log(2);
  for(i=1;i<=Ntrials+Ntrials_test;i++){if (locator[i]>0) temp+=(double)locator[i]*log((double)locator[i])/log(2)/(double)Nspikes_overall;}
  cout<<"information per spike (empirically) ="<<temp<<endl;  
  cout<<"from now on all information is given in units of information per spike"<<endl;

  if (ispike==0) ispike=temp;
  ispike*=log(2);
  
  cout<<Ntrials+Ntrials_test<<" "<<Movie_length<<endl;

  for(i=1;i<=Ntrials_test;i++) {
    locator_test[i]=locator[i+(testrep-1)*Ntrials_test];
    for(k=1;k<=fsize;k++)   
	  stimuli_test[(i-1)*fsize+k]=stimuli[(i+(testrep-1)*Ntrials_test-1)*fsize+k];
	 
  }
  for(i=1;i<=Ntrials_test*(Nparts-testrep);i++){
    locator[i+(testrep-1)*Ntrials_test]=locator[i+testrep*Ntrials_test]; 
    for(k=1;k<=fsize;k++)  
	  stimuli[(i+(testrep-1)*Ntrials_test-1)*fsize+k]=stimuli[(i+testrep*Ntrials_test-1)*fsize+k];
  }
  if (testrep !=Nparts){
    for(i=1;i<=Ntrials_test;i++)  {
      locator[i+(Nparts-1)*Ntrials_test]=locator_test[i]; 
      for(k=1;k<=fsize;k++)  
	    stimuli[(i+(Nparts-1)*Ntrials_test-1)*fsize+k]=stimuli_test[(i-1)*fsize+k]; 
    }
  }
  
  for(i=1;i<=nlags-1;i++) {
    if ((testrep!=Nparts)&&(testrep!=1)){
      locator[(testrep-1)*Ntrials_test+i-nlags+1]=0;
    }
    locator[(Nparts-1)*Ntrials_test+i-nlags+1]=0;
  }
  cout<<"number of trials -(nlags-1+LAG_SHIFT)="<<Ntrials_short<<" Ntrials="<<Ntrials<<endl;
  Nspikes=sum(locator+LAG_SHIFT,Ntrials_short);
  Nspikes_test=sum(locator_test+LAG_SHIFT,Ntrials_test_short);
  cout<<"number of spikes ="<<Nspikes<<endl;
  cout<<"number of spikes in the test case = "<<Nspikes_test<<endl; 
  
  rbar=(double)Nspikes/(double)(Ntrials_short)/imax(locator+LAG_SHIFT,Ntrials_short);

  filenamebest = cmatrix(1,Nvec,0,255);
  filenametest = cmatrix(1,Nvec,0,255);

  sprintf(filenameLog,"%s-ND-n%u-p%u.log",prefix,Nvec,testrep); 
  cout<<filenameLog<<endl;
  for(i = 1; i <= Nvec; i++) {
    sprintf(filenamebest[i],"%s-ND-n%u-v%u-p%u.bst",prefix,Nvec,i,testrep); 
    cout<<filenamebest[i]<<endl;
    sprintf(filenametest[i],"%s-ND-n%u-v%u-p%u.dat",prefix,Nvec,i,testrep); 
    cout<<filenametest[i]<<endl;
  }
  sprintf(filenamePX,"%s-ND-n%u-p%u.pxt",prefix,Nvec,testrep); 
  cout<<filenamePX<<endl;

  if(Nvec > 1) {
    for(i = 1; i < Nvec; i++) {
      //Copy the vector files from the previous dimension
      if (Nvec == 2) {
	  sprintf(buffer,"cp %s-1D-n%u-v%u-p%u.bst %s",prefix,Nvec-1,i,testrep,filenamebest[i]);
	  cout << buffer << endl;
	  system(buffer);
	  sprintf(buffer,"cp %s-1D-n%u-v%u-p%u.dat %s",prefix,Nvec-1,i,testrep,filenametest[i]);
	  cout << buffer << endl;
	  system(buffer);
      }
	else{
	  sprintf(buffer,"cp %s-ND-n%u-v%u-p%u.bst %s",prefix,Nvec-1,i,testrep,filenamebest[i]);
	  cout << buffer << endl;
	  system(buffer);
	  sprintf(buffer,"cp %s-ND-n%u-v%u-p%u.dat %s",prefix,Nvec-1,i,testrep,filenametest[i]);
	  cout << buffer << endl;
	  system(buffer);
	}

      //Load test vectors
      fp=fopen(filenametest[i],"rb");   
      cout<<"reading from file "<<filenametest[i]<<endl;
      if ( fp == NULL)    myerror(filenametest[i]);
      fread(v+1+(i-1)*Nn,sizeof(double),Nn,fp);
      fclose(fp);
      
      //Orthogonalize vector
      for(j = 1; j < i; j++) {
	temp=dot(v+(j-1)*Nn,v+(i-1)*Nn,Nn);
	for(k = 1; k <= Nn; k++) v[k+Nn*(i-1)] -= temp*v[k+(j-1)*Nn];
      }
    
      //Normalize vector
      temp=length(v+(i-1)*Nn,Nn);
      for(j=1;j<=Nn;j++) v[j+(i-1)*Nn]/= temp;
    }

    //Create new vector from stimulus
    for(i=1;i<=Nn;i++) v[i+(Nvec-1)*Nn] += stimuli[i+fsize*(100*Nvec+(testrep-1)*Ntrials_test)%(Ntrials_short*fsize-Nn)];

    //Orthogonalize vector
    for(i = 0; i < Nvec - 1; i++) {
      temp=dot(v+i*Nn,v+(Nvec-1)*Nn,Nn);
      for(j = 1; j <= Nn; j++) v[j+(Nvec-1)*Nn] -= temp*v[j+i*Nn];
    }

    //Normalize vector
    temp=length(v+(Nvec-1)*Nn,Nn);
    for(j=1; j<=Nn; j++) v[j+(Nvec-1)*Nn] /= temp;
  } else {
    //Create first vector from STA
    meanresp(v,stimuli,Nn,fsize,locator+LAG_SHIFT,Nspikes,Ntrials_short);
    temp=length(v,Nn);
    for(i=1; i<=Nn; i++) v[i] /= temp;
  }

  //Create projections for information calculation
  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,v+i*Nn,stimuli,Nn,Ntrials_short);

  //Calculate starting information for vN
  achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);
  cout<<"information about "<<Nvec<<"  vectors "<<-achieved_value<<endl;
  temptr=0.001;//Tstart; 
  
  cout<<"information in the starting vector = "<<-achieved_value<<endl;
  
  count=0;
  overall_count=1;
  max_count=biniter[1];
  //no need for modifications neyond this point
  // the first MID v[1:Nn], the second MId is v[1+Nn:2*Nn];  
  for(i=1;i<=Nvec*Nn;i++) vbest[i]=v[i];

  fp = fopen(filenametest[Nvec],"wb");
  if(fp == NULL) cout<<"could not open "<<filenametest[1]<<endl;
  fwrite(v+1+(Nvec-1)*Nn,sizeof(double),Nn,fp);
  fclose(fp);
  
  fp = fopen(filenamebest[Nvec],"wb");
  if(fp == NULL) cout<<"could not open "<<filenamebest[1]<<endl;
  fwrite(vbest+1+(Nvec-1)*Nn,sizeof(double),Nn,fp);
  fclose(fp);

  double infomax;
  int *infoind;
  infoind=ivector(1,Nvec);

  for(k=1;k<=nsteps;k++) {
    if(k > 1) {
      //Free vectors from previous step
      free_dvector(sxt,1,size*Nn*NUM_THREADS);
      free_dvector(sx,1,size*Nn*NUM_THREADS);
      free_dvector(pxt,1,size*NUM_THREADS);
      free_dvector(px,1,size*NUM_THREADS);

      Nbins=binstep[k];
      sizeN[1] = 1;
      for(i = 1; i <= Nvec; i++) sizeN[i+1] = sizeN[i] * Nbins;
      size = sizeN[Nvec+1];

      for(j=1;j<=Nvec;j++) {
	fp=fopen(filenametest[j],"rb");
	if(fp==NULL) myerror(filenametest[j]);
	fread(vtest+1+(j-1)*Nn,sizeof(double),Nn,fp);
	fclose(fp);
      }

      sxt=dvector(1,size*Nn*NUM_THREADS);
      sx=dvector(1,size*Nn*NUM_THREADS);
      pxt=dvector(1,size*NUM_THREADS);
      px=dvector(1,size*NUM_THREADS);

      max_count+=biniter[k];
      twoinfo_test_max=twoinfo_test_current=-Nneginform_double_test(Nvec,vtest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
      cout<<"information in the starting vector along the test set= "<<twoinfo_test_max<<endl;
    } else {

      twoinfo_test_max=twoinfo_test_current=-Nneginform_double_test(Nvec,vbest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
      cout<<"information in the starting vector along the test set= "<<twoinfo_test_max<<endl;
    }
    while (overall_count<= max_count)  { 
      if(overall_count%100==0) {
	for(j=0;j<Nvec-1;j++) {
	  //Exchange v_j and v_Nvec
	  for(i = 1; i <= Nn; i++) {
	    temp=v[i+j*Nn];
	    v[i+j*Nn]=v[i+(Nvec-1)*Nn];
	    v[i+(Nvec-1)*Nn]=temp;
	    temp=vbest[i+j*Nn];
	    vbest[i+j*Nn]=vbest[i+(Nvec-1)*Nn];
	    vbest[i+(Nvec-1)*Nn]=temp;
	  }

	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,vbest+i*Nn,stimuli,Nn,Ntrials_short);
	  achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);
	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,v+i*Nn,stimuli,Nn,Ntrials_short);
	  fret=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);
	  
	  // The main optimization subroutine, the last parameters are the names of the function to calculate information along v+Nn given v, and the name of the function to calcaulte 2D derivative for v+Nn given v

	  Ntrysa_double(Nvec,v,vbest,dv,Nn,2,stimuli,locator+LAG_SHIFT,Nspikes,Ntrials_short,ftol,&temptr,itmax,&count,updateFactor,&achieved_value,filenameLog,tStart,tDown,tUp,tStopTup,&iter,&fret,&info_ave,Nneginform_double,Nderivative_double);
	  if (fret >= 0){
	    cout<<"error in trysa_double with twonegfixed optimization "<<fret;
	    cout<<" "<<v[1]<<" "<<v[Nn]<<" "<<vbest[1]<<" "<<vbest[Nn]<<" "<<dv[1]<<" "<<dv[Nn]<<endl;
	    exit(1);
	  }
	  //Exchange v_j and v_Nvec
	  for(i = 1; i <= Nn; i++) {
	    temp=v[i+j*Nn];
	    v[i+j*Nn]=v[i+(Nvec-1)*Nn];
	    v[i+(Nvec-1)*Nn]=temp;
	    temp=vbest[i+j*Nn];
	    vbest[i+j*Nn]=vbest[i+(Nvec-1)*Nn];
	    vbest[i+(Nvec-1)*Nn]=temp;
	  }
	  twoinfo_test_current=-Nneginform_double_test(Nvec,vbest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
	  if (twoinfo_test_current>twoinfo_test_max) {
	    twoinfo_test_max=twoinfo_test_current;

	    for(i=1; i<=Nvec; i++) infoind[i]=0;
	    for(i=1; i<=Nvec; i++) {
	      infomax=0.0;
	      for(j=1; j<=Nvec; j++) {
		temp=0.0;
		for(n=1; n<i; n++) {
		  if(j==infoind[n]) temp = -1.0;
		}
		if(temp >= 0.0) {
		  for(n=1; n<=Nn; n++) vtest[n+(i-1)*Nn]=vbest[n+(j-1)*Nn];
		  temp=-Nneginform_double_test(i,vtest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
		  if(infomax < temp) {
		    infoind[i]=j;
		    infomax=temp;
		  }
		}
	      }
	      for(n=1; n<=Nn; n++) vtest[n+(i-1)*Nn]=vbest[n+(infoind[i]-1)*Nn];
	    }

	    for(i = 1; i <= Nvec; i++) {
	      fp=fopen(filenametest[i],"wb");   
	      if (fp == NULL) cout<<"could not open "<<filenametest[i]<<endl; 
	      if (fp!=NULL) {
		fwrite(vtest+1+(i-1)*Nn,sizeof(double),Nn,fp);
		fclose(fp);
	      }
	    }
	  }
	}	  
	cout<<overall_count<<" "<<twoinfo_test_current<<" "<<twoinfo_test_max<<" "<<v[1+(Nvec-1)*Nn]<<endl;
	  
	overall_count++;
	  
	outf.open(filenameLog,ios::out|ios::app);
	if (! outf.is_open()) cout<<"could not open "<<filenameLog<<endl;
	if ( outf.is_open()){ 
	  outf<<"#old achieved value="<<-achieved_value<<" current (two)info_test="<<twoinfo_test_current<<" ";
	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,vbest+i*Nn,stimuli,Nn,Ntrials_short);
	  achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short); 
	  outf<<"#new achieved value="<<-achieved_value<<" max (two)info_test="<<twoinfo_test_max<<endl;
	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,v+i*Nn,stimuli,Nn,Ntrials_short);
	  outf<<(-Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short))<<" ";
	  outf<<temptr<<" "<<count<<" "<<overall_count<<" "<<iter<<endl;
	  outf.close();
	}
	for(i = 1; i <= Nvec; i++) {
	  fp=fopen(filenamebest[i],"wb");   
	  if (fp == NULL) cout<<"could not open "<<filenamebest[i]<<endl; 
	  if (fp!=NULL) {
	    fwrite(vbest+1+(i-1)*Nn,sizeof(double),Nn,fp);
	    fclose(fp);
	  }
	}
	if  (temptr < tStopTup ) {temptr=tStopTup;} 
      } else {
	for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,vbest+i*Nn,stimuli,Nn,Ntrials_short);
	achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);
	for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,v+i*Nn,stimuli,Nn,Ntrials_short);
	fret=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);

	// The main optimization subroutine, the last parameters are the names of the function to calculate information along v+Nn given v, and the name of the function to calcaulte 2D derivative for v+Nn given v

	Ntrysa_double(Nvec,v,vbest,dv,Nn,2,stimuli,locator+LAG_SHIFT,Nspikes,Ntrials_short,ftol,&temptr,itmax,&count,updateFactor,&achieved_value,filenameLog,tStart,tDown,tUp,tStopTup,&iter,&fret,&info_ave,Nneginform_double,Nderivative_double);
	if (fret >= 0){
	  cout<<"error in trysa_double with twonegfixed optimization "<<fret;
	  cout<<" "<<v[1]<<" "<<v[Nn]<<" "<<vbest[1]<<" "<<vbest[Nn]<<" "<<dv[1]<<" "<<dv[Nn]<<endl;
	  exit(1);
	}
	twoinfo_test_current=-Nneginform_double_test(Nvec,vbest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
	if (twoinfo_test_current>twoinfo_test_max) {
	  twoinfo_test_max=twoinfo_test_current;

	  for(i=1; i<=Nvec; i++) infoind[i]=0;
	  for(i=1; i<=Nvec; i++) {
	    infomax=0.0;
	    for(j=1; j<=Nvec; j++) {
	      temp=0.0;
	      for(n=1; n<i; n++) {
		if(j==infoind[n]) temp = -1.0;
	      }
	      if(temp >= 0.0) {
		for(n=1; n<=Nn; n++) vtest[n+(i-1)*Nn]=vbest[n+(j-1)*Nn];
		temp=-Nneginform_double_test(i,vtest,stimuli_test,Nn,locator_test+LAG_SHIFT,Nspikes_test,Ntrials_test_short,pstim1_test);
		if(infomax < temp) {
		  infoind[i]=j;
		  infomax=temp;
		}
	      }
	    }
	    for(n=1; n<=Nn; n++) vtest[n+(i-1)*Nn]=vbest[n+(infoind[i]-1)*Nn];
	  }

	  for(i = 1; i <= Nvec; i++) {
	    fp=fopen(filenametest[i],"wb");   
	    if (fp == NULL) cout<<"could not open "<<filenametest[i]<<endl; 
	    if (fp!=NULL) {
	      fwrite(vtest+1+(i-1)*Nn,sizeof(double),Nn,fp);
	      fclose(fp);
	    }
	  }
	}
	
	cout<<overall_count<<" "<<twoinfo_test_current<<" "<<twoinfo_test_max<<" "<<v[1+(Nvec-1)*Nn]<<endl;
	
	overall_count++;
        
	outf.open(filenameLog,ios::out|ios::app);
	if (! outf.is_open()) cout<<"could not open "<<filenameLog<<endl;
	if ( outf.is_open()){ 
	  outf<<"#old achieved value="<<-achieved_value<<" current (two)info_test="<<twoinfo_test_current<<" ";
	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,vbest+i*Nn,stimuli,Nn,Ntrials_short);
	  achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short); 
	  outf<<"#new achieved value="<<-achieved_value<<" max (two)info_test="<<twoinfo_test_max<<endl;
	  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,v+i*Nn,stimuli,Nn,Ntrials_short);
	  outf<<(-Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short))<<" ";
	  outf<<temptr<<" "<<count<<" "<<overall_count<<" "<<iter<<endl;
	  outf.close();
	}
	for(i = 1; i <= Nvec; i++) {
	  fp=fopen(filenamebest[i],"wb");   
	  if (fp == NULL) cout<<"could not open "<<filenamebest[i]<<endl; 
	  if (fp!=NULL) {
	    fwrite(vbest+1+(i-1)*Nn,sizeof(double),Nn,fp);
	    fclose(fp);
	  }
	}
	if  (temptr < tStopTup ) {temptr=tStopTup;} 
      }
    }
  }

  //writing input-output function
  for(i = 1; i <= Nvec; i++) {
    fp=fopen(filenamebest[i],"wb");   
    if (fp == NULL) cout<<"could not open "<<filenamebest[i]<<endl; 
    if (fp!=NULL) {
      fwrite(vbest+1+(i-1)*Nn,sizeof(double),Nn,fp);
      fclose(fp);
    }
  }

  fp=fopen(filenamePX,"wb");
  if ( fp ==NULL ) myerror(filenamePX);
  for(i = 0; i < Nvec; i++) Proj1d_double(pstim1+i*Ntrials_short,vbest+i*Nn,stimuli,Nn,Ntrials_short);
  achieved_value=Nneginform_double(Nvec,0.0,locator+LAG_SHIFT,Nspikes,Ntrials_short);
  for(i = 1; i <= Nvec; i++) {
    std=sqrt(var(pstim1+(i-1)*Ntrials_short,Ntrials_short));
    xminN[i] = min(pstim1+(i-1)*Ntrials_short,Ntrials_short);
    xmaxN[i] = max(pstim1+(i-1)*Ntrials_short,Ntrials_short);
    xminN[i]/=std;
    xmaxN[i]/=std;
    for(j=1;j<=Nbins;j++) x1[j]=xminN[i]+(double)(j-1)/(double)(Nbins)*(xmaxN[i]-xminN[i]);
    fwrite(x1+1,sizeof(double),Nbins,fp);
  }
  for(i=1;i<=size;i++) pxt[i]*=rbar;
  fwrite(px+1,sizeof(double),size,fp);
  fwrite(pxt+1,sizeof(double),size,fp);
  fclose(fp);
  
  cout<<"successfully wrote into the file for testrep: "<<testrep<<endl;

  free_ivector(infoind,1,Nvec);

  free_cmatrix(filenametest,1,Nvec,0,255);
  free_cmatrix(filenamebest,1,Nvec,0,255);

  free_ivector(sizeN,1,Nvec+1);
  free_ivector(indN,1,Nvec);
  free_dvector(xminN,1,Nvec);
  free_dvector(xmaxN,1,Nvec);
  free_dvector(stepN,1,Nvec);
  free_dvector(x1,1,Nbins);

  free_ivector(locator_test,1,Ntrials_test);   
  free_dvector(stimuli_test,1,Ntrials_test*fsize);
  free_dvector(df,1,Nn); 
  free_dvector(xt,1,Nn);
  free_dvector(sxt,1,size*Nn*NUM_THREADS); 
  free_dvector(sx,1,size*Nn*NUM_THREADS); 
  free_dvector(pxt,1,size*NUM_THREADS); 
  free_dvector(px,1,size*NUM_THREADS);
  free_dvector(pstim4,1,Ntrials_short);
  free_dvector(pstim3,1,Ntrials_short);
  free_dvector(pstim2,1,Ntrials_short);
  free_dvector(pstim1_test,1,Ntrials_test_short*Nvec);
  free_dvector(pstim1,1,Ntrials_short*Nvec);

  free_ivector(locator,1,Movie_length);
  free_dvector(stimuli,1,Movie_length*fsize);

  free_dvector(vbest,1,Nn*Nvec);
  free_dvector(vtest,1,Nn*Nvec);
  free_dvector(c,1,5);
  free_dvector(dv,1,Nn); 
  free_dvector(v,1,Nn*Nvec);

  delete cf;

  cout<<"Fin"<<endl;

  return 0; 
}

