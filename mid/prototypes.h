#include <stdio.h>

/*
 * Function prototypes
 */

void Ntrysa_double(int, double [], double [], double [], unsigned long, int, double *, int *, unsigned long, unsigned long, double, double *, int, int *, double, double *, char *, double, double, double, double, int *, double *, double *, double (*)(int, double, int *, unsigned long, unsigned long), void (*)(int, double, double *, double *, double *, double *, unsigned long, int *, unsigned long, unsigned long));

double Proj1d_double(double *, double *, double *, unsigned long, unsigned long);

double chooseN_double(int, double *, double *, double *, unsigned long, int *,unsigned long, unsigned long,double (*)(int, double, int *,unsigned long, unsigned long));

double NdbrentsaN_double(int, double, double, double, double *, double (*)(int, double, int *, unsigned long, unsigned long),double (*)(int, double, double *, double *, unsigned long, double *,int *, unsigned long, unsigned long), double, double *, double *, double *, unsigned long, double *,int *, unsigned long, unsigned long);

void Nderivative_double(int, double, double *, double *, double *, double *, unsigned long, int *, unsigned long, unsigned long);

void Nderivative1d_double(int, double, double *, double *, double *, double *, unsigned long, int *, unsigned long, unsigned long);

double Ndf1dimN(int, double, double *, double *, unsigned long, double *, int *, unsigned long, unsigned long);

void NdlinminsaN_double(int, double [], double [], unsigned long, int, double *, int *, unsigned long, unsigned long, double *, double *, double (*)(int, double, int *, unsigned long, unsigned long), void (*)(int, double, double *, double *, double *, double *, unsigned long, int *, unsigned long, unsigned long));

void NmnbrakN(int, double *, double *, double *, double *, double *, double *,double (*)(int, double, int *,unsigned long, unsigned long), double *, double *, unsigned long,double *, int *, unsigned long, unsigned long);

double Nf1dimN(int, double, int *, unsigned long, unsigned long);

double Nneginform1d_double(int, double, int *, unsigned long, unsigned long);

double Nneginform_double_test(int, double *,double *, unsigned long, int *, unsigned long, unsigned long, double *);

double Nneginform_double(int, double, int *, unsigned long, unsigned long);

double neginform_double(double *,double *, unsigned long, int *, unsigned long, unsigned long);

double twoneginform_double(double *,double *,unsigned long, int *,unsigned long, unsigned long);

void runChecks(FILE *, unsigned long, int, int, int, int, int, int, int, int, int);

void readinframes_to_double(char *,double *, unsigned long,int, int, int, int, int, int, int, int);

void readinframes_as_double(char *,double *, unsigned long, int, int, int, int, int, int, int, int);

int countSpikeFrames(char *);

void readinspikes(char *, int *, int *, int, int);

double ran2(long *);

