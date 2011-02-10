/*
 * Prototypes for functions defined in myutil.c
 */

void myerror(const char *error_text);

void erroropen(char error_text[]);

void putzeros(char * buffer);

unsigned long sum(int v[], unsigned long n);

unsigned long csum(unsigned char v[], unsigned long n);

double dsum(double v[], unsigned long n);

void rlft2(double **data, double *speq, unsigned long nn1,
	       unsigned long nn2, int isign);

double mean(double v[],unsigned long n);

double infosmall(int *locator, unsigned long ntrials);

double varsmall(int *locator, unsigned long ntrials);

double infosmall_intracell(int *locator1, int *locator2, unsigned long ntrials);

double mean(double **v,unsigned long nv,unsigned long nh,int offset1,int offset2);

double mean(int v[], unsigned long n);

unsigned short mean(unsigned short v[], unsigned long n);

double length(double v[],unsigned long n);

double dot(double *v1,double *v2,unsigned long n);

double det(double *f1,double *f2,double *v1,double *v2, unsigned long n);

int factrl(int n);

void print(double **m,unsigned long nv,unsigned long nh);

void print(unsigned long n,unsigned long m[]);

void print(unsigned long n,double *m);

double max(double ra[],unsigned long n);

double min(double ra[],unsigned long n);

int imax(int ra[],unsigned long n);

int imin(int ra[],unsigned long n);

double var(double *v,unsigned long size);

double excess(double *v,unsigned long size);

void meanstim(double *sta,signed char *mi,unsigned long n,unsigned long fsize,  unsigned long Ntrials);

void meanstim(double *sta,double *mi,unsigned long n,unsigned long fsize,  unsigned long Ntrials);

void meanresp(double *sta,signed char *mi,unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials);

void  meanresp(double *sta,double *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials);

void  meanresp_char(double *sta,unsigned char *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials);

void  meanresportho(double *sta2,double *sta,double *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials);

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);

char **cmatrix(long nrl, long nrh, long ncl, long nch);

