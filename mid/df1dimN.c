extern double (*nrdfundN)(int Nvec, double a, double [], double [], double [], double *stimuli, unsigned long n, int *locator, unsigned long Nspikes, unsigned long Ntrials);
extern double *xt,*df;


double Ndf1dimN(int Nvec, double x, double *p, double *xi, unsigned long n, double *stimuli, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  unsigned long j;
  double df1 = 0.0;

  for(j = 1; j <= n; j++) xt[j] = p[j + (Nvec-1)*n] + x * xi[j];
  (*nrdfundN)(Nvec, x, p, xt, df, stimuli, n, locator, Nspikes, Ntrials);
  for(j  = 1; j <= n; j++) df1 += df[j]*xi[j];
  return df1;
}

