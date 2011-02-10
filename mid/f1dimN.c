extern double (*nrfuncdN)(int NVec, double a, int *locator, unsigned long Nspikes, unsigned long Ntrials);

double Nf1dimN(int Nvec, double x, int *locator, unsigned long Nspikes, unsigned long Ntrials) 
{
  double f;
  f = (*nrfuncdN)(Nvec,x,locator,Nspikes,Ntrials);
  return f;
}

